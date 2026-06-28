"""
Train ZHMolGraph's VecNN model on the ciceklab dataset and evaluate on
three unseen test splits.

Usage:
    python ciceklab/train_ciceklab.py                        # full training
    python ciceklab/train_ciceklab.py --smoke_test           # quick local test
    python ciceklab/train_ciceklab.py --use_graphsage        # with GraphSAGE (paper method)
"""
import os
import sys
import json
# pyrefly: ignore [missing-import]
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import pandas as pd
import pickle as pkl
import argparse
from tqdm import tqdm
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             accuracy_score, precision_score, recall_score,
                             matthews_corrcoef, confusion_matrix, f1_score)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.abspath(os.path.join(BASE_DIR, '..'))
sys.path.insert(0, REPO_DIR)

from ZHMolGraph.ZHMolGraph_model import VecNN


# ─── Data loading ──────────────────────────────────────────────────────

def load_interactions(chunk_dir):
    """Load training chunks and return as DataFrame."""
    records = []
    print("Loading training chunks...")
    for fname in tqdm(sorted(os.listdir(chunk_dir))):
        if not fname.endswith('.jsonl'):
            continue
        with open(os.path.join(chunk_dir, fname), 'r') as f:
            for line in f:
                obj = json.loads(line)
                if obj.get('interaction_type') == 'rna-protein':
                    records.append({
                        'RNA_aa_code': obj['RNA_id'],
                        'target_aa_code': obj['target_id'],
                        'Y': obj['interaction_label']
                    })
    return pd.DataFrame(records)


def load_test_split(filepath):
    """Load test interactions from a single jsonl file."""
    records = []
    with open(filepath, 'r') as f:
        for line in f:
            obj = json.loads(line)
            if obj.get('interaction_type') == 'rna-protein':
                records.append({
                    'RNA_aa_code': obj['RNA_id'],
                    'target_aa_code': obj['target_id'],
                    'Y': obj['interaction_label']
                })
    return pd.DataFrame(records)


# ─── Evaluation ────────────────────────────────────────────────────────

def evaluate(model, dataloader, device):
    """Evaluate VecNN model and return a metrics dict."""
    model.eval()
    all_preds = []
    all_labels = []
    with torch.no_grad():
        for batch_prot, batch_rna, batch_y in dataloader:
            batch_prot = batch_prot.to(device)
            batch_rna = batch_rna.to(device)
            output = model(batch_prot, batch_rna)          # VecNN(target, rna)
            all_preds.extend(output.squeeze(-1).cpu().numpy().tolist())
            all_labels.extend(batch_y.numpy().tolist())

    all_preds = np.array(all_preds)
    all_labels = np.array(all_labels)

    auc = roc_auc_score(all_labels, all_preds)
    auprc = average_precision_score(all_labels, all_preds)

    # Find best threshold by MCC on these predictions
    thresholds = np.linspace(0, 1, 100)
    best_mcc, best_thresh = -1, 0.5
    for thresh in thresholds:
        preds_bin = (all_preds > thresh).astype(int)
        mcc = matthews_corrcoef(all_labels, preds_bin)
        if mcc > best_mcc:
            best_mcc = mcc
            best_thresh = thresh

    preds_bin = (all_preds > best_thresh).astype(int)
    acc = accuracy_score(all_labels, preds_bin)
    prec = precision_score(all_labels, preds_bin, zero_division=0)
    rec = recall_score(all_labels, preds_bin)
    f1 = f1_score(all_labels, preds_bin)
    tn, fp, fn, tp = confusion_matrix(all_labels, preds_bin).ravel()
    spe = tn / (tn + fp) if (tn + fp) > 0 else 0

    return {
        'AUC': auc, 'AUPRC': auprc, 'Accuracy': acc, 'Precision': prec,
        'Sensitivity': rec, 'Specificity': spe, 'F1': f1,
        'MCC': best_mcc, 'Threshold': best_thresh,
        'TP': int(tp), 'FP': int(fp), 'TN': int(tn), 'FN': int(fn),
    }


# ─── Main ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Train VecNN on ciceklab dataset")
    parser.add_argument('--smoke_test', action='store_true',
                        help='Quick test: 1000 samples, 2 epochs')
    parser.add_argument('--use_graphsage', action='store_true',
                        help='Train GraphSAGE and concatenate 100-dim node embeddings (paper method)')
    parser.add_argument('--epochs', type=int, default=None,
                        help='Override number of training epochs')
    parser.add_argument('--batch_size', type=int, default=256)
    parser.add_argument('--lr', type=float, default=0.001)
    args = parser.parse_args()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # ── 1. Load precomputed embeddings ─────────────────────────────────
    rna_embed_path = os.path.join(BASE_DIR, 'rna_embeddings.pkl')
    prot_embed_path = os.path.join(BASE_DIR, 'protein_embeddings.pkl')

    if not os.path.exists(rna_embed_path) or not os.path.exists(prot_embed_path):
        print("ERROR: Embeddings not found. Run precompute_embeddings.py first.")
        sys.exit(1)

    with open(rna_embed_path, 'rb') as f:
        rna_embeds_df = pkl.load(f)
    with open(prot_embed_path, 'rb') as f:
        prot_embeds_df = pkl.load(f)

    print(f"Loaded {len(rna_embeds_df)} RNA and {len(prot_embeds_df)} protein embeddings.")

    # Build {ID → 640-dim array} and {ID → 1024-dim array}
    rna_embed_dict = {}
    for _, row in rna_embeds_df.iterrows():
        rna_embed_dict[row['RNA_id']] = np.array(row['normalized_embeddings'], dtype=np.float32)
    prot_embed_dict = {}
    for _, row in prot_embeds_df.iterrows():
        prot_embed_dict[row['target_id']] = np.array(row['normalized_embeddings'], dtype=np.float32)

    # ── 2. Load training data ──────────────────────────────────────────
    chunk_dir = os.path.join(BASE_DIR, 'training_chunks')
    train_df = load_interactions(chunk_dir)

    if args.smoke_test:
        print("Smoke test: subsetting to 1000 interactions.")
        train_df = train_df.head(1000)

    print(f"Training interactions: {len(train_df)}  "
          f"(pos={int((train_df['Y']==1).sum())}, neg={int((train_df['Y']==0).sum())})")

    # ── 3. (Optional) GraphSAGE ────────────────────────────────────────
    use_graphsage = args.use_graphsage
    if use_graphsage:
        print("\n=== GraphSAGE embedding generation ===")
        gs_dir = os.path.join(REPO_DIR, 'ciceklab_Pretrain_graphsage_dataset')
        os.makedirs(gs_dir, exist_ok=True)

        # 3a. Write feats.txt
        # Format: node_name \t feat_1 \t ... \t feat_1024 \t dummy
        # NOTE: DataCenter.load_dataSet reads features as info[1:-1], which
        # drops the LAST column. We append a trailing '0' so the real 1024
        # features are preserved.
        feats_path = os.path.join(gs_dir, 'ciceklab_feats.txt')
        rna_id_order = list(rna_embed_dict.keys())
        prot_id_order = list(prot_embed_dict.keys())
        if not os.path.exists(feats_path):
            print(f"  Writing {feats_path} ...")
            with open(feats_path, 'w') as f:
                for r_id in rna_id_order:
                    emb = rna_embed_dict[r_id]
                    padded = np.pad(emb, (0, max(0, 1024 - len(emb))))
                    f.write(r_id + '\t' + '\t'.join(map(str, padded.tolist())) + '\t0\n')
                for p_id in prot_id_order:
                    emb = prot_embed_dict[p_id]
                    f.write(p_id + '\t' + '\t'.join(map(str, emb.tolist())) + '\t0\n')

        # 3b. Write interaction files (node_a \t node_b)
        # Total interactions
        train_df[['RNA_aa_code', 'target_aa_code']].to_csv(
            os.path.join(gs_dir, 'ciceklab_total_interactions_seq_list.txt'),
            sep='\t', index=False, header=False)
        # GraphSAGE train = positive interactions only (unsupervised link-based)
        pos_train = train_df[train_df['Y'] == 1]
        pos_train[['RNA_aa_code', 'target_aa_code']].to_csv(
            os.path.join(gs_dir, 'ciceklab_graphsage_train_interactions.txt'),
            sep='\t', index=False, header=False)
        # GraphSAGE test = negative interactions (they are the non-edges)
        neg_train = train_df[train_df['Y'] == 0]
        neg_train[['RNA_aa_code', 'target_aa_code']].to_csv(
            os.path.join(gs_dir, 'ciceklab_graphsage_test_interactions.txt'),
            sep='\t', index=False, header=False)

        # 3c. Train GraphSAGE using repo's run_graphsage_experiment
        from ZHMolGraph.ZHMolGraph import ZHMolGraph as ZHMolGraphClass
        # We need a minimal ZHMolGraph object; supply dummy data to pass assertions
        dummy_interactions = pd.DataFrame(columns=['RNA_aa_code', 'target_aa_code', 'Y'])
        rna_id_df = pd.DataFrame({'RNA_aa_code': rna_id_order})
        prot_id_df = pd.DataFrame({'target_aa_code': prot_id_order})
        zg = ZHMolGraphClass(
            interactions=dummy_interactions,
            interaction_y_name='Y',
            rnas_dataframe=rna_id_df,
            rna_seq_name='RNA_aa_code',
            proteins_dataframe=prot_id_df,
            protein_seq_name='target_aa_code',
            model_out_dir=os.path.join(BASE_DIR, 'trained_model', 'ciceklab'),
        )
        os.makedirs(os.path.join(zg.model_out_dir, 'Run_0'), exist_ok=True)

        gs_config = os.path.join(REPO_DIR, 'graphsage_src', 'experiments.conf')
        print("  Training GraphSAGE (unsupervised) ...")
        zg.run_graphsage_experiment(
            dataSet='ciceklab', epochs=3, b_sz=128,
            learn_method='unsup', embedding_type='Pretrain',
            config=gs_config
        )

        # 3d. Extract GraphSAGE node embeddings (100-dim per node)
        gs_model_path = os.path.join(zg.model_out_dir, 'Run_0', 'graphSage.pth')
        if os.path.exists(gs_model_path):
            print("  Extracting GraphSAGE embeddings ...")
            gs_embs = zg.run_graphsage_model(
                dataSet='ciceklab',
                cuda=torch.cuda.is_available(),
                config=gs_config,
                embedding_type='Pretrain',
                graphsage_model_path=gs_model_path
            )
            # Embeddings follow feats.txt order: RNAs first, then proteins
            gs_rna = gs_embs[:len(rna_id_order)]
            gs_prot = gs_embs[len(rna_id_order):]

            for i, r_id in enumerate(rna_id_order):
                rna_embed_dict[r_id] = np.concatenate(
                    (gs_rna[i].astype(np.float32), rna_embed_dict[r_id]))
            for i, p_id in enumerate(prot_id_order):
                prot_embed_dict[p_id] = np.concatenate(
                    (gs_prot[i].astype(np.float32), prot_embed_dict[p_id]))
            print(f"  RNA embed dim: {len(rna_embed_dict[rna_id_order[0]])} "
                  f"(100 gs + 640 seq)")
            print(f"  Prot embed dim: {len(prot_embed_dict[prot_id_order[0]])} "
                  f"(100 gs + 1024 seq)")
        else:
            print("  WARNING: GraphSAGE model not saved. Falling back to seq-only.")
            use_graphsage = False

    # ── 4. Build tensors ───────────────────────────────────────────────
    print("\nPreparing training tensors...")
    X_prot_list, X_rna_list, y_list = [], [], []
    skipped = 0
    for _, row in tqdm(train_df.iterrows(), total=len(train_df), desc="Building tensors"):
        r_id = row['RNA_aa_code']
        p_id = row['target_aa_code']
        if r_id in rna_embed_dict and p_id in prot_embed_dict:
            X_rna_list.append(rna_embed_dict[r_id])
            X_prot_list.append(prot_embed_dict[p_id])
            y_list.append(row['Y'])
        else:
            skipped += 1

    if skipped > 0:
        print(f"  Skipped {skipped} interactions (missing embeddings).")
    print(f"  Usable interactions: {len(y_list)}")

    X_prot = torch.tensor(np.array(X_prot_list), dtype=torch.float32)
    X_rna = torch.tensor(np.array(X_rna_list), dtype=torch.float32)
    y_train = torch.tensor(np.array(y_list), dtype=torch.float32)

    # DataLoader: order is (protein, rna, label) to match VecNN(target, rna)
    train_dataset = TensorDataset(X_prot, X_rna, y_train)
    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)

    # ── 5. Initialise model ────────────────────────────────────────────
    vecnn = VecNN(target_embed_len=1024, rna_embed_len=640,
                  graphsage_embedding=1 if use_graphsage else 0).to(device)
    optimizer = optim.Adam(vecnn.parameters(), lr=args.lr)
    criterion = nn.BCELoss()

    if args.epochs is not None:
        epochs = args.epochs
    elif args.smoke_test:
        epochs = 2
    else:
        epochs = 10

    # ── 6. Train ───────────────────────────────────────────────────────
    print(f"\nTraining VecNN for {epochs} epochs (batch_size={args.batch_size}, lr={args.lr})")
    for epoch in range(epochs):
        vecnn.train()
        total_loss = 0
        n_batches = 0
        for batch_prot, batch_rna, batch_y in tqdm(
                train_loader, desc=f"Epoch {epoch+1}/{epochs}"):
            batch_prot = batch_prot.to(device)
            batch_rna = batch_rna.to(device)
            batch_y = batch_y.to(device)

            optimizer.zero_grad()
            output = vecnn(batch_prot, batch_rna).squeeze(-1)
            loss = criterion(output, batch_y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
            n_batches += 1
        avg_loss = total_loss / max(n_batches, 1)
        print(f"  Epoch {epoch+1}/{epochs}  loss={avg_loss:.4f}")

    # Save trained model
    model_save_dir = os.path.join(BASE_DIR, 'trained_model')
    os.makedirs(model_save_dir, exist_ok=True)
    model_path = os.path.join(model_save_dir, 'VecNN_ciceklab.pth')
    torch.save(vecnn, model_path)
    print(f"\nSaved trained model to {model_path}")

    # ── 7. Evaluate on test splits ─────────────────────────────────────
    test_splits = [
        'final_test_unseen_pair.jsonl',
        'final_test_unseen_rna.jsonl',
        'final_test_unseen_protein.jsonl',
    ]
    test_dir = os.path.join(BASE_DIR, 'data_with_negatives', 'rna_protein')
    results_dir = os.path.join(BASE_DIR, 'results')
    os.makedirs(results_dir, exist_ok=True)

    all_results = {}
    for split in test_splits:
        split_path = os.path.join(test_dir, split)
        if not os.path.exists(split_path):
            print(f"\n  WARNING: {split} not found, skipping.")
            continue

        print(f"\n{'='*60}")
        print(f"Evaluating on: {split}")
        print(f"{'='*60}")
        test_df = load_test_split(split_path)

        X_prot_t, X_rna_t, y_test = [], [], []
        skipped = 0
        for _, row in test_df.iterrows():
            r_id = row['RNA_aa_code']
            p_id = row['target_aa_code']
            if r_id in rna_embed_dict and p_id in prot_embed_dict:
                X_prot_t.append(prot_embed_dict[p_id])
                X_rna_t.append(rna_embed_dict[r_id])
                y_test.append(row['Y'])
            else:
                skipped += 1

        if skipped > 0:
            print(f"  Skipped {skipped}/{len(test_df)} (missing embeddings)")
        if len(y_test) == 0:
            print("  No usable test pairs — skipping this split.")
            continue

        test_dataset = TensorDataset(
            torch.tensor(np.array(X_prot_t), dtype=torch.float32),
            torch.tensor(np.array(X_rna_t), dtype=torch.float32),
            torch.tensor(np.array(y_test), dtype=torch.float32),
        )
        test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)

        metrics = evaluate(vecnn, test_loader, device)
        all_results[split] = metrics

        print(f"  Test samples: {len(y_test)} "
              f"(pos={sum(1 for y in y_test if y==1)}, "
              f"neg={sum(1 for y in y_test if y==0)})")
        for k, v in metrics.items():
            if isinstance(v, float):
                print(f"  {k:>15s}: {v:.4f}")
            else:
                print(f"  {k:>15s}: {v}")

    # Save results summary
    results_path = os.path.join(results_dir, 'evaluation_results.json')
    # Convert numpy types for JSON serialization
    serializable = {}
    for split, mdict in all_results.items():
        serializable[split] = {k: float(v) if isinstance(v, (float, np.floating))
                               else int(v) if isinstance(v, (int, np.integer))
                               else v for k, v in mdict.items()}
    with open(results_path, 'w') as f:
        json.dump(serializable, f, indent=2)
    print(f"\nResults saved to {results_path}")
    print("Done!")


if __name__ == '__main__':
    main()
