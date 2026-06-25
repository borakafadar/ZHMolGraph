import os
import sys
import json
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
                             matthews_corrcoef, confusion_matrix)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from ZHMolGraph.ZHMolGraph_model import VecNN
from ZHMolGraph.ZHMolGraph import ZHMolGraph

def load_interactions(chunk_dir):
    """Load training chunks and return as DataFrame."""
    records = []
    print("Loading training chunks...")
    for file in tqdm(os.listdir(chunk_dir)):
        if file.endswith('.jsonl'):
            with open(os.path.join(chunk_dir, file), 'r') as f:
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

def evaluate(model, dataloader, device):
    """Evaluate VecNN model."""
    model.eval()
    all_preds = []
    all_labels = []
    with torch.no_grad():
        for batch_rna, batch_prot, batch_y in dataloader:
            batch_rna = batch_rna.to(device)
            batch_prot = batch_prot.to(device)
            
            output = model(batch_prot, batch_rna)
            all_preds.extend(output.squeeze().cpu().numpy().tolist())
            all_labels.extend(batch_y.numpy().tolist())
            
    all_preds = np.array(all_preds)
    all_labels = np.array(all_labels)
    
    auc = roc_auc_score(all_labels, all_preds)
    auprc = average_precision_score(all_labels, all_preds)
    
    # Calculate best MCC threshold
    thresholds = np.linspace(0, 1, 100)
    best_mcc = -1
    best_thresh = 0.5
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
    
    tn, fp, fn, tp = confusion_matrix(all_labels, preds_bin).ravel()
    spe = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    return {
        'AUC': auc, 'AUPRC': auprc, 'Accuracy': acc, 'Precision': prec,
        'Sensitivity': rec, 'Specificity': spe, 'MCC': best_mcc, 'Threshold': best_thresh
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smoke_test', action='store_true', help='Run a quick smoke test on a subset of data')
    parser.add_argument('--use_graphsage', action='store_true', help='Use GraphSAGE embeddings as stated in paper')
    args = parser.parse_args()

    base_dir = os.path.dirname(__file__)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # 1. Load precomputed embeddings
    rna_embed_path = os.path.join(base_dir, 'rna_embeddings.pkl')
    prot_embed_path = os.path.join(base_dir, 'protein_embeddings.pkl')
    
    if not os.path.exists(rna_embed_path) or not os.path.exists(prot_embed_path):
        print("Embeddings not found. Please run precompute_embeddings.py first.")
        return
        
    with open(rna_embed_path, 'rb') as f:
        rna_embeds_df = pkl.load(f)
    with open(prot_embed_path, 'rb') as f:
        prot_embeds_df = pkl.load(f)
        
    print(f"Loaded {len(rna_embeds_df)} RNA embeddings and {len(prot_embeds_df)} Protein embeddings.")

    # Convert to dict for fast lookup
    # Note: RNA_id column added in precompute script stores the original ID
    rna_embed_dict = {row['RNA_id']: row['normalized_embeddings'] for _, row in rna_embeds_df.iterrows()}
    prot_embed_dict = {row['target_id']: row['normalized_embeddings'] for _, row in prot_embeds_df.iterrows()}

    # 2. Load dataset
    chunk_dir = os.path.join(base_dir, 'training_chunks')
    train_df = load_interactions(chunk_dir)
    
    if args.smoke_test:
        print("Smoke test enabled: Subsetting dataset...")
        train_df = train_df.head(1000)
        
    print(f"Loaded {len(train_df)} training interactions.")
    
    # Decide if we use GraphSAGE based on arguments
    use_graphsage = args.use_graphsage
    
    # Check if we need to generate GraphSAGE interactions and features
    if use_graphsage:
        print("Generating GraphSAGE datasets (feats and interactions)...")
        gs_dir = os.path.join(base_dir, '..', 'ciceklab_Pretrain_graphsage_dataset')
        os.makedirs(gs_dir, exist_ok=True)
        
        # 1. Generate feats.txt
        feats_path = os.path.join(gs_dir, 'ciceklab_feats.txt')
        if not os.path.exists(feats_path):
            print("Creating ciceklab_feats.txt...")
            with open(feats_path, 'w') as f:
                # RNAs (padded to 1024)
                for r_id, emb in rna_embed_dict.items():
                    padded = np.pad(emb, (0, 1024 - len(emb)))
                    row = [r_id] + padded.tolist()
                    f.write('\t'.join(map(str, row)) + '\n')
                # Proteins (already 1024)
                for p_id, emb in prot_embed_dict.items():
                    row = [p_id] + emb.tolist()
                    f.write('\t'.join(map(str, row)) + '\n')
                    
        # 2. Generate interaction txts
        # Total interactions
        train_df[['RNA_aa_code', 'target_aa_code']].to_csv(
            os.path.join(gs_dir, 'ciceklab_total_interactions_seq_list.txt'),
            sep='\t', index=False, header=False
        )
        
        # Positive interactions for train
        pos_train = train_df[train_df['Y'] == 1]
        pos_train[['RNA_aa_code', 'target_aa_code']].to_csv(
            os.path.join(gs_dir, 'ciceklab_graphsage_train_interactions.txt'),
            sep='\t', index=False, header=False
        )
        
        # Negative interactions plus all test (dummy test for GraphSAGE training)
        # We will just put a dummy test list here, since GraphSAGE is unsupervised on structure
        neg_train = train_df[train_df['Y'] == 0]
        neg_train[['RNA_aa_code', 'target_aa_code']].to_csv(
            os.path.join(gs_dir, 'ciceklab_graphsage_test_interactions.txt'),
            sep='\t', index=False, header=False
        )
        
        print("Training GraphSAGE model...")
        sys.path.append(os.path.join(base_dir, '..'))
        # Create ZHMolGraph object to use its run_graphsage_experiment method
        # It requires setting up properties first
        zg = ZHMolGraph()
        zg.model_out_dir = os.path.join(base_dir, 'trained_model', 'ciceklab')
        os.makedirs(os.path.join(zg.model_out_dir, 'Run_0'), exist_ok=True)
        
        zg.run_graphsage_experiment(
            dataSet='ciceklab_Pretrain', epochs=3, b_sz=128, 
            learn_method='unsup', embedding_type='Pretrain',
            config=os.path.join(base_dir, '..', 'graphsage_src', 'experiments.conf')
        )
        print("Finished GraphSAGE training. Extracting embeddings...")
        
        # Extract embeddings for the training set nodes
        gs_model_path = os.path.join(zg.model_out_dir, 'Run_0', 'graphSage.pth')
        gs_embs = zg.run_graphsage_model(
            dataSet='ciceklab_Pretrain', cuda=torch.cuda.is_available(), 
            config=os.path.join(base_dir, '..', 'graphsage_src', 'experiments.conf'),
            embedding_type='Pretrain', graphsage_model_path=gs_model_path
        )
        
        # GraphSAGE returns embeddings in the order nodes appear in feats.txt
        # Our feats.txt wrote RNAs first, then Proteins
        num_rnas = len(rna_embed_dict)
        num_prots = len(prot_embed_dict)
        
        gs_rna_embs = gs_embs[:num_rnas]
        gs_prot_embs = gs_embs[num_rnas:]
        
        # Update rna_embed_dict and prot_embed_dict with concatenated embeddings
        # We need to map them back correctly
        for i, (r_id, emb) in enumerate(rna_embed_dict.items()):
            rna_embed_dict[r_id] = np.concatenate((gs_rna_embs[i], emb))
            
        for i, (p_id, emb) in enumerate(prot_embed_dict.items()):
            prot_embed_dict[p_id] = np.concatenate((gs_prot_embs[i], emb))

    vecnn = VecNN(target_embed_len=1024, rna_embed_len=640, graphsage_embedding=1 if use_graphsage else 0).to(device)
    optimizer = optim.Adam(vecnn.parameters(), lr=0.001)
    criterion = nn.BCELoss()
    
    # 3. Prepare Training DataLoader
    print("Preparing training tensors...")
    X_rna = []
    X_prot = []
    y_train = []
    for _, row in tqdm(train_df.iterrows(), total=len(train_df)):
        r_id = row['RNA_aa_code']
        p_id = row['target_aa_code']
        if r_id in rna_embed_dict and p_id in prot_embed_dict:
            X_rna.append(rna_embed_dict[r_id])
            X_prot.append(prot_embed_dict[p_id])
            y_train.append(row['Y'])
            
    X_rna = torch.tensor(np.array(X_rna), dtype=torch.float32)
    X_prot = torch.tensor(np.array(X_prot), dtype=torch.float32)
    y_train = torch.tensor(np.array(y_train), dtype=torch.float32)
    
    train_dataset = TensorDataset(X_rna, X_prot, y_train)
    train_loader = DataLoader(train_dataset, batch_size=256, shuffle=True)
    
    # 4. Train Model
    epochs = 2 if args.smoke_test else 10
    print("Starting training...")
    for epoch in range(epochs):
        vecnn.train()
        total_loss = 0
        for batch_rna, batch_prot, batch_y in tqdm(train_loader, desc=f"Epoch {epoch+1}/{epochs}"):
            batch_rna = batch_rna.to(device)
            batch_prot = batch_prot.to(device)
            batch_y = batch_y.to(device)
            
            optimizer.zero_grad()
            output = vecnn(batch_prot, batch_rna)
            loss = criterion(output.squeeze(), batch_y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        print(f"Epoch {epoch+1} Loss: {total_loss/len(train_loader):.4f}")
        
    # 5. Evaluate on Test Splits
    test_splits = [
        'final_test_unseen_pair.jsonl',
        'final_test_unseen_rna.jsonl',
        'final_test_unseen_protein.jsonl'
    ]
    test_dir = os.path.join(base_dir, 'data_with_negatives', 'rna_protein')
    
    for split in test_splits:
        print(f"\nEvaluating on {split}...")
        test_df = load_test_split(os.path.join(test_dir, split))
        X_rna_t, X_prot_t, y_test = [], [], []
        for _, row in test_df.iterrows():
            r_id = row['RNA_aa_code']
            p_id = row['target_aa_code']
            if r_id in rna_embed_dict and p_id in prot_embed_dict:
                X_rna_t.append(rna_embed_dict[r_id])
                X_prot_t.append(prot_embed_dict[p_id])
                y_test.append(row['Y'])
                
        if len(y_test) == 0:
            print("No matching IDs found for this split.")
            continue
            
        test_dataset = TensorDataset(
            torch.tensor(np.array(X_rna_t), dtype=torch.float32),
            torch.tensor(np.array(X_prot_t), dtype=torch.float32),
            torch.tensor(np.array(y_test), dtype=torch.float32)
        )
        test_loader = DataLoader(test_dataset, batch_size=256, shuffle=False)
        
        metrics = evaluate(vecnn, test_loader, device)
        for k, v in metrics.items():
            print(f"{k}: {v:.4f}")

if __name__ == '__main__':
    main()
