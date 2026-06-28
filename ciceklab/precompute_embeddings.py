"""
Precompute RNA-FM and ProtTrans embeddings for the ciceklab dataset.

This script:
1. Scans training_chunks and data_with_negatives to find all unique RNA/protein IDs
2. Extracts their sequences from rna.fa / protein.fa
3. Uses ZHMolGraph's get_rnafm_embeddings and get_ProtTrans_embeddings to compute
   normalized embeddings, and saves them as pickle files.

The resulting pkl files contain DataFrames with columns:
  - RNA_aa_code / target_aa_code (the sequence)
  - normalized_embeddings (the 640/1024-dim vector)

We also save an ID-to-sequence mapping so train_ciceklab.py can look up embeddings by ID.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import pickle as pkl
from tqdm import tqdm

# Parse FASTA without BioPython dependency (avoids install issues)
def parse_fasta(fasta_file):
    """Simple FASTA parser returning dict of {id: sequence}."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:]  # full header after '>'
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            seqs[current_id] = ''.join(current_seq)
    return seqs


def get_unique_ids(chunk_dir, test_dir):
    """Scan all chunks and test files to find unique RNA and target IDs."""
    unique_rnas = set()
    unique_proteins = set()

    print("Scanning training chunks for unique IDs...")
    for file in tqdm(sorted(os.listdir(chunk_dir))):
        if file.endswith('.jsonl'):
            with open(os.path.join(chunk_dir, file), 'r') as f:
                for line in f:
                    obj = json.loads(line)
                    if obj.get('interaction_type') == 'rna-protein':
                        unique_rnas.add(obj['RNA_id'])
                        unique_proteins.add(obj['target_id'])

    print("Scanning test/validation files for unique IDs...")
    for file in tqdm(sorted(os.listdir(test_dir))):
        if file.endswith('.jsonl'):
            with open(os.path.join(test_dir, file), 'r') as f:
                for line in f:
                    obj = json.loads(line)
                    if obj.get('interaction_type') == 'rna-protein':
                        unique_rnas.add(obj['RNA_id'])
                        unique_proteins.add(obj['target_id'])

    return unique_rnas, unique_proteins


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    chunk_dir = os.path.join(base_dir, 'training_chunks')
    test_dir = os.path.join(base_dir, 'data_with_negatives', 'rna_protein')

    # ── 1. Collect unique IDs ──────────────────────────────────────────
    unique_rna_ids, unique_protein_ids = get_unique_ids(chunk_dir, test_dir)
    print(f"Found {len(unique_rna_ids)} unique RNA IDs and "
          f"{len(unique_protein_ids)} unique protein IDs.")

    # ── 2. Extract sequences from FASTA ────────────────────────────────
    print("Parsing rna.fa ...")
    all_rna_seqs = parse_fasta(os.path.join(base_dir, 'rna.fa'))
    print(f"  Total sequences in rna.fa: {len(all_rna_seqs)}")

    print("Parsing protein.fa ...")
    all_protein_seqs = parse_fasta(os.path.join(base_dir, 'protein.fa'))
    print(f"  Total sequences in protein.fa: {len(all_protein_seqs)}")

    # Filter to only sequences we need
    rna_seqs = {k: v for k, v in all_rna_seqs.items() if k in unique_rna_ids}
    protein_seqs = {k: v for k, v in all_protein_seqs.items() if k in unique_protein_ids}

    missing_rna = unique_rna_ids - set(rna_seqs.keys())
    missing_prot = unique_protein_ids - set(protein_seqs.keys())
    if missing_rna:
        print(f"  WARNING: {len(missing_rna)} RNA IDs not found in rna.fa")
    if missing_prot:
        print(f"  WARNING: {len(missing_prot)} protein IDs not found in protein.fa")

    print(f"Matched {len(rna_seqs)} RNA sequences and {len(protein_seqs)} protein sequences.")

    if len(rna_seqs) == 0 or len(protein_seqs) == 0:
        print("ERROR: No sequences found. Check that FASTA IDs match JSONL IDs.")
        return

    # ── 3. Save ID→sequence mapping ────────────────────────────────────
    id_map_path = os.path.join(base_dir, 'id_to_seq_map.pkl')
    with open(id_map_path, 'wb') as f:
        pkl.dump({'rna': rna_seqs, 'protein': protein_seqs}, f)
    print(f"Saved ID→sequence mapping to {id_map_path}")

    # ── 4. Build DataFrames for ZHMolGraph ─────────────────────────────
    # The ZHMolGraph class expects:
    #   rnas_dataframe with column 'RNA_aa_code' containing sequences
    #   proteins_dataframe with column 'target_aa_code' containing sequences
    # Keep insertion order so we can map embeddings back to IDs later.
    rna_id_list = list(rna_seqs.keys())
    rna_seq_list = [rna_seqs[rid] for rid in rna_id_list]
    protein_id_list = list(protein_seqs.keys())
    protein_seq_list = [protein_seqs[pid] for pid in protein_id_list]

    rna_df = pd.DataFrame({'RNA_aa_code': rna_seq_list})
    protein_df = pd.DataFrame({'target_aa_code': protein_seq_list})

    # The constructor also needs an interactions DF with the seq columns.
    # We pass an empty one just to satisfy the assertions.
    dummy_interactions = pd.DataFrame(columns=['RNA_aa_code', 'target_aa_code', 'Y'])

    sys.path.insert(0, os.path.join(base_dir, '..'))
    from ZHMolGraph.ZHMolGraph import ZHMolGraph

    vecnn_obj = ZHMolGraph(
        interactions=dummy_interactions,
        interaction_y_name='Y',
        rnas_dataframe=rna_df,
        rna_seq_name='RNA_aa_code',
        proteins_dataframe=protein_df,
        protein_seq_name='target_aa_code',
        model_out_dir=os.path.join(base_dir, 'trained_model'),
        debug=True
    )

    # ── 5. Compute RNA-FM embeddings ───────────────────────────────────
    # Call WITHOUT prediction_interactions so the method uses self.rna_list
    # and computes + stores mean/std normalization constants on self.
    rna_out_path = os.path.join(base_dir, 'rna_embeddings.pkl')
    if not os.path.exists(rna_out_path):
        print(f"\nComputing RNA-FM embeddings for {len(rna_df)} RNAs ...")
        print("  (This can take hours on CPU. Use a GPU if possible.)")
        vecnn_obj.get_rnafm_embeddings(
            prediction_interactions=None,
            embedding_dimension=640,
            replace_dataframe=True
        )
        # After this call:
        #   vecnn_obj.rnas_dataframe has columns [RNA_aa_code, normalized_embeddings]
        #   Each row's normalized_embeddings is a 640-dim numpy array.
        # We add the original ID back so train_ciceklab.py can look up by ID.
        result_df = vecnn_obj.rnas_dataframe.copy()
        result_df['RNA_id'] = rna_id_list
        with open(rna_out_path, 'wb') as f:
            pkl.dump(result_df, f)
        print(f"Saved RNA embeddings to {rna_out_path}")
    else:
        print(f"RNA embeddings already exist at {rna_out_path}, skipping.")

    # ── 6. Compute ProtTrans embeddings ────────────────────────────────
    prot_out_path = os.path.join(base_dir, 'protein_embeddings.pkl')
    if not os.path.exists(prot_out_path):
        print(f"\nComputing ProtTrans embeddings for {len(protein_df)} proteins ...")
        print("  (This can take hours on CPU. Use a GPU if possible.)")
        vecnn_obj.get_ProtTrans_embeddings(
            prediction_interactions=None,
            embedding_dimension=1024,
            replace_dataframe=True
        )
        result_df = vecnn_obj.proteins_dataframe.copy()
        result_df['target_id'] = protein_id_list
        with open(prot_out_path, 'wb') as f:
            pkl.dump(result_df, f)
        print(f"Saved protein embeddings to {prot_out_path}")
    else:
        print(f"Protein embeddings already exist at {prot_out_path}, skipping.")

    print("\nDone! You can now run train_ciceklab.py.")


if __name__ == '__main__':
    main()
