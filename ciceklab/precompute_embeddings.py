import os
import sys
import json
import numpy as np
import pandas as pd
import pickle as pkl
from tqdm import tqdm
from Bio import SeqIO

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from ZHMolGraph import ZHMolGraph

def get_unique_ids(chunk_dir, test_dir):
    """Scan all chunks and test files to find unique RNA and target IDs"""
    unique_rnas = set()
    unique_proteins = set()
    
    print("Scanning training chunks for unique IDs...")
    for file in tqdm(os.listdir(chunk_dir)):
        if file.endswith('.jsonl'):
            with open(os.path.join(chunk_dir, file), 'r') as f:
                for line in f:
                    obj = json.loads(line)
                    if obj.get('interaction_type') == 'rna-protein':
                        unique_rnas.add(obj['RNA_id'])
                        unique_proteins.add(obj['target_id'])
                        
    print("Scanning test files for unique IDs...")
    for file in tqdm(os.listdir(test_dir)):
        if file.endswith('.jsonl'):
            with open(os.path.join(test_dir, file), 'r') as f:
                for line in f:
                    obj = json.loads(line)
                    if obj.get('interaction_type') == 'rna-protein':
                        unique_rnas.add(obj['RNA_id'])
                        unique_proteins.add(obj['target_id'])
                        
    return unique_rnas, unique_proteins

def extract_sequences(fasta_file, target_ids):
    """Extract sequences from FASTA file for given IDs"""
    seqs = {}
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Some IDs have "::" or other variations but the jsonl ID should match the FASTA ID exactly
            if record.id in target_ids:
                seqs[record.id] = str(record.seq)
    return seqs

def main():
    base_dir = os.path.dirname(__file__)
    chunk_dir = os.path.join(base_dir, 'training_chunks')
    test_dir = os.path.join(base_dir, 'data_with_negatives', 'rna_protein')
    
    # 1. Get unique IDs
    unique_rnas, unique_proteins = get_unique_ids(chunk_dir, test_dir)
    print(f"Found {len(unique_rnas)} unique RNAs and {len(unique_proteins)} unique Proteins in dataset.")
    
    # 2. Extract sequences
    print("Extracting RNA sequences...")
    rna_seqs = extract_sequences(os.path.join(base_dir, 'rna.fa'), unique_rnas)
    print(f"Extracted {len(rna_seqs)} RNA sequences.")
    
    print("Extracting Protein sequences...")
    protein_seqs = extract_sequences(os.path.join(base_dir, 'protein.fa'), unique_proteins)
    print(f"Extracted {len(protein_seqs)} Protein sequences.")
    
    # 3. Create DataFrames for ZHMolGraph
    rna_df = pd.DataFrame([{'RNA_aa_code': seq} for seq in rna_seqs.values()])
    protein_df = pd.DataFrame([{'target_aa_code': seq} for seq in protein_seqs.values()])
    
    # Check if we have sequences to compute
    if rna_df.empty or protein_df.empty:
        print("Error: No sequences found. Please check FASTA ID matching.")
        return
        
    print(f"Preparing to compute embeddings for {len(rna_df)} RNAs and {len(protein_df)} Proteins.")
    
    # 4. Initialize ZHMolGraph
    # Dummy paths as we don't have interactions yet
    vecnn_object = ZHMolGraph.ZHMolGraph(
        interactions_location=None,
        interactions=pd.DataFrame(columns=['RNA_aa_code', 'target_aa_code', 'Y']),
        interaction_y_name='Y',
        rnas_dataframe=rna_df,
        rna_seq_name='RNA_aa_code',
        proteins_dataframe=protein_df,
        protein_seq_name='target_aa_code',
        model_out_dir=os.path.join(base_dir, 'trained_model', 'ZHMolGraph_VecNN_model_RPI_ciceklab'),
        debug=False
    )
    
    # 5. Compute RNA embeddings
    rna_out_path = os.path.join(base_dir, 'rna_embeddings.pkl')
    if not os.path.exists(rna_out_path):
        print("Computing RNA-FM embeddings (this may take a long time)...")
        rna_embeds_df = vecnn_object.get_rnafm_embeddings(prediction_interactions=rna_df, replace_dataframe=False)
        # Store back the mapping from ID to sequence so we can easily join during training
        rna_embeds_df['RNA_id'] = list(rna_seqs.keys())
        with open(rna_out_path, 'wb') as f:
            pkl.dump(rna_embeds_df, f)
        print(f"Saved RNA embeddings to {rna_out_path}")
    else:
        print(f"RNA embeddings already exist at {rna_out_path}")
        
    # 6. Compute Protein embeddings
    prot_out_path = os.path.join(base_dir, 'protein_embeddings.pkl')
    if not os.path.exists(prot_out_path):
        print("Computing ProtTrans embeddings (this may take a long time)...")
        prot_embeds_df = vecnn_object.get_ProtTrans_embeddings(prediction_interactions=protein_df, replace_dataframe=False)
        prot_embeds_df['target_id'] = list(protein_seqs.keys())
        with open(prot_out_path, 'wb') as f:
            pkl.dump(prot_embeds_df, f)
        print(f"Saved Protein embeddings to {prot_out_path}")
    else:
        print(f"Protein embeddings already exist at {prot_out_path}")

if __name__ == '__main__':
    main()
