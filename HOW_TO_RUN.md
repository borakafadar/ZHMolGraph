# ZHMolGraph Evaluation on ciceklab Dataset

## Overview

This document explains the scripts added to train and evaluate the `ZHMolGraph` model on the custom `ciceklab` RNA-protein interaction dataset.

## Changes Made

### New files created

| File | Purpose |
|------|---------|
| `ciceklab/precompute_embeddings.py` | Compute RNA-FM (640-dim) and ProtTrans (1024-dim) embeddings for all unique sequences |
| `ciceklab/train_ciceklab.py` | Train VecNN and evaluate on 3 unseen test splits |
| `run_smoke_test.py` | Wrapper: quick local test (1000 samples, 2 epochs) |
| `run_full_training.py` | Wrapper: full training with GraphSAGE |

### Existing files modified

| File | Change |
|------|--------|
| `graphsage_src/experiments.conf` | Added `ciceklab_Pretrain_*` paths for GraphSAGE dataset |

## Pipeline

The ZHMolGraph pipeline has two stages:

1. **Embedding extraction** — RNA-FM for RNA sequences (640-dim) and ProtTrans for protein sequences (1024-dim). This is the compute-heavy step.
2. **Model training** — Train the VecNN (1D-CNN) classifier that takes RNA + protein embeddings as input and predicts interaction probability.

Optionally, GraphSAGE (unsupervised GNN) can be used to generate 100-dim structure-aware node embeddings from the interaction graph, which get concatenated with the sequence embeddings before VecNN training. This is the full method described in the paper.

## How to Run

### Step 1: Precompute Embeddings

This step extracts LLM embeddings for every unique RNA and protein in your dataset. It only needs to be run once. Subsequent runs are skipped if the output files already exist.

```bash
# From the repo root:
python ciceklab/precompute_embeddings.py
```

**Output files** (saved in `ciceklab/`):
- `rna_embeddings.pkl` — DataFrame with columns `[RNA_aa_code, normalized_embeddings, RNA_id]`
- `protein_embeddings.pkl` — DataFrame with columns `[target_aa_code, normalized_embeddings, target_id]`
- `id_to_seq_map.pkl` — Mapping from IDs to sequences

> **Note:** This can take many hours on CPU. Use a GPU server if possible.

### Step 2: Train and Evaluate

**Smoke Test (local, CPU-friendly):**
```bash
python run_smoke_test.py
```
This trains on 1,000 interactions for 2 epochs. Use this to verify the pipeline works before running the full training.

**Full Training (GPU server, with GraphSAGE):**
```bash
python run_full_training.py
```
This trains the full pipeline including GraphSAGE unsupervised graph embeddings, as described in the ZHMolGraph paper.

**Direct usage with custom options:**
```bash
python ciceklab/train_ciceklab.py --epochs 20 --batch_size 512 --lr 0.0005
python ciceklab/train_ciceklab.py --smoke_test
python ciceklab/train_ciceklab.py --use_graphsage --epochs 15
```

### Step 3: Results

After evaluation, results are saved to:
- **Console output** — Full metrics for each test split
- **`ciceklab/results/evaluation_results.json`** — Machine-readable results
- **`ciceklab/trained_model/VecNN_ciceklab.pth`** — Saved model checkpoint

The model is evaluated on all three test splits:
- `final_test_unseen_pair.jsonl` — Unseen RNA-protein pairs
- `final_test_unseen_rna.jsonl` — Unseen RNA sequences
- `final_test_unseen_protein.jsonl` — Unseen protein sequences

## Environment

No additional packages are required beyond what's in the existing `requirements.txt`. The scripts use: `torch`, `pandas`, `tqdm`, `transformers`, `scikit-learn`, `numpy`, `fm` (RNA-FM).
