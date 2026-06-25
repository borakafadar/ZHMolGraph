import os
import sys

def main():
    print("Running Full Training for ciceklab dataset...")
    # Executes the full training loop with GraphSAGE embeddings
    cmd = f"{sys.executable} ciceklab/train_ciceklab.py --use_graphsage"
    print(f"Executing: {cmd}")
    os.system(cmd)

if __name__ == '__main__':
    main()
