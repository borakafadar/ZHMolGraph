import os
import sys

def main():
    print("Running Smoke Test for ciceklab dataset...")
    # The smoke test could just run train_ciceklab on a very small subset of data,
    # or limit the number of epochs.
    # To implement this properly, we can call train_ciceklab.py with arguments,
    # but for simplicity, we will instruct the user to just run precompute_embeddings
    # and then run train_ciceklab with a smoke_test flag.
    
    cmd = f"{sys.executable} ciceklab/train_ciceklab.py --smoke_test"
    print(f"Executing: {cmd}")
    os.system(cmd)

if __name__ == '__main__':
    main()
