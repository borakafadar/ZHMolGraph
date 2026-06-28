"""Full training: train VecNN with GraphSAGE on the complete ciceklab dataset."""
import os
import sys

def main():
    print("=" * 60)
    print("ZHMolGraph Full Training  (with GraphSAGE)")
    print("=" * 60)
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    cmd = f'"{sys.executable}" ciceklab/train_ciceklab.py --use_graphsage'
    print(f"Running: {cmd}\n")
    os.chdir(repo_dir)
    exit_code = os.system(cmd)
    sys.exit(exit_code >> 8)

if __name__ == '__main__':
    main()
