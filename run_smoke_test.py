"""Smoke test: train VecNN on a tiny subset of the ciceklab data."""
import os
import sys

def main():
    print("=" * 60)
    print("ZHMolGraph Smoke Test  (1000 samples, 2 epochs)")
    print("=" * 60)
    # Run from repo root so graphsage_src imports work
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    cmd = f'"{sys.executable}" ciceklab/train_ciceklab.py --smoke_test'
    print(f"Running: {cmd}\n")
    os.chdir(repo_dir)
    exit_code = os.system(cmd)
    sys.exit(exit_code >> 8)  # propagate exit code

if __name__ == '__main__':
    main()
