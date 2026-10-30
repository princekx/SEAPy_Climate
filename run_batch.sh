#!/bin/bash
#SBATCH --job-name=seapy_batch
#SBATCH --output=seapy_batch.out
#SBATCH --error=seapy_batch.err
#SBATCH --time=06:00:00
#SBATCH --mem=80gb
#SBATCH --cpus-per-task=8
#/data/apps/sss/environments/default-2026_06_10/bin/python main.py --runs obs control experiment --mjo-levels 1 2 3 4 --force
/data/apps/sss/environments/default-2026_06_10/bin/python main.py --runs obs control experiment --html-only