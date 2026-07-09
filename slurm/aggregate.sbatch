#!/bin/bash
# Combine all per-sample tables, build the LoB model and the plots. Submitted by slurm/submit.sh
# with a dependency so it only runs after every per-sample task has succeeded.
#SBATCH --job-name=lob_aggregate
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/aggregate_%j.out
#SBATCH --error=logs/aggregate_%j.err
# --- adjust for your cluster ---
# #SBATCH --partition=YOUR_PARTITION
# #SBATCH --account=YOUR_ACCOUNT

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

# module load python/3.9
# source .venv/bin/activate

python3 scripts/aggregate_cohort.py
python3 scripts/run_lob.py
python3 scripts/plot_regions.py
python3 scripts/artefact_by_run.py --top 100
python3 scripts/plot_artefact_by_run.py
