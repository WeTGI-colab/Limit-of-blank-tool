#!/bin/bash
# One SLURM array task per sample. Submitted by slurm/submit.sh (which sets --array).
#SBATCH --job-name=lob_pileup
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/pileup_%A_%a.out
#SBATCH --error=logs/pileup_%A_%a.err
# --- adjust for your cluster ---
# #SBATCH --partition=YOUR_PARTITION
# #SBATCH --account=YOUR_ACCOUNT

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

# make Python available (adjust to your cluster: module load / conda activate / venv)
# module load python/3.9
# source .venv/bin/activate

# pick a python that can import pysam (the pile-up needs it); override with PY=/path/to/python
PY="${PY:-python3}"
if ! "$PY" -c 'import pysam' 2>/dev/null; then
  for cand in /usr/bin/python3.11 python3.11; do
    if command -v "$cand" >/dev/null 2>&1 && "$cand" -c 'import pysam' 2>/dev/null; then PY="$cand"; break; fi
  done
fi
"$PY" -c 'import pysam' 2>/dev/null || { echo "ERROR: no python with pysam found (set PY=)"; exit 1; }

"$PY" scripts/pileup_one_sample.py --index "${SLURM_ARRAY_TASK_ID}"
