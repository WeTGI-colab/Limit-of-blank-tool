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

# pick a python with the full stack the aggregate needs (pandas, numpy AND scipy -- run_lob uses
# scipy); a conda 'base' with only pandas is not enough. Override with PY=/path/to/python
PY="${PY:-python3}"
if ! "$PY" -c 'import pandas, numpy, scipy' 2>/dev/null; then
  for cand in /usr/bin/python3.11 python3.11; do
    if command -v "$cand" >/dev/null 2>&1 && "$cand" -c 'import pandas, numpy, scipy' 2>/dev/null; then PY="$cand"; break; fi
  done
fi
"$PY" -c 'import pandas, numpy, scipy' 2>/dev/null || { echo "ERROR: no python with pandas+numpy+scipy found (set PY=)"; exit 1; }

# raw data tables only -- plots and the PDF are generated locally from these
"$PY" scripts/aggregate_cohort.py            # results/cohort_alt.tsv
"$PY" scripts/export_vcf_variants.py         # results/vcf_variants{,_summary}.tsv
"$PY" scripts/run_lob.py                     # results/lob_table.tsv (Gaussian + beta-binomial)
"$PY" scripts/artefact_by_run.py --top 100   # results/artefact_by_run.tsv

# bundle the tables to bring back
tar -czf results/deliverables.tar.gz \
  results/cohort_alt.tsv results/lob_table.tsv results/vcf_variants.tsv \
  results/vcf_variants_summary.tsv results/artefact_by_run.tsv
echo "Wrote results/deliverables.tar.gz -- scp this back."
