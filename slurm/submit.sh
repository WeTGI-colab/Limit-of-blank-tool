#!/bin/bash
# Submit the whole pile-up -> aggregate -> LoB pipeline to SLURM:
#   1. an array job with one task per sample in data/manifest.tsv
#   2. an aggregation job that runs only after the array succeeds (afterok dependency)
#
# Prerequisites (run once, on a login node or interactively):
#   python3 scripts/discover_samples.py --data-path /path/to/data --ids-file config/run_ids.txt
#   python3 scripts/fetch_reference.py --ref-fasta /path/to/GRCh38.fa
#
# Usage:  bash slurm/submit.sh [max_concurrent]     (default 50 concurrent tasks)
set -euo pipefail
cd "$(dirname "$0")/.."

MAXC="${1:-50}"
N=$(( $(wc -l < data/manifest.tsv) - 1 ))
if [ "$N" -lt 1 ]; then echo "Empty manifest -- run discover_samples.py first."; exit 1; fi
mkdir -p logs

echo "Submitting ${N} per-sample pile-up tasks (max ${MAXC} concurrent)..."
JID=$(sbatch --parsable --array="1-${N}%${MAXC}" slurm/pileup_array.sbatch)
echo "  array job: ${JID}"

AGG=$(sbatch --parsable --dependency="afterok:${JID}" slurm/aggregate.sbatch)
echo "  aggregation job (after array): ${AGG}"
echo "Watch with: squeue -u \$USER"
