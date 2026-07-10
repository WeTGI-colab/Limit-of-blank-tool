#!/bin/bash
# Run the whole analysis on the server, end to end, then bundle the deliverables into a single
# tarball you can scp back to your machine.
#
#   bash scripts/run_all.sh
#
# It assumes the per-sample pile-ups already exist (results/pileup/*.tsv from the SLURM array).
# If they do not, either submit the array first (bash slurm/submit.sh) or set FULL=1 to compute
# them serially here (slower):  FULL=1 bash scripts/run_all.sh
#
# A python with pysam is picked automatically (falls back to /usr/bin/python3.11); override with
# PY=/path/to/python.  Plot/report steps are best-effort: if matplotlib/reportlab are missing on
# the server they are skipped, and you regenerate them locally from the bundled TSVs.
set -uo pipefail
cd "$(dirname "$0")/.."                      # repo root

die() { echo "ERROR: $*" >&2; exit 1; }
run() { echo "+ $*"; "$@" || die "step failed: $*"; }
try() { echo "+ (optional) $*"; "$@" || echo "  ...skipped (dependency missing?) -- run it locally from the TSVs"; }

# --- pick a python that can import pysam (needed only for the pile-up step) ---
PY="${PY:-python3}"
if ! "$PY" -c 'import pysam' 2>/dev/null; then
  for cand in /usr/bin/python3.11 python3.11 python3; do
    if command -v "$cand" >/dev/null 2>&1 && "$cand" -c 'import pysam' 2>/dev/null; then PY="$cand"; break; fi
  done
fi
echo "Using python: $PY  ($("$PY" --version 2>&1))"

# --- 1. per-sample pile-ups (skip if present; each task skips its own existing output) ---
if ! ls results/pileup/*.tsv >/dev/null 2>&1; then
  if [ "${FULL:-0}" = "1" ]; then
    n=$(($(wc -l < data/manifest.tsv) - 1))
    echo "Computing $n per-sample pile-ups serially (FULL=1)..."
    for i in $(seq 1 "$n"); do run "$PY" scripts/pileup_one_sample.py --index "$i"; done
  else
    die "no per-sample pile-ups in results/pileup/ -- run 'bash slurm/submit.sh' first, or set FULL=1"
  fi
fi

# --- 2. cohort tables (these NEED the server: BAMs and VCFs live here) ---
run "$PY" scripts/aggregate_cohort.py        # results/cohort_alt.tsv
run "$PY" scripts/export_vcf_variants.py     # results/vcf_variants{,_summary}.tsv
run "$PY" scripts/run_lob.py                 # results/lob_table.tsv (Gaussian + beta-binomial)
run "$PY" scripts/artefact_by_run.py --top 100   # results/artefact_by_run.tsv

# --- 3. figures and PDF (best-effort; regenerate locally if a dependency is missing) ---
try "$PY" scripts/plot_regions.py
try "$PY" scripts/plot_artefact_by_run.py
try "$PY" scripts/make_report.py             # results/artefact_report.pdf

# --- 4. bundle everything to bring back ---
bundle="results/deliverables.tar.gz"
paths=(results/cohort_alt.tsv results/lob_table.tsv results/vcf_variants.tsv
       results/vcf_variants_summary.tsv results/artefact_by_run.tsv)
[ -d results/plots ] && paths+=(results/plots)
[ -f results/artefact_report.pdf ] && paths+=(results/artefact_report.pdf)
run tar -czf "$bundle" "${paths[@]}"

echo
echo "Done. Bring this one file back to your machine:"
echo "  $(pwd)/$bundle"
echo "From your Mac:"
echo "  scp <user>@<server>:$(pwd)/$bundle .   &&   tar -xzf deliverables.tar.gz"
echo "If the PDF/plots were skipped on the server, regenerate them locally after extracting:"
echo "  python3 scripts/plot_regions.py && python3 scripts/plot_artefact_by_run.py && python3 scripts/make_report.py"
