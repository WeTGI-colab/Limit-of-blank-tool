#!/usr/bin/env python3
"""Discover production BAM/VCF pairs and write data/manifest.tsv.

Given a base data directory and a file of run IDs (one per line, e.g. config/run_ids.txt),
locate every folder directly under the base whose name ends with an ID, recurse into it, and
pair each ``*_FINAL.bam`` with its ``*_FINAL.vcf``. The BAM carries a leading 's' that the VCF
does not (``s<name>_FINAL.bam`` <-> ``<name>_FINAL.vcf``); a single ``*_FINAL.vcf`` in the same
directory is used as a fallback. Each pair is one sample; every sample found is written to the
manifest and aggregated together downstream.

Usage:
    python3 scripts/discover_samples.py --data-path /path/to/data --ids-file config/run_ids.txt
"""
import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "manifest.tsv"


def find_vcf(bam):
    """Locate the VCF paired with a *_FINAL.bam in the same directory."""
    d, name = bam.parent, bam.name
    candidates = []
    if name.startswith("s"):
        candidates.append(d / (name[1:-4] + ".vcf"))     # drop leading 's', .bam -> .vcf
    candidates.append(d / (name[:-4] + ".vcf"))
    for c in candidates:
        if c.exists():
            return c
    loose = sorted(d.glob("*_FINAL.vcf"))
    return loose[0] if len(loose) == 1 else None


def sample_id(bam):
    core = bam.name[:-4]                                  # strip '.bam'
    if core.startswith("s"):
        core = core[1:]
    return core[:-6] if core.endswith("_FINAL") else core


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-path", required=True, help="base directory holding the run folders")
    ap.add_argument("--ids-file", default=str(REPO / "config" / "run_ids.txt"),
                    help="file of run IDs, one per line")
    args = ap.parse_args()

    base = Path(args.data_path)
    if not base.is_dir():
        sys.exit(f"data-path not found: {base}")
    ids = [ln.strip() for ln in Path(args.ids_file).read_text().splitlines() if ln.strip()]

    rows, missing_id, missing_vcf = [], [], []
    for run_id in ids:
        run_dirs = [d for d in base.iterdir() if d.is_dir() and d.name.endswith(run_id)]
        if not run_dirs:
            missing_id.append(run_id)
            continue
        for rd in run_dirs:
            for bam in sorted(rd.rglob("*_FINAL.bam")):
                vcf = find_vcf(bam)
                if vcf is None:
                    missing_vcf.append(bam.name)
                rows.append((sample_id(bam), run_id, str(bam.resolve()),
                             str(vcf.resolve()) if vcf else ""))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("sample\trun\tbam\tvcf\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print(f"{len(rows)} sample(s) across {len(ids) - len(missing_id)} run ID(s) "
          f"-> {OUT.relative_to(REPO)}")
    if missing_id:
        print(f"WARNING: no folder ending with these IDs: {', '.join(missing_id)}")
    if missing_vcf:
        print(f"WARNING: {len(missing_vcf)} BAM(s) without a paired VCF")


if __name__ == "__main__":
    main()
