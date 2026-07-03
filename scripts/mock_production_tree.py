#!/usr/bin/env python3
"""Lay the synthetic cohort out as a production-style directory tree, to validate the full
real-data path (discovery -> manifest -> pile-up -> LoB) end to end without production access.

It mirrors the production layout described by the lab:
  <data>/<something ending with a run ID>/<sample subfolder>/
        s<name>_FINAL.bam   (+ .bai)
        <name>_FINAL.vcf         (paired VCF, without the leading 's')

Symlinks point back to data/synthetic so no data is duplicated. Run discovery afterwards:
  python3 scripts/discover_samples.py --data-path data/prod_like --ids-file config/run_ids.txt
"""
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SYN = REPO / "data" / "synthetic"
OUT = REPO / "data" / "prod_like"
IDS_FILE = REPO / "config" / "run_ids.txt"
N_RUNS = 4                                  # spread samples across this many run folders


def main():
    samples = sorted(p.name for p in SYN.glob("sample_*") if (p / f"{p.name}.bam").exists())
    if not samples:
        sys.exit("No synthetic samples found; run scripts/make_synthetic_data.py first.")
    ids = [ln.strip() for ln in IDS_FILE.read_text().splitlines() if ln.strip()][:N_RUNS]

    if OUT.exists():
        shutil.rmtree(OUT)
    for i, sname in enumerate(samples):
        run_id = ids[i % len(ids)]
        core = f"260626_M01875_1267_000000000-{run_id}_W{2614000 + i}-002"
        sub = OUT / f"project_{run_id}" / f"{core}"       # run folder ends with the ID
        sub.mkdir(parents=True, exist_ok=True)
        src_bam = (SYN / sname / f"{sname}.bam").resolve()
        src_bai = (SYN / sname / f"{sname}.bam.bai").resolve()
        src_vcf = (SYN / sname / f"{sname}.vcf").resolve()
        (sub / f"s{core}_FINAL.bam").symlink_to(src_bam)
        (sub / f"s{core}_FINAL.bam.bai").symlink_to(src_bai)
        (sub / f"{core}_FINAL.vcf").symlink_to(src_vcf)
        print(f"  {sname} -> project_{run_id}/{core}/  (s{core}_FINAL.bam)")

    print(f"\nBuilt {len(samples)} samples across {min(len(ids), len(samples))} run folders "
          f"in {OUT.relative_to(REPO)}")
    print("Next: python3 scripts/discover_samples.py --data-path data/prod_like "
          "--ids-file config/run_ids.txt")


if __name__ == "__main__":
    main()
