#!/usr/bin/env python3
"""Retrieve the GRCh38 reference sequence for each amplicon interval in the panel BED.

The full human reference is not required for local development. Instead, the exact
sequence spanning each amplicon (chromStart..chromEnd, i.e. including primer flanks) is
fetched from the UCSC REST API and stored as a compact multi-FASTA. This provides the
reference bases needed to (a) synthesise reads and (b) determine the non-reference alleles
during pile-up, while keeping the repository free of large genomic assets.

Output: reference/regions.fa  (one record per unique amplicon interval,
         header = "chrom:chromStart-chromEnd", 0-based half-open, matching the BED).
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BED = REPO / "bed" / "panel.bed"
OUT = REPO / "reference" / "regions.fa"
UCSC = "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"


def unique_intervals(bed_path):
    """Collapse the forward/reverse primer rows to unique (chrom, start, end) intervals."""
    seen = {}
    for line in bed_path.read_text().splitlines():
        if not line.strip():
            continue
        chrom, start, end = line.split("\t")[:3]
        key = (chrom, int(start), int(end))
        seen[key] = None
    return sorted(seen, key=lambda k: (k[0], k[1]))


def fetch(chrom, start, end):
    url = UCSC.format(chrom=chrom, start=start, end=end)
    with urllib.request.urlopen(url, timeout=30) as resp:
        payload = json.load(resp)
    dna = payload["dna"].upper()
    if len(dna) != end - start:
        raise ValueError(f"{chrom}:{start}-{end} expected {end - start} bp, got {len(dna)}")
    return dna


def main():
    intervals = unique_intervals(BED)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for chrom, start, end in intervals:
        dna = fetch(chrom, start, end)
        records.append(f">{chrom}:{start}-{end}\n{dna}")
        print(f"  {chrom}:{start}-{end}  {len(dna)} bp", file=sys.stderr)
        time.sleep(0.2)  # be polite to the public API
    OUT.write_text("\n".join(records) + "\n")
    print(f"Wrote {len(records)} intervals to {OUT.relative_to(REPO)}", file=sys.stderr)


if __name__ == "__main__":
    main()
