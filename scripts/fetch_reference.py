#!/usr/bin/env python3
"""Extract the GRCh38 reference sequence for each amplicon interval in the panel BED.

The full reference is not stored in the repository. Region sequences are obtained either from
a local reference FASTA (recommended on the analysis server -- use the same GRCh38 build that
Pisces uses) or, as a fallback, from the UCSC REST API when the machine has internet access.

    # from a local FASTA (must have a .fai index alongside it)
    python3 scripts/fetch_reference.py --ref-fasta /path/to/GRCh38.fa

    # or from UCSC (needs internet)
    python3 scripts/fetch_reference.py

Output: reference/regions.fa -- one record per unique amplicon interval, header
        "chrom:chromStart-chromEnd" (panel convention, 0-based half-open, matching the BED).
"""
import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob.samples import matching_contig  # noqa: E402

BED = REPO / "bed" / "panel.bed"
OUT = REPO / "reference" / "regions.fa"
UCSC = "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"


def unique_intervals(bed_path):
    """Collapse forward/reverse primer rows to unique (chrom, start, end) intervals."""
    seen = {}
    for line in bed_path.read_text().splitlines():
        if not line.strip():
            continue
        chrom, start, end = line.split("\t")[:3]
        seen[(chrom, int(start), int(end))] = None
    return sorted(seen, key=lambda k: (k[0], k[1]))


def fetch_ucsc(chrom, start, end):
    url = UCSC.format(chrom=chrom, start=start, end=end)
    with urllib.request.urlopen(url, timeout=30) as resp:
        payload = json.load(resp)
    return payload["dna"].upper()


def fetch_local(fasta, chrom, start, end):
    contig = matching_contig(fasta.references, chrom)
    if contig is None:
        sys.exit(f"{chrom} not found in FASTA (available e.g. {list(fasta.references)[:3]}...)")
    return fasta.fetch(contig, start, end).upper()      # 0-based half-open, matches BED


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-fasta", help="local GRCh38 FASTA (with .fai); omit to use UCSC")
    args = ap.parse_args()

    intervals = unique_intervals(BED)
    OUT.parent.mkdir(parents=True, exist_ok=True)

    fasta = None
    if args.ref_fasta:
        import pysam
        fasta = pysam.FastaFile(args.ref_fasta)          # builds/uses the .fai index

    records = []
    for chrom, start, end in intervals:
        dna = fetch_local(fasta, chrom, start, end) if fasta else fetch_ucsc(chrom, start, end)
        if len(dna) != end - start:
            sys.exit(f"{chrom}:{start}-{end} expected {end - start} bp, got {len(dna)}")
        records.append(f">{chrom}:{start}-{end}\n{dna}")
        print(f"  {chrom}:{start}-{end}  {len(dna)} bp", file=sys.stderr)
        if not fasta:
            time.sleep(0.2)                              # be polite to the public API

    OUT.write_text("\n".join(records) + "\n")
    print(f"Wrote {len(records)} intervals to {OUT.relative_to(REPO)}"
          f"  (source: {'FASTA ' + args.ref_fasta if fasta else 'UCSC'})", file=sys.stderr)


if __name__ == "__main__":
    main()
