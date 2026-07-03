"""Sample manifest handling and chromosome-name normalisation.

Production BAM/VCF live outside the repository; a manifest (data/manifest.tsv) lists, one row
per sample, the sample id, its run id, and absolute paths to the BAM and paired VCF. Both the
synthetic generator and the production discovery script emit this manifest, so the rest of the
pipeline is agnostic to where the reads came from.

Chromosome naming may differ between the panel BED ('chr17') and the sequencer output ('17');
helpers here reconcile the two.
"""
from pathlib import Path


def load_manifest(path):
    """Return a list of {sample, run, bam, vcf} dicts."""
    rows = []
    for line in Path(path).read_text().splitlines()[1:]:
        if not line.strip():
            continue
        f = line.split("\t")
        rows.append({"sample": f[0], "run": f[1], "bam": f[2],
                     "vcf": f[3] if len(f) > 3 else ""})
    return rows


def matching_contig(references, chrom):
    """Map a panel chromosome ('chr17') to the contig name present in a BAM ('chr17' or '17')."""
    refs = set(references)
    if chrom in refs:
        return chrom
    plain = chrom[3:] if chrom.startswith("chr") else chrom
    for cand in (plain, "chr" + plain):
        if cand in refs:
            return cand
    return None


def normalize_chrom(chrom):
    """Force the panel convention (chr-prefixed) on a chromosome name."""
    return chrom if chrom.startswith("chr") else "chr" + chrom
