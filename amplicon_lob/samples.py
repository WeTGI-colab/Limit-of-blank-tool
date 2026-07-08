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


def load_masks(manifest):
    """{sample: set((chrom, pos1, alt))} of genuine (FILTER=PASS) calls to exclude.

    Only FILTER=PASS records are masked -- those are the genuine patient variants. Records the
    caller has flagged (non-PASS, e.g. Pisces ``high_diff_MBQ`` oxidation artefacts) are the
    background we are modelling, so they are deliberately KEPT in the blank. Multi-allelic ALTs
    are split; chromosome names are normalised to the panel convention.
    """
    masks = {}
    for row in manifest:
        keys = set()
        vcf = row.get("vcf")
        if vcf and Path(vcf).exists():
            for line in Path(vcf).read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split("\t")
                if len(f) < 7 or f[6] != "PASS":
                    continue                       # keep caller-flagged artefacts in the blank
                chrom, pos = normalize_chrom(f[0]), int(f[1])
                for alt in f[4].split(","):
                    keys.add((chrom, pos, alt))
        masks[row["sample"]] = keys
    return masks
