"""Panel definitions: amplicon intervals and reference bases (genomic coordinates)."""
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BED = REPO / "bed" / "panel.bed"
REGIONS = REPO / "reference" / "regions.fa"


def load_amplicons(bed=BED):
    """Unique amplicons (forward/reverse primer rows collapsed) with insert (thick) coords."""
    amps = {}
    for line in bed.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        name = f[3].rstrip("_").rsplit("_", 1)[0]        # 'TP53_009_F_' -> 'TP53_009'
        amps[name] = {
            "name": name, "gene": name.split("_")[0], "chrom": f[0],
            "start": int(f[1]), "end": int(f[2]),
            "thick_start": int(f[6]), "thick_end": int(f[7]),
        }
    return list(amps.values())


def thick_regions(amps):
    """Insert intervals (chrom, start0, end0) -- primer flanks excluded."""
    return [(a["chrom"], a["thick_start"], a["thick_end"]) for a in amps]


def gene_by_chrom(amps):
    return {a["chrom"]: a["gene"] for a in amps}


def load_reference_bases(regions=REGIONS):
    """Return {(chrom, pos1): reference_base} from the compact region FASTA.

    Record headers are 'chrom:start-end' (0-based half-open, matching the BED).
    """
    bases, header = {}, None
    chrom = start = None
    for line in regions.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:].strip()
            loc, coords = header.split(":")
            chrom = loc
            start = int(coords.split("-")[0])
            seq_pos = 0
        elif header:
            for i, b in enumerate(line.strip()):
                bases[(chrom, start + seq_pos + i + 1)] = b   # 1-based genomic
            seq_pos += len(line.strip())
    return bases
