#!/usr/bin/env python3
"""Generate a synthetic, ground-truth cohort of amplicon BAM/VCF files over the panel BED.

Design rationale
----------------
The objective is to develop and *validate* a positional background-error (Limit of Blank)
model. Validation requires known truth, which real data cannot provide. Each sample is
therefore simulated over the exact panel intervals (GRCh38) with three superimposed signals:

  1. Baseline stochastic error   -- a low, uniform per-base substitution rate (the random
                                     "floor" of the instrument).
  2. Systematic artefact sites   -- specific (position, alternate base) loci with an elevated,
                                     strand-biased error rate. These emulate the reproducible
                                     miscalls the assay is known to produce and are the events
                                     the detector must flag.
  3. Genuine low-VAF variants    -- true mutations carried by a subset of samples at ~1.5-2%
                                     VAF with balanced strand support, recorded in each
                                     carrier's VCF. After VCF masking these must NOT be flagged.

Reads are stored on the forward reference strand (BAM convention); read orientation is encoded
in the FLAG (0x10) so that strand-biased artefacts produce an asymmetric forward/reverse alt
count, exactly as in real data.

Outputs (under data/synthetic/):
  sample_XX/sample_XX.bam(.bai)   aligned reads at true GRCh38 coordinates
  sample_XX/sample_XX.vcf         genuine variants carried by that sample
  ground_truth.tsv                the full truth set (artefacts + real variants)
"""
import argparse
import array
import os
import numpy as np
import pysam
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BED = REPO / "bed" / "panel.bed"
REF = REPO / "reference" / "regions.fa"
OUTDIR = REPO / "data" / "synthetic"

# True GRCh38 chromosome lengths for the panel contigs (BAM @SQ header).
CHROM_LEN = {
    "chr2": 242193529, "chr5": 181538259, "chr13": 114364328,
    "chr15": 101991189, "chr17": 83257441,
}

BASES = np.array(list("ACGT"))
IDX = {b: i for i, b in enumerate(BASES)}
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}  # realistic dominant substitution

BASELINE_ERR = 5e-4

# Ground truth, expressed relative to each amplicon's thick (insert) start so positions are
# always valid. alt = transition of the reference base at that locus (resolved at runtime).
ARTEFACT_SITES = [
    {"amplicon": "TP53_009", "offset": 40, "strand": "+", "rate": 0.015},
    {"amplicon": "FLT3_004", "offset": 10, "strand": "-", "rate": 0.020},
    {"amplicon": "TP53_005", "offset": 100, "strand": "+", "rate": 0.012},
    {"amplicon": "NPM1_002", "offset": 50, "strand": "-", "rate": 0.010},
]
REAL_VARIANTS = [
    {"amplicon": "IDH1_001", "offset": 26, "vaf": 0.020, "samples": [0, 1, 2, 3, 4]},
    {"amplicon": "TP53_007", "offset": 80, "vaf": 0.015, "samples": [2, 3, 6, 7]},
    {"amplicon": "IDH2_003", "offset": 20, "vaf": 0.018, "samples": [0, 5, 8]},
]


def load_amplicons():
    """Return {amplicon_base_name: dict(chrom, start, end, thick_start, thick_end)}."""
    amps = {}
    for line in BED.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        chrom, start, end, name = f[0], int(f[1]), int(f[2]), f[3]
        base = name.rstrip("_").rsplit("_", 1)[0]  # strip trailing F_/R_ -> e.g. TP53_009
        amps[base] = {
            "chrom": chrom, "start": start, "end": end,
            "thick_start": int(f[6]), "thick_end": int(f[7]),
        }
    return amps


def load_reference():
    """Return {'chrom:start-end': sequence} from the compact region FASTA."""
    seqs, header = {}, None
    for line in REF.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:].strip()
            seqs[header] = []
        elif header:
            seqs[header].append(line.strip())
    return {h: "".join(s) for h, s in seqs.items()}


def amplicon_seq(amp, refseqs):
    """Full amplicon reference sequence (forward strand, includes primer flanks)."""
    return refseqs[f"{amp['chrom']}:{amp['start']}-{amp['end']}"]


def resolve_truth(amps, refseqs):
    """Attach genomic position, reference and alternate base to each truth entry."""
    def resolve(entry):
        amp = amps[entry["amplicon"]]
        seq = amplicon_seq(amp, refseqs)
        gpos = amp["thick_start"] + entry["offset"]          # genomic, 0-based
        ref = seq[gpos - amp["start"]]                        # base into the amplicon record
        return {**entry, "chrom": amp["chrom"], "pos": gpos,
                "ref": ref, "alt": TRANSITION[ref]}
    return [resolve(e) for e in ARTEFACT_SITES], [resolve(e) for e in REAL_VARIANTS]


def build_header():
    return {"HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": c, "LN": CHROM_LEN[c]} for c in
                   sorted(CHROM_LEN, key=lambda x: int(x[3:]))]}


def simulate_amplicon(amp, seq, depth, artefacts, real_here, rng):
    """Return (read_matrix[depth, L] of base indices, strands[depth]) for one amplicon.

    Reads span the thick (insert) interval only; primer flanks are not sequenced.
    """
    ts, te = amp["thick_start"], amp["thick_end"]
    L = te - ts
    template = np.array([IDX[b] for b in seq[ts - amp["start"]: te - amp["start"]]])
    reads = np.tile(template, (depth, 1))

    # strand: first half forward (flag 0), second half reverse (flag 16), then shuffle
    strands = np.zeros(depth, dtype=int)
    strands[depth // 2:] = 16
    rng.shuffle(strands)

    # Per-base quality (Phred) and per-read mapping quality. True (reference/genuine) bases are
    # high quality; artefacts and random errors carry lower base quality, so the two QC regimes
    # (raw vs Pisces-like minbq20/minmq30) see different amounts of signal.
    bq = np.full((depth, L), 37, dtype=np.uint8)

    # 1) baseline stochastic error: random substitution, low base quality
    err = rng.random((depth, L)) < BASELINE_ERR
    if err.any():
        shift = rng.integers(1, 4, size=err.sum())
        reads[err] = (reads[err] + shift) % 4
        bq[err] = rng.integers(2, 16, size=err.sum())          # errors are low-quality

    # 2) systematic strand-biased artefacts, base quality straddling the minbq=20 cutoff
    for a in artefacts:
        if a["chrom"] != amp["chrom"] or not (ts <= a["pos"] < te):
            continue
        col = a["pos"] - ts
        on_strand = strands == (16 if a["strand"] == "-" else 0)
        hit = on_strand & (rng.random(depth) < a["rate"])
        reads[hit, col] = IDX[a["alt"]]
        bq[hit, col] = rng.integers(8, 32, size=int(hit.sum()))  # ~half below minbq=20

    # 3) genuine variants (balanced strand, high base quality) carried by this sample
    for v in real_here:
        if v["chrom"] != amp["chrom"] or not (ts <= v["pos"] < te):
            continue
        col = v["pos"] - ts
        hit = rng.random(depth) < v["vaf"]
        reads[hit, col] = IDX[v["alt"]]
        bq[hit, col] = 37                                       # real variants pass the filter

    # mapping quality: mostly high, a minority below the minmq=30 filter
    mq = np.full(depth, 60, dtype=int)
    n_low = int(0.08 * depth)
    if n_low:
        mq[rng.choice(depth, n_low, replace=False)] = 20

    return reads, strands, bq, mq


def write_sample(idx, amps, refseqs, artefacts, real_variants, depth, tmpdir):
    sname = f"sample_{idx:02d}"
    sdir = OUTDIR / sname
    sdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1000 + idx)
    real_here = [v for v in real_variants if idx in v["samples"]]

    header = build_header()
    unsorted = tmpdir / f"{sname}.unsorted.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
        rid = {c: i for i, c in enumerate(
            sorted(CHROM_LEN, key=lambda x: int(x[3:])))}
        for base, amp in amps.items():
            seq = amplicon_seq(amp, refseqs)
            reads, strands, bq, mq = simulate_amplicon(
                amp, seq, depth, artefacts, real_here, rng)
            L = reads.shape[1]
            for r in range(depth):
                seg = pysam.AlignedSegment()
                seg.query_name = f"{base}:{r}"
                seg.flag = int(strands[r])
                seg.reference_id = rid[amp["chrom"]]
                seg.reference_start = amp["thick_start"]
                seg.mapping_quality = int(mq[r])
                seg.cigar = [(0, L)]
                seg.query_sequence = "".join(BASES[reads[r]])
                seg.query_qualities = array.array("B", bq[r].tolist())
                bam.write(seg)
    pysam.sort("-o", str(sdir / f"{sname}.bam"), str(unsorted))
    pysam.index(str(sdir / f"{sname}.bam"))
    os.remove(unsorted)

    # per-sample VCF of genuine variants
    with open(sdir / f"{sname}.vcf", "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="Designed variant allele frequency">\n')
        for c in sorted(CHROM_LEN, key=lambda x: int(x[3:])):
            vcf.write(f"##contig=<ID={c},length={CHROM_LEN[c]}>\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in real_here:
            vcf.write(f"{v['chrom']}\t{v['pos'] + 1}\t.\t{v['ref']}\t{v['alt']}"
                      f"\t.\tPASS\tVAF={v['vaf']}\n")
    return sname


def write_ground_truth(artefacts, real_variants):
    path = OUTDIR / "ground_truth.tsv"
    with open(path, "w") as fh:
        fh.write("kind\tamplicon\tchrom\tpos_1based\tref\talt\tstrand\trate_or_vaf\tsamples\n")
        for a in artefacts:
            fh.write(f"artefact\t{a['amplicon']}\t{a['chrom']}\t{a['pos'] + 1}\t{a['ref']}"
                     f"\t{a['alt']}\t{a['strand']}\t{a['rate']}\tALL\n")
        for v in real_variants:
            s = ";".join(f"sample_{i:02d}" for i in v["samples"])
            fh.write(f"real\t{v['amplicon']}\t{v['chrom']}\t{v['pos'] + 1}\t{v['ref']}"
                     f"\t{v['alt']}\t.\t{v['vaf']}\t{s}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples", type=int, default=10)
    ap.add_argument("--depth", type=int, default=5000)
    args = ap.parse_args()

    amps = load_amplicons()
    refseqs = load_reference()
    artefacts, real_variants = resolve_truth(amps, refseqs)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_ground_truth(artefacts, real_variants)

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        for i in range(args.samples):
            name = write_sample(i, amps, refseqs, artefacts, real_variants,
                                args.depth, Path(td))
            print(f"  wrote {name}  ({args.depth}x over {len(amps)} amplicons)")

    # manifest so the downstream pipeline treats synthetic and production data identically
    with open(OUTDIR.parent / "manifest.tsv", "w") as fh:
        fh.write("sample\trun\tbam\tvcf\n")
        for i in range(args.samples):
            name = f"sample_{i:02d}"
            fh.write(f"{name}\tsynthetic\t{OUTDIR / name / (name + '.bam')}"
                     f"\t{OUTDIR / name / (name + '.vcf')}\n")
    print(f"Cohort: {args.samples} samples in {OUTDIR.relative_to(REPO)} "
          f"(manifest: data/manifest.tsv)")
    print("Truth set:")
    print((OUTDIR / 'ground_truth.tsv').read_text())


if __name__ == "__main__":
    main()
