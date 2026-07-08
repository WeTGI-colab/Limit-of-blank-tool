"""Strand-resolved per-base pile-up over amplicon intervals.

For every reference position the four nucleotide counts are tallied separately for the two
amplicon directions. In this assay the reads are single-end and the sequencing direction is
NOT the BAM reverse flag (which is relative to the genome); it is the amplicon that produced
the read, recorded in the ``CO:Z:Amplicon: <name>_F_/_R_`` tag -- ``_F_`` = forward amplicon,
``_R_`` = reverse amplicon. Strand resolution matters: a genuine variant is supported by both
directions, whereas a systematic sequencing artefact is typically direction-biased.

Read-level quality filtering mirrors routine clinical practice (minimum base quality, minimum
mapping quality, and exclusion of secondary/supplementary/duplicate/QC-fail alignments).
"""
import pysam

from amplicon_lob.samples import matching_contig

BASES = "ACGT"


def read_is_reverse(aln):
    """True for a reverse-amplicon read.

    Uses the amplicon direction from the ``CO`` tag (``..._R_`` = reverse, ``..._F_`` = forward)
    when present; otherwise falls back to the BAM reverse flag (e.g. synthetic test data).
    """
    if aln.has_tag("CO"):
        return str(aln.get_tag("CO")).strip().endswith("_R_")
    return aln.is_reverse


def pileup_bam(bam_path, regions, min_bq=20, min_mq=20):
    """Return {(chrom, pos0): {base: [forward_count, reverse_count]}} over the intervals.

    Results are keyed by the panel chromosome name (chr-prefixed) regardless of how the BAM
    names its contigs; an unindexed BAM is indexed on the fly.
    """
    result = {}
    with pysam.AlignmentFile(str(bam_path)) as bam:
        if not bam.has_index():
            pysam.index(str(bam_path))
            bam = pysam.AlignmentFile(str(bam_path))
        for chrom, start, end in regions:
            contig = matching_contig(bam.references, chrom)
            if contig is None:
                continue
            for col in bam.pileup(contig, start, end, truncate=True, stepper="all",
                                  ignore_orphans=False, min_base_quality=min_bq,
                                  max_depth=100_000_000):
                counts = {b: [0, 0] for b in BASES}
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip:
                        continue
                    aln = pr.alignment
                    if (aln.is_secondary or aln.is_supplementary or aln.is_duplicate
                            or aln.is_qcfail or aln.mapping_quality < min_mq):
                        continue
                    base = aln.query_sequence[pr.query_position]
                    if base in counts:
                        counts[base][1 if read_is_reverse(aln) else 0] += 1
                result[(chrom, col.reference_pos)] = counts
    return result
