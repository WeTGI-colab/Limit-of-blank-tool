"""Strand-resolved per-base pile-up over amplicon intervals.

For every reference position in the requested intervals the four nucleotide counts are
tallied separately for the forward and reverse read strands. Strand resolution is essential:
a genuine low-VAF variant is supported on both strands, whereas a systematic sequencing
artefact is typically strand-biased, and that asymmetry is what later separates the two.

Read-level quality filtering mirrors routine clinical practice (minimum base quality, minimum
mapping quality, and exclusion of secondary/supplementary/duplicate/QC-fail alignments).
"""
import pysam

from amplicon_lob.samples import matching_contig

BASES = "ACGT"


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
                                  ignore_orphans=False, min_base_quality=min_bq):
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
                        counts[base][1 if aln.is_reverse else 0] += 1
                result[(chrom, col.reference_pos)] = counts
    return result
