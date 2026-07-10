#!/usr/bin/env python3
"""Export every variant called across the cohort's VCFs, for the record.

Walks every VCF listed in data/manifest.tsv and writes two tables:

  results/vcf_variants.tsv          long: one row per (sample, variant) -- sample, run, gene,
                                    locus, ref, alt, FILTER, QUAL, VAF and depth (parsed from the
                                    Pisces FORMAT/sample column: VF/AF and DP).
  results/vcf_variants_summary.tsv  per (chrom, pos, ref, alt): how many samples called it, how
                                    many as PASS vs caller-filtered, the set of FILTERs seen, mean
                                    VAF/depth and which runs carried a PASS call.

This is a plain catalogue of what the caller reported (all FILTERs, PASS and non-PASS). It is
independent of the blank/LoB model -- it does not mask anything -- so it can be cross-checked
against the flagged artefact sites (a site that is a systematic artefact but is repeatedly emitted
as PASS is the clinical risk). Multi-allelic ALTs are split; per-allele FORMAT fields (VF/AD) are
indexed by allele where present.
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import panel, samples  # noqa: E402

MANIFEST = REPO / "data" / "manifest.tsv"
RESULTS = REPO / "results"


def _fmt_field(fmt, sample_col):
    """Return {key: value} for the (first-sample) FORMAT column, or {} if absent."""
    if not fmt or not sample_col:
        return {}
    return dict(zip(fmt.split(":"), sample_col.split(":")))


def _info_field(info):
    """Return {key: value} for the INFO column (flags map to '')."""
    d = {}
    for kv in info.split(";"):
        if not kv or kv == ".":
            continue
        k, _, v = kv.partition("=")
        d[k] = v
    return d


def _per_allele(value, alt_index, n_alt):
    """Pick allele-specific value from a comma list. AD carries ref first (n_alt+1 values)."""
    if value is None:
        return None
    parts = value.split(",")
    if len(parts) == n_alt:                 # one value per ALT (e.g. VF)
        return parts[alt_index]
    if len(parts) == n_alt + 1:             # ref first, then ALTs (e.g. AD)
        return parts[alt_index + 1]
    return value if len(parts) == 1 else None


def parse_vcf(path):
    """Yield dicts for each (record, alt-allele) in a VCF, without pysam (works on any VCF text)."""
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 8:
            continue
        chrom = samples.normalize_chrom(f[0])
        pos, ref, alts, qual, filt = int(f[1]), f[3], f[4].split(","), f[5], f[6]
        info = _info_field(f[7])
        fmt = _fmt_field(f[8] if len(f) > 8 else "", f[9] if len(f) > 9 else "")
        for i, alt in enumerate(alts):
            # variant frequency and depth: prefer the per-sample FORMAT (Pisces VF/AD/DP),
            # fall back to INFO (VAF/AF/DP) for callers that report there.
            vf = _per_allele(fmt.get("VF", fmt.get("AF")), i, len(alts))
            ad = _per_allele(fmt.get("AD"), i, len(alts))
            dp = fmt.get("DP") or info.get("DP")
            vaf = None
            if vf is not None:
                vaf = float(vf)
            elif ad is not None and dp not in (None, "", "0"):
                vaf = int(ad) / int(dp)
            elif info.get("VAF", info.get("AF")) is not None:
                iv = _per_allele(info.get("VAF", info.get("AF")), i, len(alts))
                vaf = float(iv) if iv is not None else None
            yield {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                   "filter": filt, "qual": qual,
                   "vaf": vaf, "depth": int(dp) if dp not in (None, "", ".") else None}


def main():
    manifest = samples.load_manifest(MANIFEST)
    gene_of = panel.gene_by_chrom(panel.load_amplicons())

    rows, missing = [], 0
    for r in manifest:
        vcf = r.get("vcf")
        if not vcf or not Path(vcf).exists():
            missing += 1
            continue
        for v in parse_vcf(vcf):
            v.update({"sample": r["sample"], "run": r["run"],
                      "gene": gene_of.get(v["chrom"], "")})
            rows.append(v)
    if not rows:
        sys.exit("No variants parsed -- check the VCF paths in the manifest.")

    cols = ["sample", "run", "gene", "chrom", "pos", "ref", "alt", "filter", "qual", "vaf", "depth"]
    long = pd.DataFrame(rows)[cols].sort_values(["chrom", "pos", "ref", "alt", "sample"])
    RESULTS.mkdir(parents=True, exist_ok=True)
    long.to_csv(RESULTS / "vcf_variants.tsv", sep="\t", index=False)

    long["is_pass"] = long["filter"] == "PASS"

    def roll(g):
        filters = sorted(set(g["filter"]))
        pass_runs = sorted(set(g.loc[g["is_pass"], "run"]))
        return pd.Series({
            "gene": g["gene"].iloc[0],
            "n_samples": g["sample"].nunique(),
            "n_pass": int(g["is_pass"].sum()),
            "n_nonpass": int((~g["is_pass"]).sum()),
            "filters": ",".join(filters),
            "mean_vaf": g["vaf"].mean(skipna=True),
            "max_vaf": g["vaf"].max(skipna=True),
            "mean_depth": g["depth"].mean(skipna=True),
            "n_runs_pass": len(pass_runs),
        })
    summary = (long.groupby(["chrom", "pos", "ref", "alt"]).apply(roll, include_groups=False)
               .reset_index().sort_values(["n_pass", "n_samples"], ascending=False))
    summary.to_csv(RESULTS / "vcf_variants_summary.tsv", sep="\t", index=False)

    n_vcf = len(manifest) - missing
    print(f"Parsed {n_vcf} VCFs ({missing} missing) -> {len(long)} variant calls, "
          f"{len(summary)} distinct (pos, ref, alt).")
    print(f"  results/vcf_variants.tsv          (per-sample calls)")
    print(f"  results/vcf_variants_summary.tsv  (per-site rollup, PASS vs filtered)")
    top = summary[summary["n_pass"] > 0].head(10)
    if len(top):
        print("\nMost recurrent PASS variants:")
        with pd.option_context("display.width", 200, "display.max_columns", None):
            print(top[["chrom", "pos", "ref", "alt", "gene", "n_samples", "n_pass",
                       "n_nonpass", "mean_vaf"]].to_string(index=False))


if __name__ == "__main__":
    main()
