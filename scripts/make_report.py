#!/usr/bin/env python3
"""Build a concise PDF report of the positional Limit-of-Blank artefact analysis.

Reads results/cohort_alt.tsv, results/lob_table.tsv and results/artefact_by_run.tsv, plus the
per-site VCF-status findings (which come from the VCF check, encoded below), and writes
results/artefact_report.pdf -- methods, results, and the sites geneticists should review.
"""
import datetime as dt
import sys
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                                Image, PageBreak)

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob.betabin import betabin_lob_from_moments  # noqa: E402

RESULTS = REPO / "results"
OUT = RESULTS / "artefact_report.pdf"

# Per-site VCF status (from the VCF check across all 287 samples). PASS on the artefact allele
# = a false positive that would be reported. status: danger / filtered / not_called / homopolymer
VCF_STATUS = [
    ("chr17:7676093", "G>T", 55, "danger",     "55 samples call G>T as PASS"),
    ("chr17:7675246", "G>T", 41, "danger",     "41 PASS (+8 filtered high_diff_MBQ)"),
    ("chr17:7676102", "G>A", 22, "danger",     "22 samples call G>A as PASS"),
    ("chr17:7675046", "T>G", 5,  "danger",     "5 PASS (+1 high_diff_MBQ)"),
    ("chr17:7675050", "C>A", 1,  "danger",     "1 PASS"),
    ("chr17:7676098", "A>T", 0,  "filtered",   "91 high_diff_MBQ (artefact, filtered)"),
    ("chr17:7676101", "A>T", 0,  "filtered",   "117 high_diff_MBQ (real insertions PASS separately)"),
    ("chr5:171410504", "C>A", 0, "filtered",   "37 high_diff_MBQ"),
    ("chr5:171410501", "C>A", 0, "filtered",   "27 high_diff_MBQ"),
    ("chr17:7676123", "C>G", 0,  "filtered",   "14 high_diff_MBQ"),
    ("chr5:171410586", "G>T", 0, "filtered",   "6 high_diff_MBQ"),
    ("chr17:7675224", "G>T", 0,  "filtered",   "5 high_diff_MBQ"),
    ("chr17:7673845", "G>A", 0,  "filtered",   "1 high_diff_MBQ"),
    ("chr17:7673662", "C>A", 0,  "not_called", "never called by Pisces"),
    ("chr5:171410523", "C>T", 0, "not_called", "never called by Pisces"),
    ("chr17:7675048", "C>G", 0,  "not_called", "never called by Pisces"),
    ("chr5:171410580", "C>A", 0, "not_called", "never called by Pisces"),
    ("chr17:7676605", "C>A", 0,  "not_called", "never called by Pisces"),
    ("chr17:7676608", "A>T", 0,  "not_called", "never called by Pisces"),
    ("chr5:171410509", "C>T", 1, "homopolymer", "low-complexity/homopolymer; indel artefacts (R5x9)"),
]
STATUS_COLOR = {"danger": colors.HexColor("#c0392b"), "filtered": colors.HexColor("#2e7d32"),
                "not_called": colors.HexColor("#2e7d32"), "homopolymer": colors.HexColor("#b8860b")}


def ensure_betabin(l):
    """Add lob_vaf_betabin / bb_rho if absent, from the stored per-site summary moments.

    run_lob.py writes both columns from the raw per-sample counts. When the report is built from a
    lob_table produced by an older run (or one staged without the per-sample pile-ups), the same
    method-of-moments fit is reconstructed from blank_mean_vaf, blank_sd_vaf and mean_depth -- the
    identical estimator, fed the stored moments instead of the raw counts.
    """
    if "lob_vaf_betabin" in l.columns:
        return l
    # mean per-sample depth: use the column if present, else derive from pooled_n / n_blank
    depth = l["mean_depth"] if "mean_depth" in l.columns else l["pooled_n"] / l["n_blank"]
    # centre on the depth-weighted pooled rate, matching run_lob's per-sample betabin_lob
    m = l["pooled_rate"] if "pooled_rate" in l.columns else l["blank_mean_vaf"]
    l = l.assign(_depth=depth, _m=m)
    out = l.apply(lambda r: betabin_lob_from_moments(
        r["_m"], r["blank_sd_vaf"] ** 2,
        1.0 / r["_depth"] if r["_depth"] else 0.0), axis=1)
    l["lob_vaf_betabin"] = [o[0] for o in out]
    l["bb_rho"] = [o[1] for o in out]
    return l


def stats():
    c = pd.read_csv(RESULTS / "cohort_alt.tsv", sep="\t")
    l = ensure_betabin(pd.read_csv(RESULTS / "lob_table.tsv", sep="\t"))
    pos = c.drop_duplicates("pos")
    # depth by amplicon direction: use the fwd/rev split if present, else reconstruct it from the
    # total depth and the forward strand fraction (older cohort_alt schema).
    if "mean_depth_fwd_filt" in pos.columns:
        dfwd, drev = pos["mean_depth_fwd_filt"], pos["mean_depth_rev_filt"]
    else:
        frac = pos["strand_frac_fwd_filt"].fillna(0.5)
        dfwd, drev = pos["mean_depth_filt"] * frac, pos["mean_depth_filt"] * (1 - frac)
    depth = (dfwd + drev).median()
    v = l["blank_mean_vaf"] * 100
    vclean = l[~l["flagged"]]["blank_mean_vaf"] * 100
    sub = (l[~l["flagged"]].assign(p=lambda d: d["blank_mean_vaf"] * 100)
           .groupby("sub")["p"].mean().sort_values(ascending=False))
    return {
        "n": int(c["n_samples_filt"].max()), "depth": depth,
        "fwd": dfwd.median(), "rev": drev.median(),
        "balanced": bool(min(dfwd.median(), drev.median())
                         / max(dfwd.median(), drev.median(), 1) >= 0.6),
        "npos": len(pos), "nflag": int(l["flagged"].sum()),
        "noise_mean": v.mean(), "noise_sd": v.std(), "noise_med": v.median(),
        "noise_p95": vclean.quantile(0.95), "sub": sub, "lob": l,
    }


def main():
    S = stats()
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle("H", parent=styles["Heading2"], textColor=colors.HexColor("#1a3c5a"),
                              spaceBefore=10, spaceAfter=4))
    body = ParagraphStyle("B", parent=styles["BodyText"], fontSize=9.5, leading=13)
    small = ParagraphStyle("S", parent=body, fontSize=8, textColor=colors.HexColor("#555555"))
    title = ParagraphStyle("T", parent=styles["Title"], fontSize=17, textColor=colors.HexColor("#1a3c5a"))

    E = []
    E.append(Paragraph("Position-specific sequencing artefacts in the amplicon assay", title))
    E.append(Paragraph(f"Positional Limit-of-Blank (CLSI EP17) analysis &mdash; internal QC / validation. "
                       f"Cohort: {S['n']} patient samples · AML panel (TP53, NPM1, FLT3, IDH1, IDH2) · "
                       f"GRCh38 · {dt.date.today().isoformat()}.", small))
    E.append(Spacer(1, 6))

    E.append(Paragraph("1. Purpose", styles["H"]))
    E.append(Paragraph(
        "We report somatic mutations at very low VAF (~1&ndash;2%). The assay produces "
        "<b>systematic</b> sequencing artefacts (not random error) at specific positions, base "
        "substitutions and amplicon directions. These sit in the same 1&ndash;2% range as genuine "
        "low-level mutations and can be mistaken for them. This analysis characterises that "
        "background so that a per-position noise ceiling can distinguish artefact from signal.", body))

    E.append(Paragraph("2. Methods", styles["H"]))
    E.append(Paragraph(
        f"Across {S['n']} patient BAMs over the fixed amplicon panel: "
        "(1) at every panel position we count non-reference reads per <b>amplicon direction "
        "(forward/reverse, from the read amplicon tag)</b> and per substitution, under two QC "
        "regimes &mdash; <i>raw</i> (minbq1/minmq1) and the caller's filter <i>(minbq20/minmq30)</i>. "
        "(2) From each sample's VCF we remove the genuine <b>FILTER=PASS</b> variants but keep the "
        "caller-flagged <b>non-PASS</b> records &mdash; those are the artefacts we model; what remains "
        "is the &lsquo;blank&rsquo;. (3) We aggregate across all samples, flag positions whose blank "
        "sits systematically above the random floor, and express each as a <b>positional Limit of "
        "Blank</b> (mean + 1.645&middot;SD, CLSI EP17). (4) We also break the signal down per "
        "sequencing run. This is not a new method (cf. AmpliSolve, iDES); the contribution is "
        "assay-specific validation on our own data.", body))

    E.append(Paragraph("3. Results", styles["H"]))
    bal = ("the two directions are well balanced, so strand comparison is valid across the panel"
           if S["balanced"] else
           "the forward/reverse coverage is uneven, so strand comparisons should be read with care")
    E.append(Paragraph(
        f"<b>Coverage.</b> Median depth &asymp; {S['depth']:,.0f}&times; per position "
        f"(forward amplicon {S['fwd']:,.0f}&times;, reverse {S['rev']:,.0f}&times; &mdash; {bal}).", body))
    E.append(Paragraph(
        f"<b>Background noise floor.</b> Across {S['npos']:,} positions the blank non-reference VAF "
        f"averages <b>{S['noise_mean']:.3f}% &plusmn; {S['noise_sd']:.3f}%</b> (median "
        f"{S['noise_med']:.3f}%; 95th percentile {S['noise_p95']:.3f}%). The large SD is not uniform "
        "scatter &mdash; it is driven by a small number of <b>systematic artefact sites</b> that rise far "
        "above the ~0.01% floor (up to ~7.5%). In other words, the noise is position-specific, which "
        "is exactly why a single flat VAF threshold is inadequate.", body))
    E.append(Paragraph(f"Typical floor by substitution (clean sites): " +
                       ", ".join(f"{k} {v:.3f}%" for k, v in S["sub"].head(6).items()) + ".", small))

    # flagged artefact table
    E.append(Paragraph(f"<b>Systematic artefact sites.</b> {S['nflag']} positions were flagged as "
                       "systematically noisy. Almost all are strongly biased to one amplicon "
                       "direction (an artefact of that primer/strand), confirming they are technical, "
                       "not biological. The strongest:", body))
    f = S["lob"][S["lob"]["flagged"]].sort_values("blank_mean_vaf", ascending=False).head(10)
    rows = [["Position", "Sub", "Blank VAF", "Strand (fwd frac)", "Positional LoB"]]
    for _, r in f.iterrows():
        rows.append([f"{r['chrom']}:{r['pos']}", r["sub"], f"{r['blank_mean_vaf']*100:.3f}%",
                     f"{r['strand_frac_fwd']:.2f}", f"{r['lob_vaf']*100:.3f}%"])
    t = Table(rows, hAlign="LEFT", colWidths=[3.2*cm, 1.4*cm, 2.2*cm, 3.0*cm, 2.6*cm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1a3c5a")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f2f6fa")])]))
    E.append(t)
    E.append(Paragraph("Strand fraction ~1 = forward-amplicon artefact, ~0 = reverse-amplicon; "
                       "~0.5 (both) would suggest a genuine variant.", small))

    E.append(Paragraph(
        "<b>Batch (per-run) effects.</b> Some artefacts are position-level (elevated across most "
        "runs), while others spike in a single sequencing run. One run in particular (HCH3FF) shows a "
        "cohort-wide G&gt;T/A&gt;T (oxidation) elevation across many positions &mdash; a run/batch QC "
        "failure. The heatmap below reads on two axes: a bright <b>row</b> = a toxic position (all "
        "runs); a bright <b>column</b> = a bad run.", body))
    if (RESULTS / "plots" / "artefact_by_run.png").exists():
        E.append(Image(str(RESULTS / "plots" / "artefact_by_run.png"), width=17*cm, height=11*cm))

    E.append(PageBreak())
    E.append(Paragraph("4. What geneticists should review", styles["H"]))
    E.append(Paragraph(
        "Pisces already filters the low-base-quality artefacts (FILTER <i>high_diff_MBQ</i>) or does "
        "not call them &mdash; those are safe. The concern is the <b>high-base-quality artefacts that "
        "Pisces passes as FILTER=PASS</b>: these would be <b>reported as genuine low-VAF mutations</b>. "
        "The check below (all 287 VCFs) found artefacts reported as PASS at the following TP53 "
        "positions:", body))
    dr = [x for x in VCF_STATUS if x[3] == "danger"]
    rows = [["Position", "Sub", "PASS calls", "Note"]]
    for site, sub, npass, _, note in dr:
        rows.append([site, sub, str(npass), note])
    t = Table(rows, hAlign="LEFT", colWidths=[3.2*cm, 1.4*cm, 2.2*cm, 8.0*cm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#c0392b")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#fbeeec")])]))
    E.append(t)
    E.append(Paragraph(
        "<b>Action:</b> these are the positions where the positional Limit of Blank adds value &mdash; a "
        "candidate call here should exceed the site-specific LoB, not the flat 1% threshold. "
        "Geneticists should review whether any issued reports contain low-VAF calls at "
        "<b>chr17:7676093 G&gt;T, chr17:7675246 G&gt;T or chr17:7676102 G&gt;A</b>, and treat samples "
        "from run <b>HCH3FF</b> with particular caution.", body))

    E.append(Paragraph("Full VCF status of the 20 top artefact positions", styles["H"]))
    rows = [["Position", "Sub", "PASS", "VCF status", "Detail"]]
    for site, sub, npass, status, note in VCF_STATUS:
        rows.append([site, sub, str(npass), status.replace("_", " "), note])
    t = Table(rows, hAlign="LEFT", colWidths=[3.0*cm, 1.3*cm, 1.2*cm, 2.4*cm, 8.1*cm])
    style = [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1a3c5a")),
             ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTSIZE", (0, 0), (-1, -1), 7.5),
             ("GRID", (0, 0), (-1, -1), 0.3, colors.grey)]
    for i, (_, _, _, status, _) in enumerate(VCF_STATUS, start=1):
        style.append(("TEXTCOLOR", (3, i), (3, i), STATUS_COLOR[status]))
    t.setStyle(TableStyle(style))
    E.append(t)
    E.append(Paragraph("PASS = calls of the artefact allele reported as genuine variants. "
                       "green = filtered/not-called (safe); red = reported (danger); "
                       "amber = homopolymer/indel region.", small))

    E.append(Paragraph("5. Strand diagnostic", styles["H"]))
    E.append(Paragraph("Each point is a position&times;alternate-base in the blank. As the noise VAF "
                       "rises, points move to strand-fraction 0 or 1 (one amplicon direction), not to "
                       "0.5 &mdash; confirming the elevated blank is technical artefact, not balanced "
                       "biological signal.", body))
    if (RESULTS / "plots" / "strand_vs_vaf.png").exists():
        E.append(Image(str(RESULTS / "plots" / "strand_vs_vaf.png"), width=11*cm, height=7.3*cm))

    E.append(Paragraph("6. Interpretation and recommendation", styles["H"]))
    E.append(Paragraph(
        "<b>Our reading of the data.</b> The assay shows two distinct, separable failure modes. "
        "(a) <b>Position-specific artefacts</b> &mdash; a handful of loci are intrinsically noisy in "
        "essentially every run (bright rows in the heatmap), dominated by oxidation/deamination "
        "substitutions (G&gt;T, C&gt;T, A&gt;T) and biased to a single amplicon direction. "
        "(b) <b>Run/batch effects</b> &mdash; an entire flow-cell can be elevated across many positions "
        "(a bright column; run HCH3FF), consistent with oxidative DNA damage in that library batch. "
        "The strand diagnostic confirms both are technical: elevated blank signal is strand-biased "
        "(fraction &rarr; 0 or 1), never balanced like a genuine variant would be.", body))
    E.append(Paragraph(
        "<b>What matters clinically.</b> The artefacts split by how the caller already handles them. "
        "Low-base-quality artefacts are filtered by Pisces (<i>high_diff_MBQ</i>) or never called "
        "&mdash; these are safe. A minority are <b>high base quality</b>: they survive the quality "
        "filter and are emitted as <b>FILTER=PASS</b>, so they would be reported as real low-VAF "
        "mutations. Those (chr17:7676093 G&gt;T, 7675246 G&gt;T, 7676102 G&gt;A) are the genuine "
        "residual risk, and are precisely why a positional Limit of Blank is needed <i>on top of</i> "
        "the existing quality filters.", body))
    E.append(Paragraph(
        "<b>Recommendation.</b> (1) Apply the positional LoB at the flagged sites &mdash; a candidate "
        "call must exceed the site-specific ceiling, not a flat 1% threshold. (2) Add a per-run "
        "artefact QC gate to catch damaged batches such as HCH3FF before sign-out. (3) Handle the "
        "low-complexity/homopolymer site (chr5:171410509) separately: its artefacts are indel "
        "slippage, not a substitution, and need a different rule.", body))

    E.append(Paragraph("7. Planned improvements from the literature", styles["H"]))
    E.append(Paragraph(
        "The model used here is deliberately simple &mdash; a Gaussian per-position ceiling "
        "(mean + 1.645&middot;SD). The published prior art (all built on the same idea: a background "
        "model from normals, applied per position) offers concrete refinements. None change the "
        "design; they improve the <i>statistics</i> and the <i>deployment</i>. In priority order:", body))
    rows = [["Source", "What it adds beyond our current model", "Status here"]]
    imp = [
        ("Beta-binomial site model\n(deepSNV / shearwater)",
         "Replace the Gaussian mean+1.645·SD ceiling with a beta-binomial fitted per site: it "
         "models between-sample overdispersion instead of one pooled rate, so the LoB is properly "
         "calibrated in the low-count tail where a Gaussian under- or over-shoots.",
         "Planned — main\nstatistical upgrade"),
        ("iDES / CAPP-Seq\n‘polishing’ (Newman 2016)",
         "Deploy the LoB as in-silico polishing: subtract each position’s own background before "
         "thresholding a candidate, rather than comparing against a flat 1% cut-off.",
         "Planned —\ndeployment step"),
        ("TNER (2018)",
         "Tri-nucleotide-context error rates: pool sparse positions by their 3-mer context so "
         "low-coverage sites borrow statistical strength from similar contexts.",
         "Planned —\nrefinement"),
        ("Mutect2 PoN +\nread-orientation model",
         "Formalise the strand bias we already observe into an orientation-bias filter (oxidation "
         "G>T, deamination C>T), and version the flagged sites as a reusable panel-of-normals "
         "blacklist shipped with the pipeline.",
         "Partly — diagnostic\ndone, filter planned"),
        ("AmpliSolve (2019)",
         "Closest precedent: position-, strand- and nucleotide-specific background from a panel of "
         "normals, calling SNVs to ~1% VAF. Our per-position / per-strand / per-substitution blank "
         "is the same construction — confirms the design.",
         "Implemented\n(this analysis)"),
    ]
    cell = ParagraphStyle("cell", parent=body, fontSize=7.5, leading=9.5)
    for src, add, st in imp:
        rows.append([Paragraph(src.replace("\n", "<br/>"), cell),
                     Paragraph(add, cell), Paragraph(st.replace("\n", "<br/>"), cell)])
    t = Table(rows, hAlign="LEFT", colWidths=[3.6*cm, 8.6*cm, 3.6*cm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1a3c5a")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTSIZE", (0, 0), (-1, 0), 8),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.grey), ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f2f6fa")])]))
    E.append(t)
    E.append(Paragraph(
        "<b>Bottom line:</b> the highest-value next step is the beta-binomial per-site LoB &mdash; it "
        "directly hardens the noise ceiling that everything else depends on. The rest (context "
        "pooling, orientation filter, in-silico polishing) are incremental refinements layered on "
        "the same per-position background.", small))

    E.append(Paragraph("8. Two LoB estimators: Gaussian vs beta-binomial", styles["H"]))
    E.append(Paragraph(
        "To make the first of those improvements concrete, both ceilings are now computed for "
        "every site. The <b>Gaussian</b> LoB is mean + 1.645&middot;SD of the per-sample blank "
        "fractions. The <b>beta-binomial</b> LoB fits a Beta distribution to the per-sample "
        "non-reference counts by method of moments &mdash; separating ordinary binomial sampling "
        "from genuine between-sample spread (<b>&rho;</b>, the overdispersion: 0 = pure sampling, "
        "&rarr;1 = strong run/batch variability) &mdash; and takes the 95th percentile of that "
        "fitted rate. The flagged set is unchanged (flagging uses the fold-change / t-test, not the "
        "ceiling); what changes is the <i>height</i> of each site&rsquo;s ceiling.", body))
    L = S["lob"]
    fl = L[L["flagged"]].sort_values("blank_mean_vaf", ascending=False).head(12)
    rows = [["Position", "Sub", "Blank VAF", "LoB Gaussian", "LoB beta-binomial", "Overdisp. ρ"]]
    tighter = 0
    for _, r in fl.iterrows():
        bb, ga = r["lob_vaf_betabin"] * 100, r["lob_vaf"] * 100
        if bb < ga:
            tighter += 1
        rows.append([f"{r['chrom']}:{r['pos']}", r["sub"], f"{r['blank_mean_vaf']*100:.3f}%",
                     f"{ga:.3f}%", f"{bb:.3f}%",
                     f"{r['bb_rho']:.3f}" if pd.notna(r["bb_rho"]) else "-"])
    t = Table(rows, hAlign="LEFT",
              colWidths=[3.0*cm, 1.2*cm, 2.2*cm, 2.8*cm, 3.4*cm, 2.2*cm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1a3c5a")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f2f6fa")])]))
    E.append(t)
    med_rho = L[L["flagged"]]["bb_rho"].median()
    E.append(Paragraph(
        f"<b>How to read it.</b> At sites with little overdispersion (&rho; near 0) the beta-binomial "
        f"ceiling sits close to the mean blank rate and is <b>tighter</b> than the Gaussian, which "
        f"inflates the ceiling through an SD that also absorbs sampling noise &mdash; a tighter ceiling "
        f"means better sensitivity to genuine low-VAF calls. At overdispersed sites (higher &rho;, "
        f"typically the run/batch-driven ones) the beta-binomial instead <b>widens</b> the ceiling to "
        f"the skewed upper tail, guarding against false positives where the Gaussian&rsquo;s symmetric "
        f"SD under-covers it. Across the {int(L['flagged'].sum())} flagged sites the beta-binomial is "
        f"tighter at {tighter} of the top {len(fl)} shown (median &rho; &asymp; {med_rho:.3f}). "
        f"Recommendation: adopt the beta-binomial ceiling as the operational LoB and keep the Gaussian "
        f"column as a transparent reference.", small))

    E.append(Spacer(1, 6))
    E.append(Paragraph("Prepared for review. Method and code are version-controlled; figures and "
                       "tables regenerate from the current cohort. Not for external distribution.",
                       small))

    SimpleDocTemplate(str(OUT), pagesize=A4, leftMargin=1.8*cm, rightMargin=1.8*cm,
                      topMargin=1.6*cm, bottomMargin=1.6*cm).build(E)
    print(f"Wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
