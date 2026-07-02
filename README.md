# Positional Limit of Blank (LoB) — Amplicon Artefact Error Model

A bioinformatics method to characterise **position-specific sequencing artefacts** in an
amplicon-based clinical assay, so that genuine low-frequency mutations (VAF ~1–2% at ~5000×)
can be distinguished from systematic background noise.

> **Status:** design / development. Built and validated locally on synthetic data (personal
> machine), then ported to the laboratory environment for use on real assay data.
> **Not intended for external publication** — this is an internal service-improvement and
> validation project.

---

## 1. What we are going to do

In our amplicon test we believe the sequencer produces **systematic errors that are not random** —
they concentrate at particular positions, base substitutions and strands (e.g. strand-biased
miscalls, homopolymer slippage, oxidation/deamination artefacts). Because we call mutations
at very low VAF (1–2%), these artefacts sit in exactly the same frequency range as the real
variants we care about, inflating the apparent VAF and causing false positives.

The plan:

1. Take **many BAM files** from the assay. Amplicons mean the regions are fixed and known
   (defined by a BED file), so every sample contributes coverage at the same loci.
2. Take each sample's **associated VCF** and use it to **remove the real mutations** (and known
   germline variants) from the pileup. What remains at each position is *blank* — background
   error, not biology.
3. Build a **per-position, per-substitution, per-strand background model** of that residual
   error across all samples (beta-binomial to capture between-sample overdispersion).
4. **Flag positions/bases where the artefactual VAF is systematically higher** than the
   random error floor of the instrument — the "noisy sites".
5. Express the result as a **positional Limit of Blank (CLSI EP17)**: instead of one flat VAF
   threshold across the panel, each position gets its own blank ceiling that a candidate
   variant must exceed to be called.

**Deliverables:** pileup → VCF-masking → positional error model → systematic-site detection →
per-position LoB table + validation report, developed against synthetic ground-truth data.

---

## 2. Literature review (prior art)

This approach is **well established** in the literature — the novelty here is not the method
but its **assay-specific validation** and the **positional-LoB (EP17) framing** using our own
clinical VCFs to define the blank.

| Work | One-line explanation | Relevance |
|---|---|---|
| **AmpliSolve** (Bioinformatics, 2019) — [bioRxiv](https://www.biorxiv.org/content/10.1101/475947v2) · [PubMed](https://pubmed.ncbi.nlm.nih.gov/31375105/) · [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6679440/) | Uses a panel of normal samples to model **position-, strand- and nucleotide-specific** background noise in amplicon deep sequencing (Ion AmpliSeq); Poisson framework; calls SNVs down to ~1% VAF. | Almost point-for-point our method. Primary reference. |
| **iDES / CAPP-Seq "background polishing"** — Newman et al., *Nat Biotechnol* 2016 — [Nature](https://www.nature.com/articles/nbt.3520) | Builds a **position-specific background error database** from healthy controls and removes variants below position-specific thresholds ("in-silico polishing"). | Canonical paper for the position-specific-threshold concept. |
| **TNER** — *BMC Bioinformatics* 2018 — [Springer](https://link.springer.com/article/10.1186/s12859-018-2428-3) | Background-error suppression for ctDNA using **tri-nucleotide context + position**. | Refinement: context-aware error rates. |
| **Beta-binomial site-specific error** (e.g. deepSNV / shearwater family) — [PoN validation, *J Mol Diagn*](https://www.sciencedirect.com/science/article/pii/S1525157821003792) · [PMC5001245](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5001245/) | Models per-site error with a beta-binomial to capture **between-sample overdispersion**, rather than a single pooled rate. | The statistical core of our model. |
| **Mutect2 Panel of Normals + read-orientation model** (GATK) | Builds a site-specific artefact blacklist from normals and models orientation-bias artefacts (OxoG/FFPE). | Same family of ideas; production-grade precedent. |

**Our contribution / where original work sits:** (a) applying and **validating** the method on
*our* assay and instrument; (b) framing the output as a **positional Limit of Blank per CLSI
EP17** rather than a flat VAF cut-off; (c) using the laboratory's **own clinical VCFs** to strip
real calls from the blank.

---

## 3. Reflection — HSST Standards of Proficiency mapping

Where this project generates portfolio evidence against the **Consultant Clinical Scientist
Standards of Proficiency**. Dissemination-related standards are intentionally excluded because
this work is **not intended for external publication**.

### Core evidence

| Standard of Proficiency (abridged) | Evidence produced by this project |
|---|---|
| Ensure the clinical scientific validation of analytical results so that complex investigations are accurately and critically evaluated | Positional **LoB/LoD validation report (CLSI EP17)** distinguishing true signal from artefact at 1–2% VAF. |
| Provide a high level of scientific expertise to complex problems in own area of specialist practice | The statistical model (beta-binomial, position/strand/substitution) applied to a genuinely complex analytical problem. |
| Continually improve the quality of clinical scientific services by directing the introduction, evaluation and application of improved scientific and operational procedures | Replacing a flat VAF threshold with a **per-position LoB filter** in the calling pipeline; new SOP. |
| Develop and apply a strategy to optimise the impact of clinical audit to deliver outcome-focused quality improvement | Characterising per-position error rates is effectively an **audit of the assay**; measured reduction in false positives. |
| Introduce and critically evaluate measures to identify, actively manage and mitigate risk to patients | Reducing false-positive / false-negative low-VAF calls is direct **diagnostic-risk mitigation**. |
| Initiate and direct research and innovation programmes to completion, evaluate outcomes and amend service provision as appropriate | The full project lifecycle: problem → method → validation → adoption into service. |
| Evaluate published research and innovation for patient benefit and make recommendations for improvements | Review of AmpliSolve / iDES / TNER → recommendation to adopt positional LoB. |
| Direct the operation of a service to ensure compliance with local, national and internationally accepted standards and guidelines | Explicit compliance with **CLSI EP17** and **ISO 15189**. |
| Lead and shape the application of advances in science, technology, research and innovation, especially in genomics and personalised / precision medicine | Ultrasensitive low-VAF detection underpins precision oncology reporting. |

### Supporting evidence

| Standard of Proficiency (abridged) | How this project contributes |
|---|---|
| Evaluate the scientific literature and work with others to develop scientific and business cases for service improvement | Business case for reducing false positives, grounded in the literature review. |
| Assess the demand and specification for evolving scientific services with users and clinical colleagues | Problem definition with the clinical users of the test. |
| Ensure the service meets service accreditation standards | Validation documentation feeds UKAS / ISO 15189 accreditation. |
| Ensure compliance with the NHS ethical and research governance framework | HRA / ethics approval for retrospective use of BAM/VCF data. |
| Communicate complex clinical scientific and technical information in a wide range of settings and formats | Presenting the method and its impact to the laboratory / MDT. |
| Through the initiation and translation of cutting-edge scientific research, bring strategic direction, innovation and continuous improvement into practice | Translating a published method into local clinical practice. |
| Participate in clinical scientific and technical teaching, training and assessment of peers | Training clinical scientists / BMS on the new QC filter. |

### Not applicable

The clinical-facing and management standards are **not** covered by this project and should be
evidenced elsewhere: direct management of complex patients in the MDT; imparting results or
prognosis to patients/families; clean and safe physical environments (H&S); financial probity;
commissioning; staff appraisal and line management; national screening programmes; and — by
choice — formal dissemination in peer-reviewed journals and at conferences.

---

## 4. Next steps

- Scaffold the Python (pysam) package: `pileup → vcf_mask → error_model → detect → positional LoB`.
- Synthetic ground-truth data generator (blank + injected strand-biased artefacts + VCF-masked
  real mutations) for local development on a personal machine.
- Validation notebook / EP17 report on the synthetic data, then port to the laboratory.
