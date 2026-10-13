# TF–BAF module: TF–BAF co-expression (TCGA) and co-dependency (DepMap) in five TF–cancer pairs — design

- Status: **approved v0.1** (review round 6: APPROVE, 0 blocking; notes
  folded in, §8). The BAF member list is frozen from annotation.
- Requested by the user on 2026-10-13. Planned with ChatGPT in the project
  review conversation; round 4 review N7 says it must stay a **separate
  module** from Step 2.
- Relation to the roadmap: this is a supporting module for the SOX9–BAF
  manuscript. It is **outside** the charter's 100% definition (Steps 0–7)
  unless the user adds it, and it does not change Step 1 or Step 2.
- Compute: Mac mini; this is a small analysis.

## 1. Question and scope

In five lineage TF–cancer pairs, which BAF (mSWI/SNF) members does the TF
track:

- (a) in mRNA expression across TCGA tumors, and
- (b) in CRISPR gene effect across DepMap models of the same lineage?

| Pair | TF | TCGA cohort | DepMap models (CRISPR, 26Q1) |
|---|---|---|---|
| CRC | SOX9 | COADREAD (COAD + READ) primary tumors, n = 379 | OncotreeLineage = Bowel, 63 cell lines |
| AML | SPI1 (PU.1) | LAML, n = 173 | OncotreePrimaryDisease = Acute Myeloid Leukemia, 32 |
| Prostate | FOXA1 | PRAD, n = 494 | OncotreeLineage = Prostate, **10** |
| Breast | GATA3 | BRCA, n = 1,092 | OncotreeLineage = Breast, 53 |
| Melanoma | SOX10 | SKCM **primary** tumors, n = 102 (new extraction, §2) | OncotreePrimaryDisease contains "Melanoma", 76 |

Out of scope:

- GTEx and TCGA normals (user: tumors only).
- Physical interaction, recruitment or causality. Expression correlation
  and co-dependency are two different kinds of association evidence;
  neither shows binding or recruitment.
- Pan-cancer discovery and other TFs.
- Organoid / 3D screens (§7 Q1): the local DepMap 26Q1 CRISPR matrix and
  ScreenSequenceMap contain only 2D screens (ScreenType 2DS/2DD). The 3D
  data (Neiswender et al. 2026) are screen-level files that are not local.

## 2. Inputs

| Input | Path | Notes |
|---|---|---|
| TCGA TPM | Step 1 `processed/v0.1/toil_<cohort>.rds` (`tpm_tf`) for TFs; **BAF genes are not TFs**, so they are read from `UCSC_Toil_TCGA_GTEx/raw/TcgaTargetGtex_rsem_gene_tpm.gz` | log2(TPM + 0.001) → linear → log2(TPM + 1) |
| SKCM samples | `TcgaTargetGTEX_phenotype.txt.gz`: 102 Primary Tumor, 366 Metastatic | one sample per patient |
| Tumor samples for the other 4 cohorts | Step 1 `results/01_tcga_tf_de/sample_manifest.tsv`, `group == "tumor"` | unchanged |
| BAF members (frozen) | `results/tfbaf/BAF_members.tsv` (from `scripts/tfbaf/00_freeze_baf_members.R`) | §3.1 |
| DepMap CRISPR | `DATA/1.Databases/DepMap3D_202608/CRISPRGeneEffect_DepMap26Q1.csv` + `Model_DepMap26Q1.csv` | Chronos gene effect; 2D cell lines |
| Purity (sensitivity) | `DATA/1.Databases/TCGA_PanCanAtlas_ABSOLUTE_purity/` | downloaded 2026-10-13, md5 = GDC |
| BRCA PAM50 (sensitivity) | `DATA/1.Databases/TCGA_PanCanAtlas_subtypes/` | `Subtype_mRNA` |

## 3. Method / contract

### 3.1 BAF members (frozen)

- **28 genes**, mapped by Ensembl ID through HGNC to GENCODE v23. BICRA and
  BICRAL are GLTSCR1 and GLTSCR1L in v23.
  - ATPase: SMARCA4, SMARCA2.
  - Core: SMARCC1, SMARCC2, SMARCD1–3, ACTL6A, ACTL6B, BCL7A–C, SMARCB1,
    SMARCE1, SS18, SS18L1.
  - cBAF-specific: ARID1A, ARID1B, DPF1–3.
  - PBAF-specific: PBRM1, ARID2, BRD7, PHF10.
  - ncBAF-specific: BRD9, BICRA, BICRAL.
- Each gene has cBAF / PBAF / ncBAF membership flags (after Mashtalir et al.
  2018; Centore et al. 2020).
- ACTB is excluded (ubiquitous, common-essential). BCL11A/B are excluded
  (Lambert TFs, BAF-associated only).
- ARID2, DPF1 and DPF3 are also in the Lambert TF list. This is recorded
  as `in_lambert_TF` and does not exclude them.
- Membership flags are a structural annotation. Some core subunits are
  context-dependent, so the flags do not show that a TF depends on one
  specific subcomplex (review N3).

### 3.2 TCGA expression correlation (primary)

- Per pair: primary tumors, one per patient, as in Step 1 (SKCM: 102
  Primary Tumor; the patient is the first 12 barcode characters).
- E = log2(TPM + 1).
- Spearman ρ(TF, BAF member) for each of the 28 members. Asymptotic
  p-value; BH across the 28 members within each pair.
- Report: n, ρ, p, q, TF median TPM, BAF member median TPM.
- **Calls** (preset, exploratory thresholds, not biological boundaries):
  - `strong` ρ ≥ 0.3 and q < 0.05;
  - `moderate` 0.2 ≤ ρ < 0.3 and q < 0.05;
  - negative calls the same way with ρ ≤ −0.2 / −0.3;
  - otherwise `none`.
- A BAF member with tumor median TPM < 1 is reported but flagged
  `low_expression`.

### 3.3 TCGA sensitivity analyses (reported, not used for calls)

1. Purity-adjusted partial Spearman (ABSOLUTE called purity; rank
   residualization as in Step 2 §3.5), on tumors with purity. n is
   reported; this includes LAML (64% coverage).
2. BRCA: luminal only (PAM50 LumA + LumB), because GATA3 is a luminal TF and
   luminal-vs-basal heterogeneity can drive correlations.
3. SKCM: metastatic tumors (n = 366) analysed separately.

### 3.4 DepMap co-dependency (primary)

- Models: in the CRISPR matrix, `ModelType = Cell Line`, lineage as in §1.
  One row per ModelID (the matrix is already model-level).
- Spearman ρ(TF gene effect, BAF gene effect) over models with both values.
  BH across the 28 members within each pair. The same call thresholds as
  §3.2 apply; positive ρ = co-dependency.
- Minimum n = 10. **Prostate (n = 10) is at the floor**; it is reported but
  labelled `underpowered` and not interpreted.
- Context, reported per pair: the TF's gene-effect distribution (median,
  fraction < −0.5) and each BAF member's median effect. A correlation
  between two non-dependent genes means little; such pairs are flagged
  `TF_not_dependency` if the TF's median effect is > −0.3. TF dependency
  prevalence (fraction of models < −0.5) is always reported alongside. A
  correlation without TF dependency is not read as functional coupling
  (review N2).
- Prostate's ρ values are reported but never used to rank or interpret
  (review N1).

### 3.5 Figures

1. **Scatter panels** (Fig TFBAF-1). The selection rule is fixed now, to
   avoid picking after looking:
   - Per pair, 3 BAF members: SMARCA4 (fixed reference), the member with
     the highest TCGA ρ, and the member with the highest DepMap ρ (or the
     next-ranked if a member repeats).
   - Each member is shown in a TCGA panel and a DepMap panel. A selected
     member that misses the call thresholds is still shown, labelled
     `nonsignificant` (review N4).
   - One point is one tumor or one model. Each panel shows ρ, q and n; the
     trend line is for display only.
2. **Summary dual-evidence dot plot** (Fig TFBAF-2):
   - y axis: the 28 BAF members, grouped ATPase / core / cBAF / PBAF /
     ncBAF; x axis: the 5 pairs, with two dots per cell (TCGA, DepMap).
   - Colour: ρ (diverging, 0-centred). Size: −log10 q. Black outline:
     call ≥ moderate. Grey: not evaluable.
   - Members supported by both data types are listed in the caption. With
     no shared support, the figure says so; it does not search for
     agreement.
   - Cross-pair comparison uses continuous ρ, n and q, not the count of
     positive calls, because power differs by n (review Q3).

## 4. Validity checks

1. **Gene mapping**: all 28 + 5 genes present in the TCGA TPM rows (by v23
   ID) and in the CRISPR matrix (by symbol + Entrez in the column name).
   Already verified for CRISPR: 33/33.
2. **TF expression**: each TF is expressed in its TCGA cohort (Step 1
   identity controls already cover SOX9, SPI1, FOXA1, GATA3). SOX10 in SKCM:
   tumor median TPM ≥ 1 expected.
3. **TF dependency** (expected biology, not a gate): SOX10 in melanoma and
   SPI1 in AML lines are known dependencies (median effect < −0.5
   expected). SOX9 in CRC (median −0.54 from the local quick check), FOXA1
   and GATA3 are reported. A miss is reported and discussed, not tuned.
4. **Co-dependency positive control**: within BAF, ρ(SMARCB1, SMARCE1) and
   ρ(ARID1A, SMARCB1) over all CRISPR cell lines are expected to be
   positive (cBAF subunits are co-essential). This checks the data
   handling; it is not a result.
5. **SKCM sample check**: 102 primary tumors, no duplicate patients.

## 5. Expected outputs

Repo (`results/tfbaf/`, small TSV + PDF):

- `BAF_members.tsv` (frozen in this PR)
- `TCGA_TF_BAF_correlation.tsv`: pair, TF, BAF, class, n, ρ, p, q, call,
  medians, `low_expression`
- `TCGA_TF_BAF_sensitivity.tsv`: purity-adjusted, BRCA luminal, SKCM
  metastatic
- `DepMap_TF_BAF_codependency.tsv`: pair, n, ρ, p, q, call, TF/BAF median
  effect, flags
- `TF_BAF_summary.tsv`: one row per pair × BAF member with both evidence
  types
- `qc/` (mapping, TF dependency, controls)
- `figures/TFBAF1_scatter.pdf`, `figures/TFBAF2_summary_dotplot.pdf`

DATA: `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/2026-10-13_tf_cancer_tfbaf_v0.1/`
(the SKCM TPM extract, per-sample tables, run records).

Results doc: `docs/TFBAF_RESULTS.md`.

## 6. Order of work after approval

1. Extract TPM for the TFs and BAF genes for the 5 cohorts (SKCM from raw),
   run the checks in §4.
2. TCGA correlations and sensitivities, then DepMap co-dependency, then
   figures.

## 7. Open questions for review

- **Q1. Organoids.** The user asked for cell lines **and** organoids. 3D
  organoid CRISPR screens (DepMap 3D, Neiswender et al. 2026) are not in
  the local 26Q1 model-level matrix. Should we (a) run cell lines now and
  add organoids later from the screen-level 3D files, after a separate
  download decision, or (b) block until they are available?
- **Q2. SKCM sample type.** Primary (n = 102) as the main analysis and
  metastatic (n = 366) as a sensitivity, or the reverse? TCGA SKCM is mostly
  metastatic, and metastatic samples add site composition (lymph node,
  skin).
- **Q3. Thresholds.** Are ρ ≥ 0.3 (strong) and ρ ≥ 0.2 (moderate) with
  q < 0.05 reasonable for both TCGA (n = 102–1,092) and DepMap (n = 10–76)?
  With n = 1,092, ρ = 0.2 is very significant; with n = 32, it is not.
- **Q4. Raw vs adjusted.** Raw Spearman as the primary TCGA statistic (the
  user asked for a simple analysis) with purity-adjusted as sensitivity.
  Is that acceptable, given that BAF subunits are broadly expressed and
  purity can drive TF–BAF correlations?

## 8. Review history

- **Round 6** (head 5ba809e): **APPROVE**, 0 blocking.
  - N1: prostate underpowered → §3.4.
  - N2: report TF dependency prevalence → §3.4.
  - N3: membership flags are annotation only → §3.1.
  - N4: show nonsignificant selected members → §3.5.
  - Q1: run cell lines now; organoids later, as a separate validation never
    merged with 2D gene effects.
  - Q2: SKCM primary is the main analysis, metastatic a sensitivity.
  - Q3: thresholds accepted; compare continuous ρ, n and q across pairs.
  - Q4: raw Spearman is primary; if the purity-adjusted value reverses or
    strongly attenuates, flag it in the results.
