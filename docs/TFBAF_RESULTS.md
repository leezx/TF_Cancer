# TF–BAF module — results (v0.1)

- Design: `docs/TFBAF_DESIGN.md` (approved, review round 6). Thresholds and
  the figure selection rule are unchanged.
- Code: `scripts/tfbaf/00_freeze_baf_members.R`, `01_extract_tcga.R`,
  `02_analyze.R`.
- Run on the Mac mini, 2026-10-13. Extraction took about 3.5 min, analysis
  under 1 min.
- Tables and figures: `results/tfbaf/`. Per-sample extract and run records:
  `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/2026-10-13_tf_cancer_tfbaf_v0.1/`.

Everything here is association: correlation across tumors (TCGA) or across
cell lines (DepMap). It does not show binding, recruitment or a functional
TF–BAF partnership.

## 1. Checks (§4)

| Check | Result |
|---|---|
| Gene mapping | 33/33 genes in the TCGA TPM rows and in the CRISPR matrix |
| SKCM samples | 102 primary, 366 metastatic; 0 duplicate patients |
| TF expression (tumor median TPM) | SOX9 72, SPI1 145, FOXA1 216, GATA3 250, SOX10 310 |
| TF dependency (median Chronos effect; fraction < −0.5) | SOX10 melanoma −1.78 (86%), SPI1 AML −0.68 (56%), SOX9 CRC −0.54 (54%); GATA3 breast −0.20 (32%) and FOXA1 prostate −0.21 (40%), both flagged `TF_not_dependency` |
| BAF co-dependency control (all 1,208 cell lines) | ρ(SMARCB1, SMARCE1) = 0.56; ρ(ARID1A, SMARCB1) = 0.48; both as expected |

## 2. TCGA tumor co-expression (main analysis)

Number of the 28 BAF members called per pair:

| Pair | n tumors | strong (ρ ≥ 0.3) | moderate | negative | purity-adjusted median ρ |
|---|---|---|---|---|---|
| CRC SOX9 | 379 | 4 | 11 | 0 | 0.22 |
| AML SPI1 | 173 | 5 | 2 | 6 | 0.05 (n = 111) |
| Prostate FOXA1 | 494 | 21 | 1 | 2 | 0.53 |
| Breast GATA3 | 1,092 | 11 | 4 | 1 | 0.22 |
| Melanoma SOX10 | 102 | 12 | 4 | 0 | 0.17 |

- **Positive correlation is broad, not subunit-specific.** In most pairs,
  most BAF members correlate positively with the TF. The extreme case is
  FOXA1 in prostate, with a median ρ of 0.53 even after purity adjustment.
- These correlations have no reference set: we do not know how the TF
  correlates with random expressed genes. So they cannot be called
  BAF-specific (open point 1).
- Examples:
  - CRC SOX9: SMARCC1 0.34, SMARCA4 0.33, BICRA 0.30; SMARCB1 0.13 (none).
  - AML SPI1: BCL7B 0.60, SMARCD3 0.46, BCL7C 0.45; negative SMARCC1
    −0.30. The purity-adjusted ρ for AML is near 0 overall; BCL7A is
    flagged `attenuated_gt50pct`.
  - Breast: the luminal-only analysis (n = 782) gives a median ρ of 0.27,
    close to the main analysis.
  - Melanoma: the metastatic set (n = 366) gives a median ρ of 0.17, lower
    than primary.

## 3. DepMap co-dependency (cell lines)

Significant calls (q < 0.05; BH within each pair's 28 members, as
designed). These are **not** corrected across all 140 tests: with a single
BH over 140 tests, both hits have q = 0.087 (review N3).

- **AML SPI1–SMARCD1**: ρ = 0.56, q = 0.022, n = 32 (strong positive).
- **CRC SOX9–DPF3**: ρ = −0.40, q = 0.035, n = 63 (strong negative).
  DPF3 is flagged low-expression in TCGA COADREAD, so this is not read as a
  mechanistic cBAF partner (review N1).
- Near misses (q 0.07–0.17): AML SPI1–BRD9 0.49; breast GATA3–SMARCD1 0.37
  and –SMARCE1 0.33; melanoma SOX10–SMARCA2 0.32. Prostate FOXA1–SMARCB1
  is 0.79 with n = 10, which is underpowered and not interpreted.
- CRC SOX9–SMARCB1: ρ = −0.01. This agrees with the earlier quick check
  (Pearson −0.035) and does not support SOX9–SMARCB1 co-dependency.

## 4. Both evidence types

**No BAF member has a statistically supported concordant association
(both TCGA and DepMap) in any pair.** This is not evidence of no functional
relationship, especially for prostate (n = 10) (review N2). The
one DepMap hit, AML SMARCD1, has TCGA ρ = 0.00. The figures show this as it
is (`figures/TFBAF2_summary_dotplot.pdf`; `figures/TFBAF1_scatter.pdf`, one
page per pair, with selected members shown even when nonsignificant).

## 5. Interpretation and limits

- Observed:
  - TF and BAF mRNA levels are broadly co-correlated across tumors.
  - In cell lines, TF–BAF co-dependency is mostly absent. The exception is
    AML SPI1–SMARCD1.
- Not established: a TF-specific BAF subunit or subcomplex partner in any
  of the five cancers.
- CRC SOX9–SMARCC1 (TCGA ρ 0.34) has no DepMap co-dependency support. This
  does not bear on the SOX9–BAF chromatin experiments (dTAG,
  accessibility, SMARCC1 occupancy): those measure a different relation,
  and co-dependency is not recruitment (review N4).
- Limits:
  - There is no background reference for expression correlation (point 1).
  - DepMap uses 2D cell lines only. Organoid screens are deferred (design
    Q1).
  - Prostate DepMap n = 10.
  - GATA3 and FOXA1 are weak dependencies in these lines, so a correlation
    with them means little functionally.

### Open points (not run; for the user/reviewer to decide)

1. **Background reference** for TCGA: the percentile of each BAF member's ρ
   among the TF's correlations with all expressed protein-coding genes in
   the same tumors. This would show whether BAF members are more correlated
   with the TF than typical genes are. It is additive and changes no call.
2. Organoid (3D) CRISPR screens, after a separate download decision.

## 6. Deviations

None from the approved design.

## 7. Review

- Round 7 (PR #5, head 5d83c5e): **APPROVE**, 0 blocking.
  - N1–N4 are folded into the wording above.
  - Open point 1 (background reference): add it as v0.2 without blocking
    the merge. Use the TF's Spearman distribution against all
    expression-filtered protein-coding genes in the same tumors and report
    the gene universe. A percentile is not set-level enrichment; that needs
    a set test that accounts for expression level.
  - Organoids are a later module.
