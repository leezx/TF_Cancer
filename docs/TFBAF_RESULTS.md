# TF–BAF module — results (v0.1, v0.2)

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

### Open points from v0.1 (both done in v0.2, §8)

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

## 8. v0.2: organoids, DepMap context, background percentile

- Design: `docs/TFBAF_V02_DESIGN.md`, approved in round 8 (PR #6).
- Code: `scripts/tfbaf/10_v02_extract.R` and `11_v02_analyze.R`.
- Run on the Mac mini on 2026-10-13. Extraction took about 2 min,
  analysis about 2 min.
- Logs: `logs/tfbaf_v02_*.log`.
- Records: `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/2026-10-13_tf_cancer_tfbaf_v0.2/`
  (`organoid/`, `background/`, `report/`).
- v0.1 tables and figures are unchanged. The v0.2 figures have their own
  files (`*_v02.pdf`).

### 8.1 Checks

| Check | Result |
|---|---|
| Organoid inventory | 178 screens / 162 organoids. Colorectal: 97 screens, 85 organoids, 83 individuals, 95 malignant screens |
| CRC units (§3.1) | 82 individuals: 68 `minLib`, 14 `v1_1`; 57 primary-derived, 25 metastasis-derived |
| Library LFC (N1) | SOX9 median LFC: `minLib` −0.44, `v1_1` −0.48 (`qc/organoid_library_LFC.tsv`) |
| Organoid positive control (all 178 screens) | ρ(SMARCB1, SMARCE1) = 0.43; ρ(ARID1A, SMARCB1) = 0.36; both positive as expected |
| Depth reference set | 295 genes depleted in ≥ 90% of screens |
| Background consistency | The background ρ for the 28 members reproduces v0.1 TCGA and DepMap ρ, and the organoid main ρ, exactly (asserted in code) |
| Background universe (N5) | TCGA 12,233–13,354 protein-coding genes (median TPM ≥ 1); DepMap 17,204–17,930 genes (600–1,326 with some missing values, scored on complete cases); organoid 16,561; 0 zero-variance genes; TF excluded |
| Organoid availability | CRC analysed. AML, prostate, breast and melanoma unavailable: there are no matching organoids |

### 8.2 Organoid CRISPR knockout-effect correlation (CRC SOX9, n = 82)

- **SOX9 is a dependency in this panel.** It is depleted in 34% of the
  units (`is_depleted`), with a median LFC of −0.45.
- **One member is called: ARID1A**, ρ = 0.35, 95% CI [0.12, 0.53],
  q = 0.031 (strong).
  - Depth-adjusted ρ = 0.36, so it is not depth-sensitive.
  - `minLib` only: 0.40 (q 0.024).
  - Primary-derived only (n = 57): 0.29 (q 0.25, no call). The direction
    is the same, but the estimate is weaker and the smaller subset loses
    significance. This alone neither shows the result is fragile nor can
    be put down to n (review N2).
- q is BH within the 28 members of this pair. It is not a global FDR
  across cancers, data types and the exploratory v0.2 analyses
  (review N5).
- SMARCB1 has ρ = 0.31 (CI [0.09, 0.50]), q = 0.059: no call in the main
  analysis. It is called in the primary-derived subset (0.42, q 0.035).
  Main SMARCB1 is not called, so it is reported, not interpreted.
- **ARID1A is rarely depleted (7%).** Its pattern is therefore
  "TF-specific dependency", not "both dependent, correlated". This is a
  **correlation of CRISPR knockout effects, not co-dependency**
  (review N1). The ρ reflects SOX9 and ARID1A knockout effects varying
  together across organoids, mostly in a range where ARID1A loss has little growth effect.
  It does not show that organoids depend on ARID1A together with SOX9.
  ARID1A is also a known tumour suppressor in CRC.
- Other v0.1-highlighted members: SMARCA4 0.17 and SMARCC1 0.19, both
  nonsignificant.

### 8.3 DepMap 2D context

- **CIs** for the two v0.1 hits:
  - AML SPI1–SMARCD1: 0.56 [0.28, 0.76].
  - CRC SOX9–DPF3: −0.40 [−0.59, −0.17].
  - The other CIs are wide at these n. For example, CRC SMARCB1 is −0.01
    [−0.27, 0.26].
- **Dependent BAF members** (≥ 25% of lines < −0.5):
  - ACTL6A, SMARCB1, SMARCE1 and DPF2 are dependent in all five lineages.
  - AML has 11 dependent members, CRC 4.
- **Patterns:**
  - Only AML SPI1–SMARCD1 is "both dependent, correlated".
  - In all pairs, the common-essential members (ACTL6A, SMARCB1, SMARCE1,
    DPF2) are "both dependent, no positive co-dependency". This is the
    case the brief warned about: a member can be essential in every line
    and still not track the TF.
  - The rest are "TF-specific dependency". All five TFs are dependent at
    ≥ 25% prevalence (GATA3 32%, FOXA1 40%).
  - GATA3 and FOXA1 also carry the v0.1 flag `TF_not_dependency`
    (median > −0.3).

### 8.4 Background percentile (descriptive; not a set-level test)

Median signed percentile of the 28 members' ρ among all genes' ρ with the
same TF:

| Pair | TCGA (bg median ρ) | DepMap 2D | Organoid |
|---|---|---|---|
| CRC SOX9 | **85** (0.10); 8 of 28 ≥ 95th | 53 | 64; 4 of 28 ≥ 95th |
| AML SPI1 | 44 (0.06) | 79 | — |
| Prostate FOXA1 | 71 (0.30) | 53 (exploratory, n = 10) | — |
| Breast GATA3 | 58 (0.16) | 61 | — |
| Melanoma SOX10 | 67 (0.16) | 40 | — |

- **CRC is the pair where the BAF correlations stand out from the
  background in TCGA.** A percentile is not a BAF gene-set enrichment
  test. It does not show that SOX9 preferentially associates with or
  regulates BAF (review N3). SMARCC1 is at the 99.2nd percentile, SMARCA4 the
  99.1st, and ARID1A the 95.1st. SMARCB1 is at the 63rd.
- In prostate, FOXA1 correlates with most genes (background median 0.30).
  Its high BAF ρ (median 0.51) is therefore only modestly above typical
  genes. This is the v0.1 "broad, not subunit-specific" finding, now
  quantified.
- The organoid ARID1A ρ is at the 99.9th percentile. So are DepMap AML
  SMARCD1 (99.9th) and the |ρ| of CRC DPF3 (99.7th).

### 8.5 Across evidence types, CRC SOX9

| Member | TCGA ρ (q; pct) | DepMap 2D ρ [CI] (q) | Organoid ρ [CI] (q) |
|---|---|---|---|
| ARID1A | 0.27 (6e-7, moderate; 95th) | 0.21 [-0.07, 0.48] (0.52) | **0.35 [0.12, 0.53] (0.031)** |
| SMARCC1 | 0.34 (5e-10, strong; 99th) | −0.01 [−0.28, 0.26] | 0.19 (0.40) |
| SMARCA4 | 0.33 (5e-10, strong; 99th) | −0.20 [−0.43, 0.05] | 0.17 (0.40) |
| SMARCB1 | 0.13 (0.014, below ρ 0.2: none; 63rd) | −0.01 [−0.27, 0.26] | 0.31 [0.09, 0.50] (0.059) |

- **ARID1A is the first member with called support from two independent
  sources:** TCGA moderate co-expression, and correlated knockout effects
  in organoids.
  DepMap 2D points the same way but is not significant.
- It is one of 28 tests, with q = 0.031 within the pair. It is weakened in
  the primary-derived subset, and ARID1A itself is rarely a dependency.
- It is a cross-platform association candidate to follow up, not a
  validated SOX9–cBAF partner. It shows neither joint dependency nor a
  direct mechanism. Co-expression and correlated knockout effects are not
  recruitment.
- SOX9–SMARCC1 and SOX9–SMARCA4 stand out in CRC tumour co-expression
  (99th percentile). They have no matching 2D or 3D knockout-effect
  correlation.

### 8.6 Deviations from the v0.2 design

- **TCGA protein-coding proxy.** The design gives the definition
  (HGNC Approved + UniProt ID, by Ensembl ID). It is a reproducible proxy,
  not the strict GENCODE protein-coding biotype (review N4).
  - This finds 19,995 of 60,498 Toil genes.
  - BICRA and BICRAL carry their GENCODE v23 names (GLTSCR1, GLTSCR1L) in
    the matrix and were mapped by `v23_symbol`.
- No other deviations.

## 9. Review (v0.2)

- Round 9 (PR #7, head c4afaee): **APPROVE**, 0 blocking.
  - N1: ARID1A is a "correlated CRISPR knockout effect", not
    co-dependency → §8.2, §8.5.
  - N2: primary-derived wording → §8.2.
  - N3: a percentile is not set enrichment → §8.4.
  - N4: HGNC + UniProt is a proxy → §8.6.
  - N5: the q scope is within-pair → §8.2.
  - The reviewer advises against further threshold or rule changes to get
    a stronger positive result. Their review was based on the results
    document only: they did not recompute or audit the code.
