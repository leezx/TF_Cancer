# Step 2: within-cancer TF co-expression and proliferation association — design

- Status: **draft for review** (design only; no compute, no new data downloaded yet)
- Answers: Q1, second half ("co-upregulated" read as within-patient
  co-expression, Task-1 §四 B); feeds Q2 (Step 5) and Q3 (Step 7).
- Builds on: Step 1 v0.1 (PR #2, frozen). Step 1 calls and thresholds are
  **not** changed here.
- Review round 3 placement: the proliferation axis belongs in this step, as
  an annotation and ranking axis, not an exclusion rule.
- Compute: Mac mini (`computingResource.md`); per-cohort matrices are at
  most about 1,000 TFs × 1,092 samples.

## 1. Scope

In scope:

1. **Proliferation association** per cohort: how strongly each TF tracks a
   predefined proliferation score across tumors, with and without purity
   adjustment.
2. **Within-cancer TF–TF co-expression** across tumor samples, adjusted for
   purity and proliferation (primary) and unadjusted (secondary).
3. **Cross-cancer recurrence** of TF–TF edges by organ system, and a
   pan-cancer recurrent-edge network with modules (Fig 3).
4. **Tumor vs normal co-expression** (descriptive, supporting) in cohorts
   with ≥ 40 TCGA normals.
5. A **prospectively defined focused Fig 1B view**: Step 1 recurrent TFs that
   are not recurrently proliferation-coupled.

Out of scope:

- Lineage-specificity and lineage-identity axes → next design, together
  with Step 3 (annotation).
- Malignant-cell validation (Step 5): bulk co-expression can still reflect
  composition after purity adjustment.
- BRCA PAM50-stratified analysis (open question Q3 below).
- Any causal, physical or regulatory interpretation. An edge means
  correlated mRNA abundance across patients, not binding, interaction or
  co-occupancy (Task-1 Fig 3 caveat).
- Pooled cross-cancer correlation. Every correlation is computed within one
  cohort; pooling cohorts creates correlations from between-cancer mean
  differences (Task-1 §四).

## 2. Inputs

| Input | Path | Notes |
|---|---|---|
| Per-cohort counts + samples | `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/processed/v0.1/toil_<cohort>.rds` | Step 1, all genes; tumor and normal samples as in the Step 1 manifest |
| Step 1 tables | `results/01_tcga_tf_de/` (`TF_reference.tsv`, `TCGA_TF_DE_all.tsv`, `TF_recurrent_upregulation.tsv`, `cohort_manifest.tsv`) | frozen v0.1 |
| Gene ID → symbol | `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/raw/gencode.v23.annotation.gene.probemap` | md5 already recorded |
| Proliferation genes | `DATA/genesets/h.all.v2026.1.Hs.symbols.gmt` (MSigDB Hallmark v2026.1) | local |
| **Tumor purity (new)** | `DATA/1.Databases/TCGA_PanCanAtlas_ABSOLUTE_purity/raw/TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | TCGA PanCanAtlas ABSOLUTE calls, GDC file 4f277128-f793-4354-a13d-30cc7fe9f6b5, 901,812 B, GDC Content-MD5 8ea2ca92c8ae58350538999dfa1174da. **Download needs user confirmation**; verify md5 and record `link.md` + `checksums.md5` |

Fallback if ABSOLUTE purity covers < 70% of a cohort's tumors: the
cohort's purity-adjusted analysis runs on the covered subset and is
flagged. No expression-derived purity estimate is introduced in this step
without a design amendment.

## 3. Method / contract

### 3.1 Samples

- Tumors: the Step 1 sample manifest `group == "tumor"` (one primary tumor
  per patient), for all 16 cohorts including LAML, GBM, LGG and PAAD.
  Co-expression needs no normal reference.
- Normals (§3.6 only): TCGA solid-tissue normals for cohorts with ≥ 40:
  BRCA 113, KIRC 72, LUAD 59, THCA 59, PRAD 51, COADREAD 51, LUSC 50,
  LIHC 50, HNSC 44. GTEx is not used for co-expression.
- Purity join: ABSOLUTE `array` (15-character barcode) = Toil sample ID;
  keep `call status == "called"`. Coverage per cohort is reported.

### 3.2 Expression scale

Per cohort and per sample group separately: DESeq2
`vst(blind = TRUE)` on the Step 1 rounded counts for all genes, computed
on that cohort's tumors only (or normals only, for §3.6). Then subset rows
to TFs and proliferation genes. Gene symbols come from the GENCODE v23
probemap; ambiguous or unmapped symbols are reported and dropped.

### 3.3 TFs analysed per cohort

A TF enters a cohort's analysis if tumor median TPM ≥ 1 and tumor
prevalence (TPM ≥ 1) ≥ 25% in that cohort, the same expression rule as
Step 1. This avoids correlations driven by zero inflation. The per-cohort
TF universe size is reported.

### 3.4 Proliferation score and TF–proliferation association

- Signature (fixed now): the union of HALLMARK_G2M_CHECKPOINT and
  HALLMARK_E2F_TARGETS (327 genes) **minus the 25 Lambert TFs it contains**
  (ATF5, CENPA, CTCF, DR1, E2F1–4, E2F8, FOXN3, HIF1A, HMGA1, HOXC10,
  MEIS1, MEIS2, MTF2, MYBL2, MYC, PURA, SMAD3, TFDP1, DNMT1, MXD3, PA2G4,
  TP53) **and minus MKI67** (held out as a check) = 301 genes. Removing TFs
  avoids scoring a TF against itself.
- Score: per cohort, the mean of per-gene z-scores (vst, across that
  cohort's tumors) over signature genes present.
- Association per TF and cohort:
  - unadjusted Spearman ρ(TF, score);
  - **purity-adjusted partial Spearman ρ (primary)**. Compute it by
    rank-transforming TF, score and purity, residualizing TF and score on
    purity rank (linear), and taking the Pearson correlation of the
    residuals. Test with df = n − 3; BH within cohort.
- **Proliferation-coupled** in a cohort: partial ρ ≥ 0.5 and q < 0.05.
- **Recurrently proliferation-coupled**: coupled in ≥ 50% of the cohorts
  where the TF is analysed, and in ≥ 3 cohorts.

These are preset annotation thresholds, not biological boundaries.
Coupled TFs are labelled, not removed.

### 3.5 TF–TF co-expression (tumors)

- For each cohort, over all pairs of analysed TFs:
  - **primary**: partial Spearman ρ adjusted for purity and proliferation
    score (the same residualization on both covariate ranks; df = n − 4);
  - secondary: unadjusted Spearman ρ.
- BH over all pairs within the cohort.
- **Edge** (preset): |partial ρ| ≥ 0.4 and q < 0.05. Positive and negative
  edges are kept separately. Recurrence and the network use positive edges
  only (co-upregulation question).
- Cohorts without purity coverage use proliferation-only adjustment for
  that cohort, flagged `purity_adjusted = FALSE`.

### 3.6 Tumor vs normal co-expression (supporting, descriptive)

For the 9 cohorts in §3.1, compute the same primary statistic in normals,
adjusted for proliferation only (no purity for normals). For each tumor
edge, report the normal ρ and a Fisher z difference test (BH within
cohort). Label an edge `tumor_specific` if normal ρ < 0.2 and the
difference q < 0.05, and `shared` if the normal edge also passes §3.5. With
n ≈ 44–113 normals, power is limited, so these labels are descriptive and
are not used for ranking.

### 3.7 Recurrence and network

- Organ systems as in Step 1 plus blood (LAML), giving 12 organ systems. A
  positive edge is present in an organ system if it passes §3.5 in any
  cohort of that system. Cohort counts are reported too.
- **Recurrent edge**: present in ≥ 3 organ systems, the same preset as
  Step 1 Fig 1B.
- Network: nodes are TFs, edges are recurrent edges, and edge weight is the
  number of organ systems. Modules come from `igraph::cluster_leiden`
  (objective = modularity, resolution 1, n_iterations = 10, seed 1).
  Modules with ≥ 5 TFs are reported. Node annotations: DBD, Step 1
  recurrence, recurrent proliferation coupling.
- **SOX9 reference (CRC only)**: SOX9's positive edges in COADREAD and their
  recurrence in other organ systems are reported. This is descriptive, not
  a selection criterion (charter non-goal: SOX9 is not assumed to be
  pan-cancer).

### 3.8 Focused Fig 1B view (prospective, additive)

Step 1 TFs with `fig1B_primary == TRUE` (C1-up in ≥ 3 organ systems) that
are **not** recurrently proliferation-coupled (§3.4). The 358-TF Step 1
Fig 1B stays the prespecified discovery figure; this view is a separate,
labelled secondary figure. LBX2 keeps its `pipeline_sensitive` flag.

## 4. Coverage / validity checks

1. **Purity coverage**: fraction of tumors with called ABSOLUTE purity per
   cohort. Below 70% is flagged and handled as in §2. Report the ρ between
   purity and proliferation score per cohort.
2. **Proliferation score check**: Spearman ρ(score, MKI67), with MKI67
   held out of the score. ≥ 0.7 is expected in every carcinoma cohort; a
   miss is investigated, not waived.
3. **Adjustment check**: E2F1–FOXM1 (both proliferation TFs) should have a
   high unadjusted ρ that is attenuated after proliferation adjustment.
   Report both values per cohort. This only shows the adjustment acts; it
   is not a result.
4. **Positive-control pairs** (preset): FOXA1–GATA3 (BRCA), CDX1–CDX2
   (COADREAD), HNF1A–HNF4A (LIHC), TP63–SOX2 (LUSC). Expected: positive
   primary partial ρ with q < 0.05. A miss is investigated.
5. **Permutation null**: per cohort, permute each TF's adjusted residual
   vector independently across samples (seed 1, one permutation), recompute
   ρ and count edges passing the §3.5 rule. The expected count is near 0.
   Report the observed/permuted edge ratio.
6. **Gene mapping**: number of signature genes and TFs mapped per cohort;
   duplicated symbols are reported.
7. **Sensitivity**: edge overlap (Jaccard) between the primary and
   unadjusted analyses per cohort, to show how much purity and
   proliferation drive co-expression.

## 5. Expected outputs

Repo (`results/02_tf_coexpression/`, small TSV):

- `cohort_qc.tsv`: n tumors, n with purity, TF universe size, ρ(purity,
  prolif), ρ(score, MKI67), permutation edge ratio, edge counts.
- `TF_proliferation_association.tsv`: TF, cohort, ρ unadjusted, partial ρ,
  q, coupled.
- `TF_proliferation_summary.tsv`: per TF n cohorts analysed, n coupled,
  `recurrently_proliferation_coupled`.
- `TF_coexpression_edges.tsv`: cohort, TF_a, TF_b, ρ unadjusted, partial
  ρ, q, sign, `purity_adjusted`, normal ρ, tumor_specific / shared label.
  Edges only.
- `TF_pair_recurrence.tsv`: pair, n organ systems, n cohorts, organ-system
  list.
- `TF_network_modules.tsv`: TF, module, degree, annotations.
- `TF_focused_fig1B.tsv` (§3.8).
- `qc/`: controls, adjustment check, permutation null, gene mapping.
- Figures: `Fig3A_<cohort>_TF_coexpression_heatmap.pdf` (top recurrent
  TFs per cohort, primary ρ), `Fig3B_recurrent_TF_network.pdf`,
  `FigS_TF_proliferation_coupling_heatmap.pdf`,
  `Fig1B_focused_nonproliferative_TFs.pdf`.

DATA (`DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/<analysis_id>/`,
analysis_id `2026-10-13_tf_cancer_step2_coexpr_v0.1`): full per-cohort
correlation matrices (`coexpr/<cohort>_{partial,unadjusted}.rds`), vst
matrices, run records.

Results doc: `docs/STEP2_TF_COEXPRESSION_RESULTS.md`.

## 6. Order of work after approval

1. Download and verify the ABSOLUTE file (after user OK).
2. Run the smoke test on COADREAD (purity coverage, score check, controls,
   permutation null). Stop if a check fails unexplained.
3. Run all 16 cohorts, then recurrence, network and figures.

## 7. Open questions for review

- **Q1. Thresholds.** Are |partial ρ| ≥ 0.4 for edges and ≥ 0.5 for
  proliferation coupling reasonable presets, given n from 66 (KICH) to
  1,092 (BRCA)? Or should the edge rule be rank-based (e.g. top 1% of
  pairs per cohort) to equalize power?
- **Q2. Adjusting for proliferation.** Adjusting co-expression for
  proliferation removes proliferation-driven modules, which are biology,
  not artifact. Is "adjusted as primary, unadjusted as secondary" the right
  split for the co-upregulation question?
- **Q3. BRCA subtypes.** PAM50 heterogeneity will drive strong
  luminal-versus-basal TF modules in BRCA. Should Step 2 add a
  within-subtype sensitivity analysis (it needs a PAM50 annotation
  download), or defer it with the Step 1 subtype follow-up?
- **Q4. Organ-system rule.** For organ systems with several cohorts (lung,
  kidney, brain), should edge presence need any cohort (as Step 1) or all
  cohorts in that system?
- **Q5. Purity source.** ABSOLUTE (DNA-based) is chosen over
  expression-derived purity (ESTIMATE) to avoid circularity with the
  expression being analysed. Acceptable?
