# Step 2: within-cancer TF co-expression and proliferation association — design

- Status: **rev 2** after review round 4 (REQUEST CHANGES, B1–B4; see §8). Design only;
  no expression compute and no new data downloaded yet.
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
- TF–BAF expression / co-dependency in five TF–cancer pairs: a separate
  analysis module, not part of this step's gates (review N7).
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
| Proliferation signature (frozen) | `results/02_tf_coexpression/proliferation_signature_genes.tsv` | built by `scripts/02_tf_coexpression/00_freeze_proliferation_signature.R` from annotation only (§3.4) |
| **BRCA PAM50 (new, §3.9)** | TCGA PanCanAtlas sample subtype table (`Subtype_mRNA` for BRCA), exact file fixed at implementation | **download needs user confirmation**; md5 + `link.md` |
| **Tumor purity (new)** | `DATA/1.Databases/TCGA_PanCanAtlas_ABSOLUTE_purity/raw/TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | TCGA PanCanAtlas ABSOLUTE calls, GDC file 4f277128-f793-4354-a13d-30cc7fe9f6b5, 901,812 B, GDC Content-MD5 8ea2ca92c8ae58350538999dfa1174da. **Download needs user confirmation**; verify md5 and record `link.md` + `checksums.md5` |

Purity missingness is handled by the two-analysis rule in §3.5 (review B1).
No expression-derived purity estimate is introduced in this step without a
design amendment.

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
`vst(blind = TRUE, nsub = 1000, fitType = "parametric")` on the Step 1
rounded counts for all genes with no prior expression filter; size factors
are estimated by `vst` on the same samples. vst is computed on that cohort's
tumors only (or normals only, for §3.6). Rows are then subset to TFs and
signature genes **by GENCODE v23 gene ID**, not by symbol. The parameters and
the DESeq2 version are written to the run record (review N1). Tumor and
normal vst values are never compared numerically; only correlations are
compared.

### 3.3 TFs analysed per cohort

A TF enters a cohort's analysis if tumor median TPM ≥ 1 and tumor
prevalence (TPM ≥ 1) ≥ 25% in that cohort, the same expression rule as
Step 1. This avoids correlations driven by zero inflation. The per-cohort
TF universe size is reported.

### 3.4 Proliferation score and TF–proliferation association

- Signature (**frozen in this PR** as
  `results/02_tf_coexpression/proliferation_signature_genes.tsv`, review B3):
  - Start from the union of HALLMARK_G2M_CHECKPOINT and HALLMARK_E2F_TARGETS
    (MSigDB v2026.1): 327 symbols.
  - Map each symbol to Ensembl through the HGNC approved-symbol table
    (`DATA/1.Databases/HGNC_gene_id_mapping/processed/v0.1/`), then map
    Ensembl to the GENCODE v23 ID. 13 symbols have older names in v23 (e.g.
    H2AX = H2AFX, KNL1 = CASC5, NSD2 = WHSC1); they are mapped by Ensembl ID,
    not by symbol.
  - Result: 327 mapped, 0 unmapped, 0 ambiguous, PAR_Y copies dropped.
  - **Exclusions are matched by Ensembl ID**:
    - the 25 Lambert TFs in the union (ATF5, CENPA, CTCF, DR1, E2F1–4, E2F8,
      FOXN3, HIF1A, HMGA1, HOXC10, MEIS1, MEIS2, MTF2, MYBL2, MYC, PURA, SMAD3,
      TFDP1, DNMT1, MXD3, PA2G4, TP53);
    - MKI67, held out as a check.
  - Final signature: 301 gene IDs. No analysed TF can enter the score under
    another symbol or ID. The implementation asserts that the signature IDs
    and the TF IDs do not overlap.
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
Coupled TFs are labelled, not removed. Two definitions are reported
separately: purity-adjusted (primary-eligible cohorts, §3.5) and
unadjusted (all cohorts). There is no "proliferation-only" variant here,
because the score is itself the associated variable and is never
residualized on itself (round 5 N1). For each, the
outputs give n_evaluable, n_coupled and n_coupled/n_evaluable (review N3).

Interpretation (review B3): the score is derived from the same expression
data. A TF–TF correlation that disappears after adjustment is
**proliferation-associated covariance**, which may be real cell-cycle
biology, not removed technical noise. Unadjusted results are kept in full.

### 3.5 TF–TF co-expression (tumors)

Three analyses, each defined in advance and never mixed (review B1):

| Analysis | Adjustment | Samples | Cohorts |
|---|---|---|---|
| **Primary** | purity + proliferation (df = n − 4) | tumors with called ABSOLUTE purity | **primary-eligible** cohorts only: purity coverage ≥ 70% of tumors **and** ≥ 50 tumors with purity |
| Sensitivity (proliferation-only) | proliferation (df = n − 3) | all tumors | all 16 cohorts |
| Unadjusted | none | all tumors | all 16 cohorts |

- Partial Spearman ρ: rank-transform each variable, residualize on the
  covariate ranks (linear), then take the Pearson correlation of the
  residuals. Test with the df shown.
- Over all pairs of analysed TFs in the cohort; BH over all pairs within
  the cohort and analysis.
- **Edge** (preset, the same rule for every analysis): |ρ| ≥ 0.4 and
  q < 0.05. Positive and negative edges are kept separately. Recurrence and
  the network use positive edges only (co-upregulation question).
- A cohort that is not primary-eligible contributes **no primary edges**.
  Its sensitivity edges are never relabelled as primary.
- Every edge row carries `analysis` ∈ {primary, prolif_only, unadjusted}.
  All three are kept in full (review Q2).

### 3.6 Tumor vs normal co-expression (supporting, descriptive)

For the 9 cohorts in §3.1, the comparison uses the **proliferation-only**
statistic in both tumors and normals, so the estimand is the same on both
sides; normals have no purity (review B2).

- Report tumor ρ, normal ρ and Δρ = tumor − normal for every tumor edge
  (proliferation-only analysis) in these cohorts.
- No Fisher z test (review B2): rank transformation and covariate
  adjustment make its variance formula uncalibrated.
- For recurrent edges (§3.7) only, add a percentile bootstrap 95% CI for
  Δρ: 1,000 resamples, seed 1, resampling tumors and normals separately,
  recomputing the residualization in each resample.
- Descriptive labels (review N6 terminology):
  - `tumor_enriched_coexpression`: the Δρ CI excludes 0 and normal ρ < 0.2;
  - `shared_coexpression`: the normal ρ also passes the §3.5 edge rule.
- With n ≈ 44–113 normals, power is limited. These labels are descriptive,
  are not used for ranking, and do not imply malignant-cell-specific
  regulation.

### 3.7 Recurrence and network

- Organ systems as in Step 1 plus blood (LAML), giving 12 organ systems.
- **Pair-level eligibility** (review B4):
  - A TF pair is *evaluable* in a cohort for an analysis if both TFs pass
    §3.3 in that cohort and the cohort is eligible for that analysis.
  - A pair is *evaluable* in an organ system if it is evaluable in ≥ 1
    cohort of that system.
  - It is *positive* in an organ system if it is a positive edge in any
    evaluable cohort of that system (review Q4: any cohort).
  - Per organ system, n_cohorts_positive / n_cohorts_evaluable and sign
    agreement across that system's cohorts are reported.
- Outputs per pair and analysis: n_organ_systems_evaluable,
  n_organ_systems_positive, and the ratio. A missing edge in an
  unevaluable system is **not** negative evidence.
- **Recurrent edge**: positive in ≥ 3 organ systems (the Step 1 Fig 1B
  preset), reported with n_positive / n_evaluable. Primary recurrence uses
  only primary-eligible cohorts. Sensitivity recurrence is computed and
  reported separately.
- Network: nodes are TFs, edges are recurrent edges, and edge weight is the
  number of organ systems. Modules come from `igraph::cluster_leiden`
  (objective = modularity, resolution 1, n_iterations = 10, seed 1).
  Modules with ≥ 5 TFs are reported with n nodes, n edges, and the top TFs by
  weighted degree. The full edge list is kept, so interpretation does not
  rest on the cluster picture alone (review N5). Node annotations: DBD,
  Step 1 recurrence, recurrent proliferation coupling.
- **SOX9 reference (CRC only)**: SOX9's positive edges in COADREAD and their
  recurrence in other organ systems are reported. This is descriptive, not
  a selection criterion (charter non-goal: SOX9 is not assumed to be
  pan-cancer).

### 3.8 Focused Fig 1B view (prospective, additive)

Step 1 TFs with `fig1B_primary == TRUE` (C1-up in ≥ 3 organ systems) that
are **not** recurrently proliferation-coupled (§3.4). The 358-TF Step 1
Fig 1B stays the prespecified discovery figure; this view is a separate,
labelled secondary figure. LBX2 keeps its `pipeline_sensitive` flag.

### 3.9 BRCA PAM50 sensitivity (review Q3)

BRCA only, not blocking other cohorts:

- (a) Primary statistic with PAM50 indicator covariates added (LumA
  reference; LumB, Her2, Basal; Normal-like excluded). df is reduced
  accordingly.
- (b) Primary statistic within LumA only, if n ≥ 100 with purity.
- Report the edge Jaccard of (a) and (b) against BRCA primary, and the
  FOXA1–GATA3 control under each. BRCA primary stays the value used for
  recurrence; (a) and (b) are sensitivity results only.

## 4. Coverage / validity checks

1. **Purity coverage**: per cohort, the fraction and number of tumors with
   called ABSOLUTE purity, and primary eligibility (§3.5). LAML
   coverage and its missingness pattern are checked explicitly (review Q5).
   Report the ρ between purity and proliferation score per cohort.
2. **Proliferation score check**: Spearman ρ(score, MKI67), with MKI67
   held out of the score. ≥ 0.7 is expected in every carcinoma cohort; a
   miss is investigated, not waived.
3. **Adjustment check**: E2F1–FOXM1 (both proliferation TFs) should have a
   high unadjusted ρ that is attenuated after proliferation adjustment.
   Report both values per cohort. This only shows the adjustment acts; it
   is not a result.
4. **Positive-control pairs** (preset biological sanity checks):
   FOXA1–GATA3 (BRCA), CDX1–CDX2 (COADREAD), HNF1A–HNF4A (LIHC), TP63–SOX2
   (LUSC). Expected: positive unadjusted ρ; the adjusted values are
   reported. They are not required to stay significant after adjustment. A
   miss is investigated; parameters are never changed to make a control
   pass (review N4).
5. **Permutation null** (diagnostic only, review N2): per cohort and
   analysis, 20 permutations (seeds 1–20). In each, permute every TF's
   residual vector independently across samples, recompute ρ, and count
   edges passing the §3.5 rule. Report the observed edge count and the mean
   and max permuted counts; no ratio is reported when the permuted count is
   0. This is a diagnostic, not a formal error-rate estimate.
6. **Gene mapping**: signature and TF IDs present per cohort, and the
   assertion that signature and TF IDs do not overlap (§3.4).
7. **Sensitivity**: edge overlap (Jaccard) between the primary and
   unadjusted analyses per cohort, to show how much purity and
   proliferation drive co-expression.

## 5. Expected outputs

Repo (`results/02_tf_coexpression/`, small TSV):

- `cohort_qc.tsv`: n tumors, n with purity, TF universe size, ρ(purity,
  prolif), ρ(score, MKI67), primary eligibility, observed and permuted
  (mean, max) edge counts per analysis, edge counts.
- `proliferation_signature_genes.tsv` (frozen in this PR).
- `TF_proliferation_association.tsv`: TF, cohort, analysis, ρ, q, coupled.
- `TF_proliferation_summary.tsv`: per TF and definition: n_evaluable,
  n_coupled, fraction, `recurrently_proliferation_coupled`.
- `TF_coexpression_edges.tsv`: cohort, analysis, TF_a, TF_b, ρ, q, sign;
  normal ρ, Δρ, bootstrap CI and the descriptive label where applicable.
  Edges only (passing §3.5).
- `TF_pair_recurrence.tsv`: pair, analysis, n_organ_systems_evaluable,
  n_organ_systems_positive, per-system cohort counts and sign agreement.
- `TF_network_modules.tsv`: TF, module, weighted degree, annotations;
  `TF_network_module_summary.tsv`: module, n nodes, n edges, top TFs.
- `BRCA_PAM50_sensitivity.tsv` (§3.9).
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

1. Download and verify the ABSOLUTE and PAM50 files (after user OK).
2. Run the smoke test on COADREAD (purity coverage, score check, controls,
   permutation null). Stop if a check fails unexplained.
3. Run all 16 cohorts, then recurrence, network and figures.

## 7. Open questions for review (round 4; answers in §8)

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

## 8. Review history

- **Round 4** (design rev 1, head 1d1636d): REQUEST CHANGES, 4 blocking
  and 7 non-blocking.
  - **B1**, primary statistic undefined under purity missingness → §3.5:
    three predefined analyses, primary eligibility (≥ 70% coverage and
    ≥ 50 tumors), no relabelling.
  - **B2**, Fisher z not calibrated → §3.6: dropped; matched estimand
    (proliferation-only on both sides); bootstrap CI for recurrent edges
    only; descriptive labels.
  - **B3**, signature auditability → §3.4: frozen by Ensembl ID via HGNC;
    exclusions by ID; no-overlap assertion; attenuation read as
    proliferation-associated covariance.
  - **B4**, pair-level denominators → §3.7: evaluable and positive organ
    systems; missing ≠ negative.
  - **N1–N6** → §3.2, §4.5, §3.4, §4.4, §3.7, §3.6.
  - **N7** → §1: TF–BAF is a separate module.
  - **Answers to the open questions:**
    - Q1: accept ρ ≥ 0.4 / 0.5 with BH; no top-1% rule.
    - Q2: adjusted primary, unadjusted kept in full.
    - Q3: add a BRCA PAM50 sensitivity → §3.9.
    - Q4: any cohort, with consistency info.
    - Q5: accept ABSOLUTE; check coverage, especially LAML.
- **Round 5** (rev 2, head 6a8f136): **APPROVE**, 0 blocking issues.
  - N1: the TF–score association has no proliferation-only variant → fixed
    in §3.4.
  - N2: `cohort_qc.tsv` still mentioned a permutation ratio → fixed in §5.
