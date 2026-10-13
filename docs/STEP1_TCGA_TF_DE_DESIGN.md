# Step 1: TCGA pan-cancer TF differential expression — design

> Written and reviewed BEFORE any compute. Status: **rev 2, addressing
> review round 1 (REQUEST CHANGES, B1–B4)** (charter §2.2); see §8. No DE has been run; all sample counts
> below are pre-compute tallies from the phenotype file, used as the
> expected values for the coverage gate in §5.

## 1. Scope

Covers Task-1 Stage 1 items 1–7 (charter Q1, first half):

1. Freeze the TF universe (`TF_reference.tsv`, Lambert fields only).
2. Freeze the cohort and sample manifests.
3. Per-cohort tumor-vs-normal DESeq2 over all genes; report every TF.
4. Build the TF × cohort matrices (log2FC, padj, tumor median TPM).
5. Count recurrent upregulation per TF.
6. Draft Fig. 1A/1B heatmaps and the Fig. 2 UpSet/dot plot (no annotation
   rows beyond the Lambert DBD family; see "not in scope"). Step 1 figures
   are labeled **"TF expression discovery"**, not "lineage-defining TF
   discovery". Lineage claims need the Step 3 annotation.

**Not in scope for Step 1** (planned later steps, see
`docs/PROJECT_SUMMARY.md`):

- Within-cancer TF co-expression and purity-adjusted partial correlation
  (Task-1 item 8) → Step 2.
- Curated functional annotation (lineage identity / developmental /
  pioneer evidence / proliferation / stress) → Step 3. Deliberately
  *after* DE so the DE is run over the whole TF universe without a
  pre-selected list (Task-1 Module 1 discovery-bias note).
- BRCA molecular subtypes (PAM50 is not in the Toil phenotype) → added
  as a Step 1 follow-up once a PAM50 source is chosen (open question Q5).
- DepMap, scRNA, ATAC → Steps 4–6.

**Deviation from Task-1**: Task-1 says to prefer GDC harmonized
STAR counts. This step uses **UCSC Toil RSEM expected_count** instead,
because (a) Toil processed TCGA and GTEx with the same pipeline, which
the TCGA-vs-GTEx contrast needs; (b) the existing GDC per-project
Rdata on Argos lost their sample-type barcodes (see the hot data list
§2.2), so they cannot separate tumor from normal. Cost: Toil COAD has
fewer primary tumors (288) than GDC (~470); see open question Q1.

## 2. Inputs

All paths are relative to `DATA/` (local SSD and Argos `/home/zz950/DATA/`
are md5-identical; see `DATA/HotData/TF_Cancer_PanCancer_TF_Atlas.md`).

| Input | Path | Note |
|---|---|---|
| Counts | `bulkRNAseq/UCSC_Toil_TCGA_GTEx/raw/TcgaTargetGtex_gene_expected_count.gz` | log2(expected_count+1), 60,498 genes × 19,109 samples, md5 `87b9b2c4…` |
| TPM | Local `bulkRNAseq/UCSC_Toil_TCGA_GTEx/raw/TcgaTargetGtex_rsem_gene_tpm.gz`; Argos `bulkRNAseq/public_bulk_DB/UCSC_Xena/` | log2(TPM+0.001); 19,131 samples (count ⊂ TPM) |
| Phenotype | `.../TcgaTargetGTEX_phenotype.txt.gz` | md5 `25cc6f2b…` |
| Gene map | `.../gencode.v23.annotation.gene.probemap` | GENCODE v23; 60,498 IDs, no duplicate IDs after stripping versions |
| TF universe | `2.PROJECTS/Cancer_Context_Dependency/TF_executors/data/raw/HumanTF_paper/DatabaseExtract_v_1.01.csv` | Lambert et al. 2018, `Is TF? == Yes` → 1,639 |

## 3. Method / contract

### 3.1 TF reference (`TF_reference.tsv`)

- Universe = Lambert v1.01 rows with `Is TF? == "Yes"` (1,639).
- Key = Lambert `Ensembl ID`, matched to the bare (version-stripped)
  GENCODE v23 ID. Assert one-to-one with no duplicate bare IDs. Pre-check:
  1,637 match; **DUX1 and DUX3 have no Ensembl ID** (they are not in
  GENCODE) and stay in the table with `in_expression_matrix = FALSE`.
- Columns: `ensembl_id, hgnc_symbol (Lambert), gencode_v23_id,
  gencode_v23_symbol, DBD, binding_mode, motif_status,
  in_expression_matrix`. A symbol mismatch between Lambert and GENCODE v23
  is recorded but never used for matching.
- No functional categories in Step 1 (deferred to Step 3).

### 3.2 Cohorts (16 analysis cohorts from 17 TCGA projects)

Sample selection, in this order:

1. `_study == "TCGA"`. This must be explicit: `detailed_category`
   "acute myeloid leukemia" also contains TARGET pediatric AML, 200 vs
   173 TCGA.
2. Tumor = `_sample_type == "Primary Tumor"`, or for LAML
   `"Primary Blood Derived Cancer - Peripheral Blood"`. Exclude Recurrent,
   Metastatic, and Additional–New Primary.
3. TCGA normal = `"Solid Tissue Normal"`.
4. Keep samples present in both the count and TPM matrices.
5. Keep one sample per patient per group (12-char barcode). Toil IDs are
   sample-level (`TCGA-XX-XXXX-01`) with no vial or aliquot field, and the
   pre-check finds no duplicates in any TCGA cohort. **If a duplicate
   appears, the run stops and the case is investigated.** No automatic
   tie-break.
6. COAD + READ form one cohort, `COADREAD`, with `project` (COAD/READ)
   as a covariate. A READ-only run is a sensitivity check (§5).

GTEx normal mapping (`_study == "GTEX"`, by `detailed_category`):

| Cohort | Tumor (pre-check) | TCGA normal | Pairs | GTEx tissue | GTEx n |
|---|---:|---:|---:|---|---:|
| COADREAD | 379 | 51 | 32 | Colon - Transverse | 166 |
| LUAD | 513 | 59 | 58 | Lung | 288 |
| LUSC | 498 | 50 | 50 | Lung | 288 |
| KIRC | 530 | 72 | 72 | Kidney - Cortex | 28 |
| KIRP | 288 | 32 | 32 | Kidney - Cortex | 28 |
| KICH | 66 | 25 | 25 | Kidney - Cortex | 28 |
| LIHC | 369 | 50 | 50 | Liver | 110 |
| LAML | 173 | 0 | 0 | **none valid** (see below) | — |
| PAAD | 178 | 4 | 4 | Pancreas | 167 |
| BRCA | 1,092 | 113 | 112 | Breast - Mammary Tissue | 178 |
| GBM | 153 | 5 | 0 | Brain - Cortex + Brain - Frontal Cortex (Ba9) | 206 |
| LGG | 509 | 0 | 0 | Brain - Cortex + Brain - Frontal Cortex (Ba9) | 206 |
| PRAD | 494 | 51 | 51 | Prostate | 100 |
| THCA | 504 | 59 | 59 | Thyroid | 279 |
| STAD | 413 | 36 | 33 | Stomach | 175 |
| HNSC | 518 | 44 | 43 | Esophagus - Mucosa (squamous epithelium proxy; flagged) | 273 |

Per-cohort normal-reference quality goes into `cohort_manifest.tsv`
(review B3). `C2 quality` applies only to the GTEx contrast:

| Cohort | C1 available | C2 quality | Lineage match / composition risk | High-confidence recurrence eligible |
|---|---|---|---|---|
| COADREAD | yes | moderate | colon region mismatch (rectum vs transverse); check READ-only sensitivity | yes (C1) |
| LUAD, LUSC | yes | moderate | bulk lung (alveolar and airway mixture) | yes (C1) |
| KIRC, KIRP, KICH | yes | low (n=28) | cortex; KICH originates in distal nephron | yes (C1) |
| LIHC | yes | moderate | — | yes (C1) |
| PAAD | **no** (4 normals) | low | whole pancreas: acinar-dominated; PDAC is ductal; high composition risk | **no**: exploratory only |
| BRCA | yes | moderate | adipose-rich mammary tissue; composition risk | yes (C1) |
| GBM, LGG | **no** | low | bulk cortex is neuron-rich, not a glial-lineage match | **no**: exploratory only |
| PRAD, THCA, STAD | yes | moderate | — | yes (C1) |
| HNSC | yes | low (proxy) | esophageal mucosa is not oral/laryngeal epithelium | yes (C1); C2 is a proxy sensitivity check only |
| LAML | no | none | — | no (DE_status = `no_valid_normal`) |

**LAML has no valid normal in Toil.** Whole Blood is mostly mature
granulocytes, not myeloid progenitors. GTEx "Bone Marrow" in Toil is
`Cells - Leukemia Cell Line (Cml)`, a CML cell line, and is excluded
outright. In Step 1, LAML gets tumor expression summaries only (median
TPM, prevalence) and **no DE**: DE fields are NA with
`DE_status = "no_valid_normal"`, and LAML is excluded from every
recurrence denominator. The later reference will be purified normal
CD34+ progenitors and myeloid stages with donor-level replication. A normal hematopoietic progenitor
reference is deferred to the scRNA step (open question Q3).

### 3.3 Contrasts

Run per cohort when eligible:

| ID | Contrast | Eligible when | DESeq2 design |
|---|---|---|---|
| C1 | TCGA tumor vs TCGA normal (unpaired) | n_normal ≥ 10 | `~ group` (COADREAD: `~ project + group`) |
| C2 | TCGA tumor vs GTEx normal | a GTEx tissue is mapped | `~ group` |
| C3 | Paired tumor vs normal | n_pairs ≥ 10 | `~ patient + group` (paired samples only) |

Evidence classes (review B1). Every DE row carries an `evidence_class`:

- **C1 = high-confidence discovery**, the only evidence used for the
  main ranking. Available in 12 cohorts / 9 organ systems.
- **C2 = exploratory, `cross_source_only`**, in every cohort. Group is
  fully confounded with source (TCGA vs GTEx collection, processing and
  RNA quality), so a C2-only call is never high-confidence without
  independent validation (the scRNA step or an external cohort). PAAD,
  GBM and LGG have only C2, so in Step 1 they contribute exploratory
  evidence only.
- **C3 = supporting sensitivity evidence** for C1 (same patients). It is
  not an independent cohort and not independent replication.
- KIRC, KIRP and KICH C2 use only 28 GTEx kidney samples. Flag it.

### 3.4 Counts and DESeq2

0. **Input-validation gate (review B2). It must pass before the full array.**
   - Scale is confirmed from the authoritative source metadata
     (`TcgaTargetGtex_gene_expected_count.json` on the Xena toil hub):
     `"unit": "log2(expected_count+1)"`, version 2016-09-03. Values are
     stored with 4 decimals.
   - Quantify the reconstruction error |(2^x−1) − round(2^x−1)| by
     expression bin (0, 1–10, 10–100, 100–1e3, >1e3 counts), per gene
     and per sample. Also report library sizes and zero fractions per
     cohort and group.
   - Known limitation, stated in the results doc: RSEM expected counts
     are estimated fragment assignments, and rounding does not restore
     quantification uncertainty. DESeq2 on rounded expected counts is
     common practice (tximport `countsFromAbundance="no"` style), not a
     guarantee.
   - **Independent-pipeline sensitivity check.** The reviewer asked for
     GDC STAR; we use **recount3** instead because the GDC per-project
     objects available to us lost tumor/normal barcodes. recount3 is an
     independent STAR-based (Monorail) count pipeline with full barcodes
     and sample types.
     - (a) COADREAD C1: recount3 data are already local
       (`DATA/bulkRNAseq/recount3_TCGA_COAD_READ_GTEx_COLON/`, COAD 546
       + READ 177 + GTEx colon 822).
     - (b) LUAD C1: needs the recount3 LUAD gene_sums, a new small
       download at implementation time, to be confirmed with the user.
     - Same DESeq2 contract on both pipelines. Compare TF results on:
       direction agreement, Tier 1/2 call overlap, and the identity of
       discordant TFs. Spearman ρ is reported but is not sufficient on
       its own.
     - **Preset acceptance** (judgment-based thresholds, fixed now):
       (i) direction agreement ≥ 95% among TFs with |shrunk LFC| > 1 in
       either pipeline; (ii) ≥ 80% of Toil Tier 1 calls are Tier 1/2 in
       recount3; (iii) every discordant Tier 1 TF is listed and examined.
       A failure blocks the full run until it is explained.
   - Q1 coverage audit: list the recount3 COAD/READ primary tumors that
     are absent from Toil (546 vs 288 COAD). Compare available sample
     annotations (project, sample type, tissue source site, stage when
     present) between Toil-present and Toil-absent samples, and report
     any imbalance.
1. Counts: `round(pmax(2^x − 1, 0))`.
2. Gene filter per contrast: keep genes with count ≥ 10 in at least
   `min(n_group1, n_group2)` samples. Report filter sensitivity: the
   number of TFs dropped per contrast, and whether the dropped TFs reach
   tumor median TPM ≥ 1, since highly imbalanced cohorts can lose
   low-abundance TFs. Run DESeq2 over **all** genes that
   pass, so size factors and dispersions use the whole transcriptome,
   then subset to TFs for reporting. TFs removed by the filter are kept
   in the output with `tested = FALSE`.
3. `DESeq()` with defaults (Wald test). `results(alpha = 0.05)` uses BH
   over all tested genes after independent filtering. That padj is the
   reported FDR. A TF-only BH is **not** recomputed. Raw p-values and
   n_tested_genes per contrast are kept.
4. Report both the MLE log2FC and the apeglm-shrunk log2FC
   (`lfcShrink(type = "apeglm")`). **All thresholds use the shrunk LFC.**
5. Expression from TPM: `2^x − 0.001`, floored at 0. Report tumor median
   TPM, reference-group median TPM, and tumor prevalence = fraction of
   tumor samples with TPM ≥ 1.

### 3.5 Calls (exploratory thresholds, preset in Task-1, not biological boundaries)

A TF is **up** in a contrast if padj < 0.05, tumor median TPM ≥ 1 and
tumor prevalence ≥ 25% (inclusive boundaries), and:

- **Tier 1**: shrunk log2FC ≥ 1
- **Tier 2**: 0.5 ≤ shrunk log2FC < 1

padj tests LFC ≠ 0. The shrunk-LFC cut is an additional effect-size
selection criterion, **not** a formal test that the true LFC exceeds the
threshold. Threshold-aware testing (`lfcThreshold`) is deferred unless
that claim is needed.

Downregulated TFs are not called, but full stats are kept in
`TCGA_TF_DE_all.tsv`. With n in the hundreds, almost every gene reaches
padj < 0.05, so the effect-size and expression filters do the real work.
This is stated in the results doc.

### 3.6 Integration and recurrence (review B4)

- TF × cohort matrices (shrunk log2FC, padj, tumor median TPM), one set
  per evidence class: `C1` (12 cohorts) and `C2_exploratory` (15
  cohorts). LAML is a column in the TPM matrix only.
- Organ systems: colorectum, lung (LUAD, LUSC), kidney (KIRC, KIRP, KICH),
  liver, pancreas, breast, brain (GBM, LGG), prostate, thyroid, stomach,
  head & neck. LAML is excluded. An organ system counts as up if at least
  one of its cohorts is up.
- `TF_recurrent_upregulation.tsv`, per TF:

| Metric | Definition |
|---|---|
| `n_organ_systems_C1_up` | organ systems with a C1 Tier 1/2 call. **Main discovery ranking** |
| `n_organ_systems_C1_tier1` | same, Tier 1 only |
| `n_organ_systems_C1_C3_supported` | C1 call with same-direction C3 (shrunk LFC > 0) where C3 is eligible |
| `n_organ_systems_C2_only_up` | up only via C2 (no C1 call, or C1 not eligible) |
| `n_organ_systems_exploratory_up` | up via C1 or C2. Secondary view |
| `n_organ_systems_discordant` | an eligible C1 and C2 disagree in direction (one up, the other shrunk LFC < 0) |
| `n_organ_systems_eligible_C1` / `_exploratory` | denominators, reported with each count |

- Fig. 1B: **primary** = TFs with `n_organ_systems_C1_up ≥ 3`;
  **supplementary** = `n_organ_systems_exploratory_up ≥ 3`. The ≥ 3 cut
  is a visualization criterion, **not** a biological definition of
  pan-cancer lineage regulation, and it is never used to drop
  lineage-specific TFs from later prioritization. The full table keeps
  every TF from 0 upward.

## 4. Compute plan (Argos, charter §2.4)

- Env: `/home/zz950/softwares/miniconda3/envs/r4p3/bin/Rscript`,
  `R_LIBS_USER=/home/zz950/softwares/R_lib_4`, DESeq2 1.38.1, apeglm
  1.20.0, ashr 2.2-54, limma 3.54.0, data.table 1.14.10. Freeze the full
  versions in `sessionInfo()`. Add a machine-readable `config.yaml` with
  input md5s, thresholds and command lines, plus `code_commit.txt`.
- `scripts/01_tcga_tf_de/`:
  - `01_build_tf_reference.R`: local or login node, light.
  - `02_build_manifests.R`: cohort and sample manifests, light.
  - `03_extract_cohort_matrices.R`: qsub, about 32 GB. Reads the 1.3 GB
    count and TPM files once and writes one `.rds` of counts and TPM per
    cohort to `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/processed/v0.1/`.
  - `04_run_deseq2.R <cohort> <contrast>`: qsub array, one task per
    eligible cohort × contrast (~40 tasks).
  - `05_integrate.R`: matrices and recurrence.
  - `06_plot_step1.R`: draft Figs. 1A/1B/2.
- Order: (1) the §3.4 input-validation gate (reconstruction QC + recount3
  COADREAD; LUAD after its download is confirmed), (2) the smoke test,
  (3) the full array.
- **Smoke test**: run KICH (66 T / 25 N, smallest) for C1/C2/C3
  end to end. Check shapes, the §5 gate and the positive controls, then
  submit the full array.
- Heavy outputs (DESeq2 objects, full all-gene tables):
  `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/<YYYY-MM-DD>_tf_cancer_step1_tcga_de_v0.1/{de,report}/`.
  This is one root per AGENTS rule 9 and includes config, commit,
  sessionInfo and logs. Only small TF-level TSVs and figures go into the
  repo under `results/01_tcga_tf_de/`.

## 5. Coverage / validity checks

1. **Sample-count gate**: realized n_tumor, n_normal, n_pairs and n_GTEx
   per cohort must equal the §3.2 table exactly. Any difference stops
   the run and is investigated before going on.
2. **TF coverage**: 1,637/1,639 TFs present. Report n_tested per contrast
   and list TFs that are untested in every cohort.
3. **Count back-transform**: per-sample library size and the
   non-integer residual; flag outliers above 3 MAD.
4. **Tissue-identity controls** (sample-label and lineage representation
   only; they do not show tumor-specific upregulation). Tumor median TPM > 1 is expected for CDX2 (COADREAD),
   NKX2-1 (LUAD, THCA), TP63 and SOX2 (LUSC, HNSC), PAX8 (KIRC, KIRP,
   THCA), HNF4A (LIHC), FOXA1 and GATA3 (BRCA), AR and NKX3-1 (PRAD),
   OLIG2 and SOX2 (GBM, LGG), SPI1 and CEBPA (LAML). A miss is
   investigated, not waived (charter §2.6).
5. **Proliferation controls**, kept separate from identity controls:
   FOXM1 and E2F1 are expected to be up
   in most carcinoma C1 contrasts. This only confirms the direction is
   correct; it is not a result.
6. **Contrast concordance**: for C1 vs C2 and C1 vs C3 per cohort, report
   Spearman ρ of shrunk TF LFC, direction agreement, and Tier 1/2 call
   overlap. Correlation alone can hide disagreement near the threshold. Low ρ (preset flag at ρ < 0.5)
   is reported and discussed, not used to drop a cohort.
7. **COADREAD**: the READ-only C1 is compared with COADREAD C1 by TF LFC ρ.

## 6. Expected outputs

Repo (`results/01_tcga_tf_de/`, small TSV):

- `TF_reference.tsv`
- `cohort_manifest.tsv` (n per group, eligible contrasts, GTEx tissue,
  normal_reference_quality, lineage_match_status, composition_risk,
  eligible_for_high_confidence_recurrence, flags)
- `sample_manifest.tsv` (sample, cohort, group, source, patient)
- `TCGA_TF_DE_all.tsv` (long format: cohort, organ_system, contrast,
  evidence_class, cross_source_only, primary_for_ranking, TF, baseMean,
  log2FC_MLE, lfcSE, log2FC_shrunk, shrinkage_method, pvalue, padj,
  n_tested_genes, tumor_median_TPM, ref_median_TPM, tumor_prevalence,
  n_tumor, n_ref, ref_source, normal_reference_quality, tested, DE_status)
- `qc/input_validation/` (reconstruction-error bins, library sizes, zero
  fractions, recount3-vs-Toil concordance, Toil-absent sample audit)
- `TCGA_TF_DE_significant.tsv` (with tier)
- `TF_x_cohort_log2FC.tsv`, `TF_x_cohort_padj.tsv`,
  `TF_x_cohort_tumorTPM.tsv`
- `TF_recurrent_upregulation.tsv`
- `qc/` (§5 gate tables)
- Draft figures: `figures/Fig1A_all_up_TFs_heatmap.pdf`,
  `figures/Fig1B_recurrent_TFs_heatmap.pdf`,
  `figures/Fig2_TF_UpSet.pdf` and `figures/Fig2_TF_dotplot.pdf`

Results doc: `docs/STEP1_TCGA_TF_DE_RESULTS.md`, written after compute.

## 7. Decisions from review round 1 (were open questions Q1–Q6)

- **Q1** Toil is the primary source. Conditions: the recount3 sensitivity
  check and the Toil-absent sample audit (§3.4 step 0).
- **Q2** The GTEx choices stand with the quality flags in §3.2. Brain
  cortex is exploratory only; esophageal mucosa is an explicitly labeled
  HNSC proxy.
- **Q3** LAML is expression-only, `DE_status = no_valid_normal`, and
  excluded from denominators.
- **Q4** Calls use the apeglm-shrunk LFC; MLE is also kept. This is not a
  threshold test (§3.5).
- **Q5** BRCA PAM50 comes from established TCGA annotations (PanCancer
  Atlas / TCGA BRCA subtype tables) with barcode-level verification. It
  is a Step 1 follow-up, not part of the main pan-BRCA run.
- **Q6** ≥ 3 organ systems is a visualization threshold only, applied
  separately to C1 and exploratory evidence (§3.6).

## 8. Review history

- **Round 1** (head 6710258, ChatGPT conversation "设计泛癌TF分析方案"):
  REQUEST CHANGES.
  - B1 cross-source confounding: accepted → evidence classes (§3.3).
  - B2 count validity: accepted. Scale verified from the Xena JSON
    metadata. GDC sensitivity replaced by recount3 (GDC barcodes not
    available to us; recount3 is an independent STAR-based pipeline),
    with preset acceptance criteria (§3.4 step 0).
  - B3 normal-reference quality: accepted → quality table (§3.2).
  - B4 recurrence by evidence quality: accepted → metric table (§3.6).
  - N1–N9 adopted. N1 is modified: there are no vial/aliquot fields in
    Toil, so a duplicate stops the run instead of a ranked tie-break.
- **Round 2** (head 5993fd6): **APPROVE**, with zero remaining blockers.
  The recount3 substitution and the preset thresholds (≥ 95% direction,
  ≥ 80% Tier 1 recovery) are accepted. These implementation clarifications
  are binding on the validation stage:
  1. Convert recount3 coverage to counts with `recount3::transform_counts`.
     Compare on a harmonized gene universe (bare Ensembl IDs).
  2. Two separate comparisons for COADREAD and LUAD:
     - (a) **matched-sample**, the primary technical validation: the same
       TCGA samples in both pipelines, with patient de-duplication
       verified;
     - (b) **full-cohort**, a secondary coverage sensitivity check.
  3. Tier 1 recovery is reported with two denominators: TFs evaluable in
     recount3, and all Toil Tier 1 TFs (non-evaluable = unrecovered).
     Missingness is reported.
  4. A threshold failure triggers investigation (sample coverage,
     annotation, low counts, normalization, shrinkage, pipeline), not
     automatic rejection. Thresholds are never changed after seeing
     results. Any exception is documented along with its effect on the
     candidate list.
  - Approved order: input QC + recount3 COADREAD → recount3 LUAD (after
    download confirmation) → KICH smoke test → full 16-cohort array.
