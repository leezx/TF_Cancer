# Step 1: TCGA pan-cancer TF differential expression — design

> Written and reviewed BEFORE any compute. Status: **draft, awaiting
> ChatGPT review** (charter §2.2). No DE has been run; all sample counts
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
   rows beyond the Lambert DBD family; see "not in scope").

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
5. Keep one sample per patient per group (12-char barcode). Pre-check
   finds no duplicates in any TCGA cohort; if any appear, keep the
   lexicographically first sample ID and log it.
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

**LAML has no valid normal in Toil.** Whole Blood is mostly mature
granulocytes, not myeloid progenitors. GTEx "Bone Marrow" in Toil is
`Cells - Leukemia Cell Line (Cml)`, a CML cell line, and is excluded
outright. In Step 1, LAML gets tumor expression summaries only (median
TPM, prevalence) and **no DE**. A normal hematopoietic progenitor
reference is deferred to the scRNA step (open question Q3).

### 3.3 Contrasts

Run per cohort when eligible:

| ID | Contrast | Eligible when | DESeq2 design |
|---|---|---|---|
| C1 | TCGA tumor vs TCGA normal (unpaired) | n_normal ≥ 10 | `~ group` (COADREAD: `~ project + group`) |
| C2 | TCGA tumor vs GTEx normal | a GTEx tissue is mapped | `~ group` |
| C3 | Paired tumor vs normal | n_pairs ≥ 10 | `~ patient + group` (paired samples only) |

- **Primary contrast** per cohort = C1 if eligible, otherwise C2. So C2
  is primary for PAAD, GBM and LGG, and LAML has none. It is recorded in
  `cohort_manifest.tsv` and used for all recurrence counts and figures.
- In C2, the source (TCGA vs GTEx) is completely confounded with group,
  so it cannot be corrected. C2 is a supplementary check wherever C1
  exists, and every C2-primary cohort is flagged in all outputs.
- KIRC, KIRP and KICH C2 use only 28 GTEx kidney samples. Flag it.

### 3.4 Counts and DESeq2

1. Counts: `round(pmax(2^x − 1, 0))`. RSEM expected counts are not
   integers, so rounding is required. Assert that max |2^x−1 − round| is
   small; report it in the QC table.
2. Gene filter per contrast: keep genes with count ≥ 10 in at least
   `min(n_group1, n_group2)` samples. Run DESeq2 over **all** genes that
   pass, so size factors and dispersions use the whole transcriptome,
   then subset to TFs for reporting. TFs removed by the filter are kept
   in the output with `tested = FALSE`.
3. `DESeq()` with defaults (Wald test). `results(alpha = 0.05)` uses BH
   over all tested genes after independent filtering. That padj is the
   reported FDR. A TF-only BH is **not** recomputed.
4. Report both the MLE log2FC and the apeglm-shrunk log2FC
   (`lfcShrink(type = "apeglm")`). **All thresholds use the shrunk LFC.**
5. Expression from TPM: `2^x − 0.001`, floored at 0. Report tumor median
   TPM, reference-group median TPM, and tumor prevalence = fraction of
   tumor samples with TPM ≥ 1.

### 3.5 Calls (exploratory thresholds, preset in Task-1, not biological boundaries)

A TF is **up** in a cohort if padj < 0.05, tumor median TPM > 1 and
tumor prevalence > 25%, and:

- **Tier 1**: shrunk log2FC > 1
- **Tier 2**: 0.5 < shrunk log2FC ≤ 1

Downregulated TFs are not called, but full stats are kept in
`TCGA_TF_DE_all.tsv`. With n in the hundreds, almost every gene reaches
padj < 0.05, so the effect-size and expression filters do the real work.
This is stated in the results doc.

### 3.6 Integration and recurrence

- Three TF × cohort matrices from each cohort's primary contrast:
  shrunk log2FC, padj, tumor median TPM. LAML is a column in the TPM
  matrix only.
- `TF_recurrent_upregulation.tsv`, per TF:
  - `n_cohorts_tier1`, `n_cohorts_tier1or2`.
  - `n_organ_systems_tier1or2`, so that one tissue is not counted
    several times. Organ systems: colorectum, lung (LUAD, LUSC), kidney
    (KIRC, KIRP, KICH), liver, pancreas, breast, brain (GBM, LGG),
    prostate, thyroid, stomach, head & neck. LAML is excluded.
  - `n_cohorts_c1c2_concordant`: up in the primary contrast and also
    shrunk LFC > 0.5 with padj < 0.05 in every other eligible contrast.
- Fig. 1B "top recurrent" = Tier 1/2 in ≥ 3 organ systems. Preset; the
  results doc reports how many TFs this selects.

## 4. Compute plan (Argos, charter §2.4)

- Env: `/home/zz950/softwares/miniconda3/envs/r4p3/bin/Rscript`,
  `R_LIBS_USER=/home/zz950/softwares/R_lib_4`, DESeq2 1.38.1. Record
  `sessionInfo()`.
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
- **Smoke test first**: run KICH (66 T / 25 N, smallest) for C1/C2/C3
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
4. **Positive expression controls**: label integrity, not DE
   expectations. Tumor median TPM > 1 is expected for CDX2 (COADREAD),
   NKX2-1 (LUAD, THCA), TP63 and SOX2 (LUSC, HNSC), PAX8 (KIRC, KIRP,
   THCA), HNF4A (LIHC), FOXA1 and GATA3 (BRCA), AR and NKX3-1 (PRAD),
   OLIG2 and SOX2 (GBM, LGG), SPI1 and CEBPA (LAML). A miss is
   investigated, not waived (charter §2.6).
5. **Proliferation sanity check**: FOXM1 and E2F1 are expected to be up
   in most carcinoma C1 contrasts. This only confirms the direction is
   correct; it is not a result.
6. **Contrast concordance**: Spearman ρ of shrunk TF LFC between C1 and
   C2, and between C1 and C3, per cohort. Low ρ (preset flag at ρ < 0.5)
   is reported and discussed, not used to drop a cohort.
7. **COADREAD**: the READ-only C1 is compared with COADREAD C1 by TF LFC ρ.

## 6. Expected outputs

Repo (`results/01_tcga_tf_de/`, small TSV):

- `TF_reference.tsv`
- `cohort_manifest.tsv` (n per group, eligible contrasts, primary
  contrast, GTEx tissue, flags)
- `sample_manifest.tsv` (sample, cohort, group, source, patient)
- `TCGA_TF_DE_all.tsv` (long format: cohort, contrast, TF, baseMean,
  log2FC_MLE, lfcSE, log2FC_shrunk, pvalue, padj, tumor_median_TPM,
  ref_median_TPM, tumor_prevalence, n_tumor, n_ref, ref_source, tested)
- `TCGA_TF_DE_significant.tsv` (with tier)
- `TF_x_cohort_log2FC.tsv`, `TF_x_cohort_padj.tsv`,
  `TF_x_cohort_tumorTPM.tsv`
- `TF_recurrent_upregulation.tsv`
- `qc/` (§5 gate tables)
- Draft figures: `figures/Fig1A_all_up_TFs_heatmap.pdf`,
  `figures/Fig1B_recurrent_TFs_heatmap.pdf`,
  `figures/Fig2_TF_UpSet.pdf` and `figures/Fig2_TF_dotplot.pdf`

Results doc: `docs/STEP1_TCGA_TF_DE_RESULTS.md`, written after compute.

## 7. Open questions for review

- **Q1** Is Toil acceptable as the primary source given its smaller COAD
  tumor count (288 vs GDC ~470)? Alternative: re-query GDC STAR-Counts
  with full barcodes for C1/C3 and keep Toil only for C2.
- **Q2** GTEx tissue choices: Colon - Transverse for COADREAD; the two
  cortex regions for GBM and LGG; **Esophagus - Mucosa as the HNSC
  proxy**.
- **Q3** LAML: is "no DE in Step 1" right, or should Whole Blood C2 run
  with a strong flag?
- **Q4** Should calls threshold on the apeglm-shrunk LFC or the MLE LFC?
- **Q5** BRCA subtypes: which PAM50 source should be used, and should
  subtype-level DE be part of Step 1 or come after it?
- **Q6** Is "≥ 3 organ systems" a reasonable preset for the Fig. 1B
  recurrent set?
