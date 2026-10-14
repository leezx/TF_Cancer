# Step 3 design — TF functional annotation and lineage specificity

Status: **draft v0.1, not reviewed**. No compute until approved.

## 1. Question

Step 1 found that "recurrently up" TFs are mostly proliferation TFs, and
that canonical lineage TFs are expressed in their tumors but often not up
versus normal. Step 2 found no pan-cancer lineage module, because a
recurrence criterion cannot see lineage-restricted programs.

Step 3 asks two questions:

1. Which TFs are lineage-restricted, and which lineage does each belong to?
2. Are lineage-restricted TFs up in the cancer of their own lineage more
   often than other TFs are? This answers Q1 directly.

All 1,639 Lambert TFs are annotated. Annotations rank and describe TFs;
they never remove a TF (Task-1: do not pre-select known pioneers).

## 2. Inputs (all local, md5 already recorded in each dataset's `link.md`)

| Input | Source |
|---|---|
| TF universe, DBD | `results/01_tcga_tf_de/TF_reference.tsv` (Lambert v1.01) |
| Adult normal tissue expression | `DATA/1.Databases/GTEx_v11_median_tpm/raw/` (68 tissues) |
| Cell-type expression | `DATA/1.Databases/HPA_single_cell_type/raw/rna_single_cell_type.tsv` |
| Tumor expression and DE | Step 1: `TF_x_cohort_tumorTPM.tsv`, `TCGA_TF_DE_all.tsv`, `cohort_manifest.tsv` |
| Cohort → GTEx tissue map | `cohort_manifest.tsv` column `gtex_tissues` (frozen in Step 1) |
| Proliferation coupling | Step 2 `TF_proliferation_association` table (PR #8) |
| Gene sets | MSigDB Hallmark v2026.1 (local `DATA/genesets/`), GO BP via `msigdbr` |

## 3. Annotation axes

### 3.1 Lineage specificity (data-driven; the main new statistic)

- **Normal tissue tau.** Tau (Yanai 2005) over GTEx tissues on
  log2(median TPM + 1). Tissues are first collapsed to organs by taking
  the maximum over sub-tissues (for example, the brain regions become
  one organ), so that organs with many sub-tissues do not dominate.
- **Lineage-restricted** = tau ≥ 0.8 and maximum organ TPM ≥ 5. Both
  cutoffs are preset (0.8 is the conventional tau cutoff for
  tissue-specific genes) and must not change after results are seen.
- **Lineage organ** = the organ with the highest TPM.
- **Cell-type tau** is computed the same way over the HPA single-cell
  types, with the top cell type recorded. It separates epithelial
  lineage TFs from stromal and immune TFs, which carry over from Step 2
  module M3.
- **Tumor-type tau** is computed over the 16 TCGA cohorts' median tumor
  TPM. It is descriptive only and is not used in the test in §4.

### 3.2 Functional categories (gene-set membership; not mutually exclusive)

| Category | Rule (frozen before compute) |
|---|---|
| Proliferation | recurrently proliferation-coupled in Step 2, or in HALLMARK_E2F_TARGETS / G2M_CHECKPOINT |
| Stress-responsive | in HALLMARK UNFOLDED_PROTEIN_RESPONSE, HYPOXIA, TNFA_SIGNALING_VIA_NFKB, or P53_PATHWAY |
| Developmental | GO BP "pattern specification process" (GO:0007389) or "cell fate commitment" (GO:0045165) |
| Stemness-associated | GO BP "stem cell population maintenance" (GO:0019827) |
| Lineage identity | lineage-restricted by §3.1 |

### 3.3 Pioneer evidence (literature; small and explicit)

- This is a curated table: TF, evidence level, PMID, one-line basis.
- The evidence levels are:
  - **experimental**: nucleosome binding in vitro, or opening closed
    chromatin shown directly;
  - **inferred**: claimed in a review or from indirect data.
- Every PMID is checked on PubMed before entry. A TF with no entry is
  "not curated", not "non-pioneer".

## 4. Pre-specified test (answers Q1)

For each TCGA cohort with a GTEx match (LAML has none and is excluded):

- **Matched lineage TFs** = lineage-restricted TFs whose lineage organ is
  the cohort's own GTEx organ.
- Compare the fraction of matched lineage TFs that are C1-up in that
  cohort with the fraction among all other expressed TFs. Report the
  odds ratio and Fisher p for each cohort, then BH across cohorts.
- Run a second readout in which "up" is replaced by "expressed in the
  tumor" (median TPM ≥ 5). This separates "kept" from "gained".

Expected interpretations, written before the result is seen:

- If lineage TFs are kept but not up, the "recurrently up" framing is
  wrong for lineage TFs. In that case Step 7 should rank by retained
  expression plus dependency, not by DE.
- If lineage TFs are enriched among up TFs, the Q1 hypothesis is
  supported for that cohort.

## 5. Positive controls (gate)

| Control | Expected lineage organ |
|---|---|
| CDX2 | colon |
| NKX2-1 | lung or thyroid |
| GATA3 | breast |
| AR | prostate |
| HNF1A or HNF4A | liver |
| PAX8 | thyroid or kidney |

All six must be lineage-restricted and assigned to the expected organ.
If fewer than five pass, stop and review the tau construction before any
results are read.

The negative control is E2F1 and FOXM1, which must not be
lineage-restricted.

## 6. Outputs

- `results/03_tf_annotation/TF_annotation.tsv`: one row per TF, with
  DBD, tau values, lineage organ, top cell type, category flags, and
  pioneer evidence.
- `results/03_tf_annotation/lineage_TF_enrichment.tsv`: the §4 test
  results per cohort.
- `results/03_tf_annotation/pioneer_curation.tsv`: the §3.3 table.
- Figure: per-cohort log2FC of matched lineage TFs versus other TFs, and
  a tau versus recurrence scatter.
- Large intermediates, if any, go to
  `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/<date>_tf_cancer_step3_annotation_v0.1/`.

## 7. Limits

- GTEx is bulk tissue. A TF expressed in a minority epithelial cell type
  can show low tau. HPA cell-type tau is reported next to it to flag
  these cases.
- The Step 1 HNSC normal is esophagus mucosa, so HNSC matching is
  approximate. It is reported but marked as such.
- Annotation is association with tissue expression, not proof of a
  lineage-determining function.
