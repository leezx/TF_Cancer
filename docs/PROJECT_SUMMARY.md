# Project Summary — TF_Cancer

> Staged deliverable report. Update at every major milestone. Detail lives
> in `Worklog.md`.

## Founding question(s)

Pan-cancer Lineage-defining TF Atlas (`Tasks/Task-1.md`):

1. **Q1**: Are lineage-defining TFs recurrently upregulated across cancers?
2. **Q2**: Are they expressed and co-expressed within malignant cells, not
   just a reflection of bulk tissue composition?
3. **Q3**: Which TFs are plausible candidates for BAF-dependent enhancer
   remodeling, and in which tractable models?

## Status

`25%` complete against the charter definition of 100% (5% → 25%, +20%,
on 2026-10-13: Step 1 closed).

The input data are inventoried and complete: Toil TCGA+GTEx counts and
TPM, DepMap 26Q1, and the Lambert TFs, md5-verified on local and Argos.
Step 1 is closed: design PR #1 and results PR #2 were merged after 3
review rounds, and all gates passed. Step 2 (within-cancer TF
co-expression + proliferation axis): design approved (PR #3); computed
2026-10-13, results in review (`docs/STEP2_TF_COEXPRESSION_RESULTS.md`).
TF–BAF module (outside the 100%): v0.1 merged (PR #4, #5); v0.2
(organoids, DepMap context, background percentile) design merged (PR #6),
results in review (PR #7).

## Roadmap and progress (100% = all rows done)

| Step | Workstream | Weight | Status | Blocker | Answers |
|---|---|---:|---|---|---|
| 0 | Data inventory + gap fill (hot data list) | 5% | **done** | — | infra |
| 1 | TCGA pan-cancer TF DE atlas (Figs 1–2) | 20% | **done** (PR #1, #2) | — | Q1 |
| 2 | Within-cancer TF co-expression network + proliferation axis (Fig 3) | 10% | computed; results in review | results review + merge | Q1 |
| 3 | Curated TF functional annotation (lineage/pioneer/etc.) | 10% | not started | — | Q1/Q3 |
| 4 | DepMap patient-vs-cell-line concordance and model table (Fig 5) | 10% | not started | Step 1 | Q2/Q3 |
| 5 | scRNA malignant-vs-normal pseudobulk validation (Fig 4) | 20% | not started | dataset selection | Q2 |
| 6 | Chromatin evidence (TCGA ATAC) | 10% | not started | Step 1 | Q3 |
| 7 | Tiered candidate prioritization, Fig 6, final results doc | 15% | not started | Steps 1–6 | Q3 |

## What's been answered

Step 1 (`docs/STEP1_TCGA_TF_DE_RESULTS.md`, frozen v0.1):

- The Toil count input is reliable for DE: recount3 agrees on direction
  (99–100%) and Tier 1 calls (99–100%) in COADREAD and LUAD.
- 917 TFs are up (C1) in ≥ 1 of 9 organ systems and 358 in ≥ 3. The most
  recurrent are proliferation TFs (E2F1, FOXM1, MYBL2, …). Q1 is **not yet
  answered**: canonical lineage TFs are expressed in their tumors but often
  not up vs normal, so "recurrently up" is not the same as
  "lineage-defining".

Step 2 (`docs/STEP2_TF_COEXPRESSION_RESULTS.md`, in review):

- Proliferation coupling is concentrated in about 16 TFs (E2F/MYB/FOXM1).
  Only 9 of the 358 Step 1 recurrent-up TFs are recurrently coupled, so a
  focused 349-TF Fig 1B view is available.
- Recurrent cross-organ TF co-expression is dominated by chr19 KZFP
  clusters, global chromatin-regulator axes, and a stromal, endothelial
  and immune module that survives purity adjustment. No pan-cancer
  lineage module was found.
- Lineage structure needs a lineage-specificity statistic (next design).
- SOX9 has no edge in COADREAD at the frozen threshold.

## What's still open

- Everything in Steps 1–7.
- Scope to confirm with the PI: the TAD hypothesis (A vs B, Task-1 §八)
  is deferred, and organoids are optional (charter non-goals).

## Steps closed so far

| Step | What | Design doc | Results doc | PR(s) |
|---|---|---|---|---|
| 1 | TCGA TF DE atlas | `docs/STEP1_TCGA_TF_DE_DESIGN.md` | `docs/STEP1_TCGA_TF_DE_RESULTS.md` | #1, #2 |

## Key locked artifacts

- Hot data list: `DATA/HotData/TF_Cancer_PanCancer_TF_Atlas.md`
- Step 1 v0.1: `results/01_tcga_tf_de/`; DESeq2 objects under
  `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/2026-10-13_tf_cancer_step1_tcga_de_v0.1/`

## Scope boundaries

See `docs/PROJECT_CHARTER.md` Part 1 "Explicit non-goals".
