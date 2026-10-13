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

`5%` complete against the charter definition of 100% (0% → 5%, +5%, on
2026-10-13).

The input data are inventoried and complete: Toil TCGA+GTEx counts and
TPM, DepMap 26Q1, and the Lambert TFs, md5-verified on local and Argos.
The Step 1 design was approved after 2 review rounds (PR #1). No analysis has run.
Infrastructure is ready; scientific progress is still 0.

## Roadmap and progress (100% = all rows done)

| Step | Workstream | Weight | Status | Blocker | Answers |
|---|---|---:|---|---|---|
| 0 | Data inventory + gap fill (hot data list) | 5% | **done** | — | infra |
| 1 | TCGA pan-cancer TF DE atlas (Figs 1–2) | 20% | design APPROVED (PR #1, 2 review rounds) | merge confirmation; recount3 LUAD download | Q1 |
| 2 | Within-cancer TF co-expression network (Fig 3) | 10% | not started | Step 1 | Q1 |
| 3 | Curated TF functional annotation (lineage/pioneer/etc.) | 10% | not started | — | Q1/Q3 |
| 4 | DepMap patient-vs-cell-line concordance and model table (Fig 5) | 10% | not started | Step 1 | Q2/Q3 |
| 5 | scRNA malignant-vs-normal pseudobulk validation (Fig 4) | 20% | not started | dataset selection | Q2 |
| 6 | Chromatin evidence (TCGA ATAC) | 10% | not started | Step 1 | Q3 |
| 7 | Tiered candidate prioritization, Fig 6, final results doc | 15% | not started | Steps 1–6 | Q3 |

## What's been answered

Nothing yet.

## What's still open

- Everything in Steps 1–7.
- Scope to confirm with the PI: the TAD hypothesis (A vs B, Task-1 §八)
  is deferred, and organoids are optional (charter non-goals).

## Steps closed so far

| Step | What | Design doc | Results doc | PR(s) |
|---|---|---|---|---|
| — | — | — | — | — |

## Key locked artifacts

- Hot data list: `DATA/HotData/TF_Cancer_PanCancer_TF_Atlas.md`

## Scope boundaries

See `docs/PROJECT_CHARTER.md` Part 1 "Explicit non-goals".
