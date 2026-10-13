# TF_Cancer

Pan-cancer lineage-defining transcription factor (TF) atlas. Asks whether
lineage-defining TFs are recurrently upregulated across cancers (TCGA vs
normal), whether that holds within malignant cells (scRNA pseudobulk), and
which TFs are plausible candidates for BAF-dependent enhancer remodeling in
tractable models (DepMap, chromatin evidence).

**Start here:**
- [`docs/PROJECT_CHARTER.md`](docs/PROJECT_CHARTER.md) — the operating
  rules this project follows (step lifecycle, data separation, compute,
  review discipline). Read this before doing any work.
- [`docs/PROJECT_SUMMARY.md`](docs/PROJECT_SUMMARY.md) — structured
  "what was done, how, why, what's been answered" summary. Read this
  next for the current state of the project.
- [`Worklog.md`](Worklog.md) — full chronological log (every review
  round, every bug caught and fixed, the progress tracker). Read this for
  detailed history or to resume active work.

## Current status

Data inventory done; Step 1 (TCGA TF DE) design under review. See
`docs/PROJECT_SUMMARY.md` for details.

## Scope

Bulk TCGA/GTEx DE, TF co-expression, DepMap concordance, scRNA validation,
and candidate prioritization. See `docs/PROJECT_CHARTER.md` Part 1 for
non-goals (wet-lab, TAD/Hi-C, organoids optional).

## Data

**No large data files live in this repo.** All raw/processed/result
data lives under
`/Volumes/Stelligen_SSD/Stelligen/DATA/<data_type>/<dataset_id>/`,
following the workspace-wide convention in `DATA/README.md`. This repo
only holds:

- `datasets/<dataset_id>/dataset.md` — manifest/metadata for each dataset
  (source URLs, sizes, what was downloaded and why, what wasn't)
- `scripts/<step>/` — analysis code, one subdirectory per pipeline step
- `results/<step>/` — real compute outputs (gene lists, calibration
  tables, small summary tables, audit docs), mirroring `scripts/`
- `docs/` — the charter, design docs, results docs, and this summary

See `datasets/` for the datasets acquired so far.

**Hard rule:** downloaded and generated data must be written directly under
`/Volumes/Stelligen_SSD/Stelligen/DATA/`; the `PR/TF_Cancer/` tree is
metadata/code/docs only. The eight hard data-placement rules are defined in
`docs/PROJECT_CHARTER.md` §2.3 and copied from
`PR/_TEMPLATE/PROJECT_CHARTER_TEMPLATE.md`.

## Compute

Heavy compute runs on Argos (see `docs/PROJECT_CHARTER.md`
§2.4). Nothing here should require running large jobs on a laptop.

## Related repos

- Input hot data list: `DATA/HotData/TF_Cancer_PanCancer_TF_Atlas.md`
- TF reference source: `DATA/2.PROJECTS/Cancer_Context_Dependency/TF_executors/`
