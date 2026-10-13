# TF_Cancer

<One-paragraph description: what this project builds/tests and why —
the founding question(s) from `docs/PROJECT_CHARTER.md`, in plain
language.>

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

<one or two sentences + a pointer: "See docs/PROJECT_SUMMARY.md for
details." Keep this section short — it goes stale fast; the summary doc
is the source of truth.>

## Scope

<what this repo covers, in a few sentences. Link out to
docs/PROJECT_CHARTER.md's "explicit non-goals / deferred scope" for
what's intentionally excluded.>

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

Heavy compute runs on <cluster name> (see `docs/PROJECT_CHARTER.md`
§2.4). Nothing here should require running large jobs on a laptop.

## Related repos

<list any pre-existing related repos this project builds on or overlaps
with, if any, with a one-line note on the relationship — mirror
TWEAKR-OncoPlacental's README.md pattern for `REPOS/...`.>
