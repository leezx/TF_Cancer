# Worklog — TF_Cancer

Chronological log. Append one entry per merged PR (or other substantial
checkpoint — a blocker surfaced and resolved, a scope decision made).
Newest entry at the bottom (or top — pick one convention on day 1 and
keep it for the life of the project; do not switch partway through).
Do not edit past entries except to fix a factual error; this is a log,
not a living doc.

Do not treat this file as decorative. It is what a new agent (or you, in
three weeks) reads to reconstruct *why* something is the way it is
without re-deriving it from the diff.

## Entry shape

Each entry should answer, briefly:

```md
## <YYYY-MM-DD> — Step <N>: <short name> (PR #<n>, merged <sha>)

**What**: one or two sentences — what this step/PR did.

**Why**: why this was the right next step (ties back to the founding
question or a dependency from a prior step).

**How**: key method/design decisions, especially anything a later
reader would otherwise have to reverse-engineer from code (e.g. "scored
GSE231559 tumor and normal populations with separate null calibrations,
not a joint pass, because X").

**Real findings**: anything actually discovered by running real
code/compute against real data — bugs caught (and how), coverage gaps,
unexpected numbers, things that matched vs. didn't match a pre-compute
estimate. Do not omit negative or unresolved findings; they're often the
most useful entries.

**Review**: N rounds on the ChatGPT reviewer tab; note anything genuinely
notable (a reviewer catch that was real, a self-caught issue, a claim
that didn't survive independent verification, a retraction of your own
prior round's conclusion). "Clean approve, no notable findings" is a
fine, honest thing to write when that's what happened — don't pad it.

**Merged**: <sha>, after explicit user confirmation on <date>.
```

## Progress tracker

Keep a running table near the top (or in a pinned section) once the
project has more than a couple of steps:

| Step | Status | PR | Notes |
|---|---|---|---|
| 1. <name> | closed | #<n> | |
| 2. <name> | closed | #<n> | |
| 3. <name> | in progress | — | blocked on <x> |

Overall: `<current>% → <new>% (+<delta>%)` against the 100% definition
in `docs/PROJECT_CHARTER.md` Part 1 — update this whenever a step
closes, per `AGENTS.md` rule 10.

---

## 2026-10-13 — Step 1: TCGA pan-cancer TF DE (PR #1 design, PR #2 implementation + results, merged 113c0c3)

**What**: DESeq2 TF differential expression over 15 TCGA cohorts (39
contrasts: C1 TCGA tumor vs normal, C2 vs GTEx, C3 paired), with LAML as
expression only. Recurrence by organ system, QC gates and draft Figs 1–2.

**Why**: answers the first half of Q1 (are TFs recurrently up across
cancers) and produces the TF × cohort asset that Steps 2–7 build on.

**How**:
- Toil counts were recovered with `round(2^x − 1)` and validated against
  recount3 (an independent pipeline) on matched samples before the full run.
- C1 is the primary evidence class; C2 is exploratory (`cross_source_only`).
- Shrunk-LFC tiers.
- Run on the Mac mini, because the Argos `all.q` slots were fully booked.

**Real findings**:
- The recount3 gate passed with a wide margin (99–100% direction, 99–100%
  Tier 1 recovery).
- LBX2 is the only discordant TF and is flagged `pipeline_sensitive`.
- Toil's COAD gap is structured by tissue source site.
- 358 TFs are C1-up in ≥ 3 of 9 organ systems. The top of the ranking is
  proliferation TFs.
- Canonical lineage TFs are expressed but often not up versus normal (KICH
  PAX8/PAX2/HNF1B are down), so Q1 is not answered by "up" alone.
- C1 vs C2 rank correlation is low in KIRC, THCA and PRAD (0.31–0.48).
- Bug caught before compute: the Argos wrapper crashed when called with no
  arguments, and `-hold_jid` released dependent jobs even after a failure.
  This was moot after the move to local.

**Review**:
- PR #1: R1 REQUEST CHANGES (B1 confounding, B2 count validity, B3 normal
  references, B4 recurrence by evidence class), then R2 APPROVE.
- PR #2: R3 APPROVE with 9 notes. Checked against the code: N2 (recount3
  uses the Toil TPM filter) also applies to the matched set, which the
  reviewer missed. Fixed in ecae021. The final-head check approved.

**Merged**: 113c0c3, after explicit user confirmation on 2026-10-13.

