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

<first real entry goes here>
