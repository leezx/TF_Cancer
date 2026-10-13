# Step 1: TCGA pan-cancer TF differential expression — results

- Design: `docs/STEP1_TCGA_TF_DE_DESIGN.md` (approved, PR #1, 2 review rounds)
- analysis_id: `2026-10-13_tf_cancer_step1_tcga_de_v0.1`
- Code: branch `step1-tcga-tf-de-impl`, scripts `scripts/01_tcga_tf_de/`; run
  records (commit, command, sessionInfo, config) under the result root
  `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/<analysis_id>/report/run_records/`
  and `processed/v0.1/run_records/`.
- Small tables and draft figures: `results/01_tcga_tf_de/`. DESeq2 objects and
  all-gene tables: `DATA/.../result/<analysis_id>/de/<source>_<set>/`.
- Run: Mac mini, 2026-10-13, `scripts/01_tcga_tf_de/run_step1_local.sh`
  `prep → validate → smoke → full` (approved order). Wall time ≈ 1 h 20 min;
  peak RSS 3.5 GB (prep), 5.6 GB (full).

Everything below is **TF expression discovery** in bulk tissue. Nothing here
shows that a TF is lineage-defining, tumor-cell-intrinsic, or functional.

## 1. Gates (all passed, thresholds as preset in the design)

| Gate | Criterion | Result |
|---|---|---|
| Sample counts (§5.1) | realized n = §3.2 table | 64/64 checks PASS |
| TF coverage (§5.2) | 1,637/1,639 in matrix | DUX1, DUX3 absent (no Ensembl ID) |
| Count back-transform (§3.4.0) | report reconstruction error | max \|y − round(y)\| = 0.5 in every bin, i.e. compatible with integer rounding (not by itself proof of exact reconstruction); max relative error 0.05% for counts > 1,000. The recount3 concordance below is the stronger evidence |
| recount3 sensitivity (§3.4.0, matched samples, gated) | direction ≥ 95% (\|LFC\| > 1); Toil Tier 1 → recount3 Tier 1/2 ≥ 80% | **PASS** — COADREAD 99.4% / 98.9% (ρ 0.987); LUAD 100% / 100% (ρ 0.993) |
| recount3, full cohort (secondary, not gated) | — | COADREAD 99.2% / 98.9% (ρ 0.967); LUAD 100% / 100% (ρ 0.993) |
| Identity controls (§5.4) | tumor median TPM > 1 | 21/21 pass (e.g. CDX2 135, NKX2-1 LUAD 77, NKX3-1 430, PAX8 THCA 342 TPM) |
| Proliferation controls (§5.5) | FOXM1, E2F1 up in most carcinoma C1 | up in 12/12 C1 cohorts except THCA FOXM1 (LFC 0.44, below Tier 2); PRAD E2F1 Tier 2 |

recount3 discordant Toil Tier 1 TFs (§3.4.0 iii): only **LBX2** (COADREAD;
Toil shrunk LFC 1.62, recount3 −0.02). recount3 baseMean is 4× higher
(122 vs 30) and flat. LBX2 overlaps the antisense gene LBX2-AS1; a likely
explanation is that recount3's coverage-based gene sums assign antisense
reads to LBX2 while RSEM does not. Not verified (strandedness, gene models and coverage not yet inspected);
LBX2 is flagged `pipeline_sensitive = TRUE` in
`TF_recurrent_upregulation.tsv` and must be checked before any downstream
prioritization.

Toil-absent audit (§3.4.0 Q1): Toil lacks 171/458 COAD and 74/166 READ
primary tumors that recount3 has, but only 3/516 LUAD. **The COAD gap is
site-structured**: tissue source site AA contributes 149 absent vs 23
present samples. Stage distributions look broadly similar
(`qc/input_validation/toil_absent_audit_*.tsv`). Toil COADREAD is therefore
a non-random subset by collection site; the recount3 full-cohort check
(ρ 0.967, direction 99.2%) suggests count-based TF DE is robust to it, but
it does not exclude site-associated biology. A CRC-central candidate should
be re-checked in the full recount3 cohort and an independent CRC dataset.

Two qualifications of the recount3 comparison: (1) the expression filters
(tumor median TPM, prevalence) use **Toil TPM for both pipelines**, in the
matched and the full set, so the comparison tests count-based DE (padj,
shrunk LFC), not an independent recount3 expression classification; in the
full set the TPM gate uses the Toil production samples. (2) Duplicate
recount3 aliquots were resolved by keeping the largest library per patient
× group (a deterministic rule, not a quality score); retained and excluded
columns are listed in
`qc/input_validation/recount3_aliquot_selection_{COADREAD,LUAD}.tsv`.

## 2. Calls per contrast (Tier 1 / Tier 2 up)

| Cohort | C1 (TCGA T vs TCGA N) | C2 (TCGA T vs GTEx) | C3 (paired) |
|---|---|---|---|
| BRCA | 93/154 | 153/186 | 67/129 |
| COADREAD | 93/167 | 127/155 | 75/148 |
| HNSC | 72/82 | 141/165 | 63/108 |
| KICH | 113/197 | 194/204 | 93/183 |
| KIRC | 99/144 | 87/189 | 69/95 |
| KIRP | 76/116 | 57/130 | 48/92 |
| LIHC | 72/165 | 46/89 | 34/128 |
| LUAD | 93/170 | 152/199 | 51/132 |
| LUSC | 172/164 | 246/203 | 118/170 |
| PRAD | 33/92 | 74/137 | 27/80 |
| STAD | 73/124 | 183/191 | 78/110 |
| THCA | 40/74 | 79/151 | 37/69 |
| GBM | — | 250/234 | — |
| LGG | — | 197/255 | — |
| PAAD | — | 182/123 | — |
| LAML | expression only (`no_valid_normal`) | | |

1,215–1,519 TFs tested per contrast. Filter sensitivity: dropped TFs with
tumor median TPM ≥ 1 are 0 in all C1 contrasts; nonzero only in GBM C2
(28), LUSC C3 (10), PAAD C2 (7), STAD C3 (3), COADREAD C3 / LUAD C3 (1).
55 TFs are untested everywhere (`qc/tf_never_tested.txt`; germline/
testis/embryonic TFs such as DUX4, NANOGNB, FIGLA, HSFY1).

As stated in the design, with n in the hundreds almost every gene reaches
padj < 0.05; the expression and effect-size filters do the selecting.

## 3. Contrast concordance (§5.6)

C1 vs C3 (paired, same patients): ρ 0.82–0.94, direction agreement
94–99% for all 12 cohorts. C1 vs C2 (cross-source GTEx reference): ρ
0.31–0.79; **flagged ρ < 0.5 in KIRC (0.31), THCA (0.39), PRAD (0.48)**;
direction agreement 73–95%. This confirms review B1: C2 mixes tumor
biology with source/batch and reference-tissue differences, so C2-only
calls stay `exploratory_cross_source`. READ-only C1 vs COADREAD C1: ρ 0.94,
up-call Jaccard 0.72.

## 4. Recurrence (C1, 9 eligible organ systems)

917 TFs are Tier 1/2 up in ≥ 1 organ system. Distribution of organ systems
up: 1: 353, 2: 206, 3: 140, 4: 96, 5: 62, 6: 29, 7: 15, 8: 11, 9: 5.
**358 TFs meet the preset Fig 1B threshold (≥ 3 organ systems)**; 754 meet
the exploratory threshold (C1 or C2).

The most recurrent TFs (up in 8–9/9 organ systems, nearly all
directionally supported by paired C3, i.e. C3 shrunk LFC > 0, not C3
significance) are dominated by cell-cycle / proliferation TFs: E2F1,
NFE2L3, HES6, HES4, SOX4, FOXM1, FOXS1, CBX2, MYBL2, SOX12, MXD3, ARID3A,
E2F3. This is the expected proliferation signal (§5.5), not lineage
biology.

## 5. Interpretation and limits

- **Observed:** a reproducible (Toil ↔ recount3) TF × cohort DE atlas for
  15 TCGA cohorts; a shared pan-cancer up-regulated TF set led by
  proliferation TFs.
- **Not established:** that any recurrently up TF is lineage-defining.
  Canonical lineage TFs are expressed in their tumors (identity controls
  21/21) but are often **not up vs normal** — e.g. KICH PAX8, PAX2, HNF1B
  are down in C1. Ranking by "up vs normal" therefore under-represents
  TFs that are maintained, rather than induced, in tumors.
- The ≥ 3-organ-system cut (preset, not changed) selects 358 TFs and does
  not discriminate on its own. The 358-TF Fig 1B stays as the prespecified
  discovery output; any focused figure will use criteria defined
  prospectively in the next design.
- Bulk tissue: composition (stroma, immune, normal epithelium) can drive
  calls; tumor-cell specificity is Step 5.

### Open points carried forward (not post-hoc changes to Step 1)

Placement agreed in review round 3; Step 1 calls and thresholds stay frozen.

1. **Step 2:** proliferation association per cohort (predefined cell-cycle
   signature, purity-aware) as an annotation and ranking axis, not an
   exclusion rule; keep the Step 1 ranking alongside it.
2. **Next design, integrated with Step 3:** a lineage-specificity axis.
   Tumor-type-restricted expression is not the same as lineage identity,
   which needs normal-lineage references and literature. Keep separate
   evidence axes; no premature composite score.
3. LBX2 pipeline sensitivity; COAD site-structured Toil gap (§1).

## 6. Deviations from the design

| Item | Design | Actual | Impact |
|---|---|---|---|
| Compute | Argos SGE array (§4) | Mac mini, sequential, 4 threads | none on results; Argos `all.q` fully booked; charter compute updated |
| DESeq2 | 1.38.1 (Argos env) | 1.50.2 (local R 4.5.2) | versions in sessionInfo |
| recount3 duplicates | not specified | keep the largest library per patient × group | COADREAD 720 → 675 columns, LUAD 599 → 575 |
| `lfcThreshold` | deferred | not used | as designed |

## 7. Review

- Round 3 (results, PR #2, same ChatGPT review conversation): **APPROVE**,
  0 blocking, 9 non-blocking notes. Verified independently against the code:
  N2 is real and also applies to the matched set (§1).
- Addressed in PR #2: N1 aliquot audit (`05b_recount3_aliquot_audit.R`,
  loader moved to `recount3_loader.R`, logic unchanged: 720 → 675 and
  599 → 575 reproduced); N2, N4, N6, N7, N8 wording above; N3
  `pipeline_sensitive` column; N9 smoke PASS marker
  (`qc/smoke_status.txt`) checked by the `full` stage. N5: no cohort will be
  rerun with a different DESeq2 version. Integration was rerun after N3;
  all other tables are unchanged.

## 8. Outputs

`results/01_tcga_tf_de/`: `TF_reference.tsv`, `cohort_manifest.tsv`,
`sample_manifest.tsv`, `deseq2_tasks.tsv`, `TCGA_TF_DE_all.tsv`,
`TCGA_TF_DE_significant.tsv`, `TF_x_cohort_{log2FC,padj}_{C1,C2}.tsv`,
`TF_x_cohort_tumorTPM.tsv`, `TF_recurrent_upregulation.tsv` (with
`pipeline_sensitive`), `qc/` (incl. `smoke_status.txt`)
(gates, concordance, controls, coverage, READ-only, input_validation/),
`figures/` (Fig1A, Fig1B primary, Fig1B supplementary C2, Fig2 UpSet,
Fig2 dot plot; drafts).
