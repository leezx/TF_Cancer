# Step 2: within-cancer TF co-expression and proliferation association — results (v0.1)

## Run details

- Design: `docs/STEP2_TF_COEXPRESSION_DESIGN.md`, approved in review R5
  (PR #3). The thresholds are those of the design, not tuned after
  looking.
- Code: `scripts/02_tf_coexpression/` (`common2.R`, `00`–`04`).
- Runs, all on the Mac mini on 2026-10-13:

| Step | Time | Peak RSS |
|---|---|---|
| Prepare | 2 min | 2.6 GB |
| Co-expression | 7 min | 0.9 GB |
| Integration, including 9 × 1,000 bootstraps | 60 min | 2.5 GB |
| Plots | < 1 min | — |

- Logs: `logs/step2_full.log`, `logs/step2_integrate.log`,
  `logs/step2_plot.log`.
- Repo tables and figures are in `results/02_tf_coexpression/`.
- The full per-cohort matrices, edge lists and large tables are in
  `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/2026-10-13_tf_cancer_step2_coexpr_v0.1/`
  (`prep/`, `coexpr/`, `report/`).

An edge is correlated mRNA abundance across tumors of one cancer. It is not
binding, interaction or co-occupancy. After ABSOLUTE purity adjustment,
bulk co-expression can still reflect the non-malignant composition (§4).

## 1. Checks (§4)

| Check | Result |
|---|---|
| Purity coverage | 88–97% in 14 cohorts. KIRC (62%) and LAML (64%) are not primary-eligible; they contribute prolif-only and unadjusted results only |
| ρ(purity, proliferation score) | −0.16 (KIRC) to 0.60 (GBM) |
| Score vs MKI67 (held out) | ≥ 0.7 in 13 of 16 cohorts. **Misses: COADREAD 0.65, HNSC 0.67, LUSC 0.54** (see below) |
| Adjustment check, E2F1–FOXM1 | Attenuated in every cohort. For example: COADREAD 0.53 → 0.16; LUAD 0.74 → 0.23; LUSC 0.47 → 0.17 (unadjusted → primary) |
| Positive controls | BRCA FOXA1–GATA3: primary 0.64 (edge). LUSC TP63–SOX2: 0.57 (edge). LIHC HNF1A–HNF4A: unadjusted 0.46 (edge) → primary 0.33 (no edge). COADREAD CDX1–CDX2: unadjusted 0.31 → primary 0.29 (no edge in any analysis) |
| Permutation null (20 per cohort and analysis) | Mean ≤ 0.05 edges, max 1 (in 5 of 46 cohort × analysis runs), against 8,130–63,674 observed positive primary edges per cohort |
| Signature / TF ID overlap | 0 (asserted). 301 signature genes are present in every cohort |

**MKI67 misses.**
- MKI67 is a single noisy marker in these cohorts. For example, in
  COADREAD TOP2A–MKI67 is 0.54, while the score vs TOP2A is 0.71.
- The score is unchanged. The design says "investigated, not waived".
- It is a deviation from the expected ≥ 0.7. It is recorded here and in
  §6.

**CDX1–CDX2 and HNF1A–HNF4A.**
- Both are positive, as expected. They miss the |ρ| ≥ 0.4 edge rule after
  adjustment.
- The design does not require controls to stay significant after
  adjustment, and no parameter was changed.

## 2. Proliferation association (§3.4)

- **Proliferation coupling is concentrated in a few TFs.**
  - Per cohort, 4 (LUSC) to 20 TFs are coupled (purity-adjusted ρ ≥ 0.5,
    q < 0.05).
  - 16 TFs are recurrently coupled (purity-adjusted; 25 unadjusted). They
    are FOXM1, DNMT1, MYBL2, ZNF367, E2F1, CENPA, E2F8, E2F2 and E2F7,
    among others.
- **Of the 358 Step 1 Fig 1B TFs, only 9 are recurrently
  proliferation-coupled.**
  - Their median purity-adjusted ρ with the score is 0.03.
  - 28 have a median ρ ≥ 0.3.
- So the proliferation TFs lead the Step 1 ranking, but most of the 358
  recurrently-up TFs do not track proliferation across tumors at this
  threshold.
- The focused Fig 1B view therefore keeps **349 TFs**
  (`TF_focused_fig1B.tsv`, `figures/Fig1B_focused_nonproliferative_TFs.pdf`).
- The proliferation heatmap
  (`figures/FigS_TF_proliferation_coupling_heatmap.pdf`) shows a band of
  moderate positive ρ below the 0.5 threshold. "Not coupled" means "below
  the preset threshold", not "no association".

## 3. TF–TF co-expression and recurrence (§3.5, §3.7)

**The network is dense.**
- Positive primary edges per cohort range from 8,130 (LUSC) to 63,674
  (THCA). Across all analyses: 406,343 primary positive edges, 542,834
  prolif-only and 589,749 unadjusted.
- Edge Jaccard of primary vs unadjusted is 0.50 (STAD) to 0.78 (THCA).
  Purity and proliferation adjustment changes a third to a half of the
  edges.

**Recurrent pairs** (positive in ≥ 3 organ systems):
- 48,976 primary, 61,435 prolif-only, 67,108 unadjusted.
- The median primary pair is positive in 4 organ systems; 651 pairs are
  positive in 11.
- At the frozen thresholds (|ρ| ≥ 0.4, ≥ 3 organ systems), the network
  has 1,154 TFs and 48,976 edges. Fig 3B is a hairball.
- The cluster picture is not informative on its own. The module table
  and the full edge lists carry the result (review N5).

**Modules** (Leiden, preset parameters; 6 with ≥ 5 TFs):

| Module | TFs | Content | Reading |
|---|---|---|---|
| 4 | 293 | 271 ZNF/KZFP; 69% on chr19 | Clustered KRAB-ZNF genes: most likely shared genomic and chromatin regulation of the chr19 clusters, not a lineage program |
| 1 | 283 | ZNF787, MBD3, ZNF524…; 23% chr19; **140 of the Fig 1B TFs** | Broadly expressed, GC-rich TFs; a global axis |
| 2 | 326 | ASH1L, BPTF, KMT2A, MGA, ZKSCAN8 | Large chromatin regulators; a global "transcriptional capacity" axis |
| 3 | 228 | TCF4, ZEB1/2, PRRX1 (stromal); ETS1, FLI1, ERG, HEYL, EPAS1 (endothelial); SPI1, IKZF1, LYL1, IRF1 (immune) | **Microenvironment composition that survives ABSOLUTE purity adjustment** |
| 5, 6 | 10, 6 | HOXA cluster, HOXC cluster | Co-located gene clusters |

**No module corresponds to a tumor-lineage TF program.**
- Lineage pairs exist as edges. FOXA1–GATA3 and TP63–SOX2 are edges.
- But they are not separated by recurrence across organ systems, because
  lineage programs are by definition not shared across organs.
- So recurrence across organ systems selects global, genomic-cluster and
  microenvironment axes. This is the main Step 2 result for the network
  question.
- Lineage-level structure belongs to the lineage-specificity axis (next
  design, out of scope here).

**SOX9 reference** (CRC; descriptive):
- **SOX9 has no edge in COADREAD primary.** Its strongest partners are
  HES1 0.36, MNX1 0.35, HOXB8 0.33 and JUN 0.31, all below the 0.4 edge
  rule. `SOX9_COADREAD_reference_edges.tsv` is therefore empty.
- In other organ systems, SOX9's recurrent primary partners are ELF3 and
  HNF1B (4 systems each), and ETV4, ETV6 and HES1 (3 each). These are
  epithelial and Notch TFs.

## 4. Tumor vs normal (§3.6; descriptive)

- 9 cohorts with ≥ 40 normals were analysed, with a bootstrap Δρ CI for
  recurrent pairs.
- `tumor_enriched_coexpression` (Δρ CI excludes 0 and normal ρ < 0.2)
  applies to 1,446 (THCA) to 8,118 (PRAD) recurrent pairs per cohort:
  COADREAD 5,893, LUAD 2,850, LUSC 1,862, KIRC 5,592, LIHC 3,193,
  BRCA 2,737 and HNSC 1,830.
- With 44–113 normals, these labels are descriptive only and are not
  ranked.

## 5. BRCA PAM50 sensitivity (§3.9)

| Analysis | n | Positive edges | Jaccard vs BRCA primary | FOXA1–GATA3 |
|---|---|---|---|---|
| Primary | 1,017 | 28,215 | — | 0.64 |
| + PAM50 covariates | 985 | 25,085 | 0.77 | 0.38 |
| LumA only | 528 | 41,028 | 0.60 | 0.40 |

PAM50 explains part of the FOXA1–GATA3 correlation, which falls from 0.64
to 0.38. Most of the BRCA network is stable. BRCA primary stays the value
used for recurrence.

## 6. Deviations from the design

1. **The large tables are in DATA, not the repo.** These are
   `TF_pair_recurrence`, `TF_coexpression_edges_recurrent_pairs` and
   `TF_tumor_vs_normal_recurrent_pairs`, together 190 MB. They are in
   `.../report/*.tsv.gz` (md5 of the decompressed files checked against
   the originals).
   - The design's `TF_coexpression_edges.tsv` is the per-cohort
     `coexpr/<cohort>_edges.tsv.gz` in DATA.
   - The repo keeps the summaries, modules, controls and QC (4.9 MB).
   - Reason: the repo rule allows small text only. The edge counts were
     far larger than the design expected.
2. **MKI67 check misses** in 3 cohorts (§1). These were investigated; the
   score is unchanged.
3. **DATA matrix naming.** The design named `coexpr/<cohort>_{partial,unadjusted}.rds`;
   the run wrote `coexpr/<cohort>.rds`, holding all analyses plus the
   normal matrix.
4. **FigS row order.** Rows are ordered by mean ρ rather than clustered.
   TFs absent from some cohorts leave NA cells, so hierarchical
   clustering fails. This is display only.

## 7. What this answers

- **Observed.**
  - Proliferation drives a small, well-defined set of TFs (E2F/MYB/FOXM1).
  - Most Step 1 recurrent-up TFs are not proliferation-coupled at
    ρ ≥ 0.5.
  - Recurrent TF co-expression across organ systems is dominated by
    KZFP genomic clusters, global chromatin-regulator axes and
    microenvironment programs.
- **Supported.**
  - Adjusting for purity and proliferation removes the cell-cycle
    covariance (E2F1–FOXM1). The controls behave.
  - A stromal, endothelial and immune module persists after purity
    adjustment. Bulk composition is not fully removed by ABSOLUTE purity.
- **Not established.**
  - Any regulatory or physical TF–TF relation.
  - A shared pan-cancer lineage TF module. The recurrence design cannot
    find one, by construction.
  - SOX9 partners in CRC: SOX9 has no edge at the frozen threshold.
- **For the next design** (lineage axis + Step 3):
  - Lineage programs need a within-organ or lineage-specificity
    statistic, not cross-organ recurrence.
  - The microenvironment module argues for a cell-type deconvolution
    covariate or a single-cell check (Step 5).
