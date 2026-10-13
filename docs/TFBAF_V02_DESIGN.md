# TF–BAF module v0.2 design: organoids and DepMap context

Status: design approved in round 8 (APPROVE, 0 blocking; notes folded in
§9). No v0.2 result has been computed. The v0.1
results (`docs/TFBAF_RESULTS.md`, PR #5) are frozen: their ρ, q, calls,
BAF list, thresholds and figure-selection rule do not change.

## 1. Why v0.2

The task brief asks for three things that v0.1 does not cover:

1. **Organoid dependency (Dataset C).** v0.1 deferred it (design Q1,
   round 6). Since then, organoid CRISPR screens have been found locally.
2. **DepMap context.** The brief asks for:
   - TF and BAF dependency variance;
   - confidence intervals with ρ kept separate from significance;
   - BAF dependency prevalence;
   - a five-way dependency pattern. It separates "co-dependent" from
     "both dependent but uncorrelated".
3. **Background reference.** In the round 7 review this was added as a
   non-blocking v0.2 item. Where does a BAF member's ρ rank among all
   genes' ρ with the same TF?

## 2. Stage 0 inventory: organoids

Source: the tumour-derived organoid biobank CRISPR screens (DepMap 3D,
2026), in
`DATA/1.Databases/DepMap3D_202608/01_tumour_derived_organoid_biobank/MOESM3_Supplementary_Tables.zip`.

- **Table 6** gives the per-screen gene-level LFC, in long format:
  `sample_ID gene ensembl_id HGNC_ID library LFC`.
- **Table 5** gives the authors' `is_depleted` call, with the same keys.
- **Table 2, sheet `clinical_data`** gives `primary_tumour_type`,
  ICD-O-3 morphology, `individual_ID` and `organoid_derived_from_metastasis`.

Checked on the raw files (read only, nothing written to `raw/`):

| Item | Value |
|---|---|
| Screens / organoids | 178 / 162 (16 organoids screened in both libraries) |
| Libraries | `minLib` 144 screens, `v1_1` 34 |
| Tumour types | Colorectal 85, Oesophageal 59, Ovarian 11, Pancreatic 4, Gastric 3 |
| Colorectal | 97 screens, 83 individuals; morphology 8140/3 89, 8480/3 5, 8144/3 1, **8263/0 tubulovillous adenoma 2 screens**; 27 screens from metastasis-derived organoids |
| Genes | all 5 TFs and all 28 frozen BAF members present in every screen |
| CRC SOX9 | median LFC −0.45; `is_depleted` in 34% of CRC screens |
| CRC ACTL6A | median LFC −1.15; depleted in 92% of screens (common-essential, as expected) |

**Feasibility per pair:**

| Pair | Organoid analysis |
|---|---|
| CRC SOX9 | **feasible**: about 81 carcinoma individuals |
| AML SPI1 | **unavailable** (no haematological organoids) |
| Prostate FOXA1 | **unavailable** (no prostate organoids) |
| Breast GATA3 | **unavailable** (no breast organoids) |
| Melanoma SOX10 | **unavailable** (no melanoma organoids) |

Unavailable pairs are reported as unavailable. Organoid RNA expression
does not replace dependency, as the brief requires.

## 3. Organoid analysis (CRC SOX9 only)

### 3.1 Unit and eligibility

The unit is one screen per individual.

1. Keep colorectal organoids with malignant morphology (ICD-O behaviour
   code /3). This drops the 8263/0 adenomas.
2. If an organoid has two screens, use the `minLib` screen (the majority
   library).
3. If an individual has more than one organoid, keep the
   primary-tumour-derived one first. Then take the lexicographically first
   `sample_ID`.

The minimum n is 10, as in v0.1.

The QC table reports how many included screens are `minLib` and how many
are `v1_1`, and the LFC distribution by library (review N1).

### 3.2 Statistics

- Spearman ρ(SOX9 LFC, BAF LFC) over individuals, for each of the 28
  members.
- BH across the 28 members.
- v0.1 call thresholds: strong ρ ≥ 0.3, moderate ρ ≥ 0.2, both with
  q < 0.05.
- Bootstrap 95% CI: 1,000 resamples of individuals, percentile method,
  seed 1.
- Dependency prevalence is the fraction of screens with the authors'
  `is_depleted` = 1. The −0.5 threshold is not used here. LFC is not on
  the Chronos gene-effect scale, so −0.5 means something different.
- Organoid values are **never merged** with the DepMap 2D gene effects.
  They are a separate evidence column.

### 3.3 Screen-depth confounder (sensitivity)

**Concern.** Raw LFC is not corrected for screen strength. A screen with
deeper overall depletion makes every essential gene more negative at
once. That can create positive "co-dependency" between any two
dependencies, including SOX9 and the BAF members that are essential.

**Pre-set sensitivity: a partial Spearman controlling for screen depth.**

- *Depth* is the median LFC, within a screen, over a reference set of
  common-essential genes.
- The reference set is the genes with `is_depleted` = 1 in ≥ 90% of all
  178 screens, all tumour types.
- The 33 TF and BAF genes are excluded from the reference set.
- The reference set is defined from the depletion calls alone, before
  any correlation is looked at.

**Reading the result.** If adjustment reverses a call or attenuates ρ by
more than half, the call is flagged `depth_sensitive`. The same rule was
used for the purity flag.

- Raw ρ, adjusted ρ, and the change (size and direction) are all reported,
  not only the flag (review N2).
- Depth may also capture real global dependency biology. A
  `depth_sensitive` association is not counted as robust organoid
  co-dependency.

### 3.4 Other sensitivities

1. `minLib` screens only.
2. Primary-tumour-derived organoids only (no metastasis-derived ones).

### 3.5 Positive control

ρ(SMARCB1, SMARCE1) and ρ(ARID1A, SMARCB1) over all 178 screens are
expected to be positive, as in v0.1 §4.4. This checks the data handling.
It is not a result.

## 4. DepMap context additions (2D, v0.1 data)

All of these are added columns. v0.1 values are not recomputed or
replaced.

1. **Variance**: SD of the gene effect, for the TF and for each BAF
   member, within the lineage.
2. **BAF prevalence**: fraction of models with gene effect < −0.5. This
   is the same rule as the v0.1 TF prevalence.
3. **Bootstrap 95% CI for ρ**: 1,000 resamples of models, percentile
   method, seed 1. Prostate stays `underpowered`.
4. **Dependency pattern** (descriptive). A gene counts as *dependent in
   the lineage* if its prevalence is ≥ 25%. For organoids, prevalence is
   the `is_depleted` fraction. Each TF–BAF row then falls in one of five
   patterns:

| Pattern | TF dependent | BAF dependent | Call |
|---|---|---|---|
| both dependent, correlated | yes | yes | ≥ moderate, positive |
| both dependent, no positive co-dependency | yes | yes | not ≥ moderate positive |
| TF-specific dependency | yes | no | any |
| BAF-specific dependency | no | yes | any |
| neither dependent | no | no | any |

- A significant negative ρ (ρ ≤ −0.2, q < 0.05) keeps its pattern and is
  flagged `negative` (review N3).
- Only "both dependent, correlated" is described as co-dependency
  consistent with joint function.
- The 25% cut-off is descriptive. It is fixed here and not tuned.

## 5. Background reference percentile

This puts each BAF member's ρ in context against all genes. It is
descriptive and does not test anything.

- Two percentiles are reported:
  - *signed*: how prominent a positive association is;
  - *absolute* (|ρ|): how unusual its strength is (for example, AML
    negative ρ). (Review N4.)
- The TF itself is excluded.
- Reported per pair: the number of background genes, plus how many were
  dropped for missing values or zero variance (review N5).
- Prostate DepMap (n = 10) percentiles are exploratory only.

- **TCGA**: per pair, Spearman ρ of the TF against every protein-coding
  gene in the Toil TPM matrix with median TPM ≥ 1 in the cohort, using the
  v0.1 samples (review Q4).
  - Protein-coding is defined as an approved HGNC record with a UniProt
    ID, matched by Ensembl gene ID (version stripped), from
    `DATA/1.Databases/HGNC_gene_id_mapping/raw/`.
  - The local GENCODE v23 probemap has no biotype column, so this proxy is
    fixed here instead of downloading new data.
  Report each BAF member's percentile among those genes.
- **DepMap**: the same, against every gene in the 26Q1 CRISPR matrix,
  over the v0.1 lineage models.
- **Organoids**: the same, against every gene in Table 6, over the §3.1
  individuals.

The question this answers: is SOX9's ρ with SMARCA4 (TCGA 0.33) unusual,
or typical of any expressed gene? The v0.1 TCGA correlations were broad.
For example, FOXA1's median ρ with BAF members was 0.53.

## 6. Figures

- **Fig TFBAF-2 (summary dot plot)**:
  - It gains a third dot, organoid, for the CRC column. In the other four
    columns the third position is grey and labelled `unavailable`.
  - The colour, size and outline rules are unchanged.
  - The caption states the organoid n and the depth-sensitivity flags.
- **Fig TFBAF-1 (scatter panels)**:
  - It gains an organoid panel for the CRC row.
  - It shows the three CRC members already selected by the v0.1 rule. No
    new selection is made.
- **New FigS TFBAF-3**: the dependency-pattern grid (pair × member,
  coloured by the five patterns), with 2D and organoid shown separately.

## 7. Outputs

Repo, under `results/tfbaf/` (small TSV + PDF):

- `Organoid_TF_BAF_codependency.tsv`: n, ρ, CI, p, q, call, prevalence,
  depth-adjusted ρ, flags.
- `Organoid_TF_BAF_sensitivity.tsv`.
- `DepMap_TF_BAF_context.tsv`: SDs, BAF prevalence, CI, pattern.
- `TF_BAF_background_percentile.tsv`.
- `qc/organoid_inventory.tsv` and `qc/organoid_positive_control.tsv`.
- The updated figures.

DATA, under one new root, next to v0.1:

- The root is
  `DATA/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/2026-10-13_tf_cancer_tfbaf_v0.2/`.
- It has the module folders `organoid/`, `depmap_context/` and
  `background/`.
- They hold the extracted organoid subset, the unit selection, and the
  run records.
- The zip is read in a stream. Nothing is written into the
  `DepMap3D_202608` folder.
- That folder is not in the standard `raw/processed/result` layout, and
  it is not in `DATA/dataset.index.md`. Registering or migrating it is
  the user's call, so it is left as it is.

Results go in `docs/TFBAF_RESULTS.md`, as a new v0.2 section.

## 8. Open questions for review

- **Q1. Unit selection.** Is the order malignant only → `minLib` →
  primary-derived → first ID right? Or should the two library screens be
  averaged when both exist? We prefer selection: LFC scales differ
  between libraries.
- **Q2. Screen-depth adjustment.** Should depth adjustment be the
  sensitivity, as proposed, or the primary organoid statistic? Raw is
  primary for consistency with v0.1. Depth is the more likely artefact in
  organoids, though.
- **Q3. The 25% prevalence cut-off for "dependent".** Is it reasonable
  for both 2D (−0.5) and organoids (`is_depleted`)?
- **Q4. Background percentile in TCGA.** Is median TPM ≥ 1 the right
  gene universe, or should it be restricted to protein-coding genes?

## 9. Review history

- **Round 8** (head of PR #6): **APPROVE**, 0 blocking.
  - N1: report the included library counts and the LFC by library → §3.1.
  - N2: report raw and adjusted ρ and the change; depth-sensitive is not
    robust → §3.3.
  - N3: the pattern is "no positive co-dependency"; negative is flagged
    separately → §4.
  - N4: signed and absolute percentiles → §5.
  - N5: exclude the TF itself; report the background n, missing values
    and zero variance; prostate is exploratory → §5.
  - Q1: selection, not averaging. Q2: raw is primary, depth-adjusted is
    the sensitivity. Q3: 25% accepted (descriptive; not comparable across
    2D and organoids).
  - Q4: protein-coding with median TPM ≥ 1 → §5 (HGNC + UniProt proxy).
