# recount3_TCGA_LUAD

- Canonical DATA path: `DATA/bulkRNAseq/recount3_TCGA_LUAD/` (local and Argos; 55-file md5 match with the COADREAD recount3 set on 2026-10-13)
- Source: recount3 (duffel.rail.bio/recount3), R recount3 1.20.0, gencode_v26, `transform_counts`; public
- Files: `raw/recount3_cache/` (gene_sums, metadata, G026 GTF; `raw/checksums.md5`), `processed/v0.1/recount3_TCGA_LUAD_gene_counts.rds` (63,856 genes × 601 samples) + coldata + `checksums.md5`
- Why: Step 1 independent-pipeline validation gate, 2nd cohort (design §3.4 step 0; review round 2). Acquired 2026-10-13 with user approval
- Also used (pre-existing, now mirrored to Argos): `DATA/bulkRNAseq/recount3_TCGA_COAD_READ_GTEx_COLON/`
