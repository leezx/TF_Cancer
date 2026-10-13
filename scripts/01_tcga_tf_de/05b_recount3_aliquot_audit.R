#!/usr/bin/env Rscript
# Review round 3 N1: write the recount3 retained/excluded aliquot table for the
# validation cohorts without rerunning DESeq2 (same loader as 04).
# Usage (repo root): Rscript scripts/01_tcga_tf_de/05b_recount3_aliquot_audit.R
source("scripts/01_tcga_tf_de/common.R")
source("scripts/01_tcga_tf_de/recount3_loader.R")
for (k in names(cfg$inputs$recount3)) invisible(load_recount3(k))
