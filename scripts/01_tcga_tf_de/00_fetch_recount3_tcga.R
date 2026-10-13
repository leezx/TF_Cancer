#!/usr/bin/env Rscript
# Fetch recount3 (Monorail) TCGA gene counts for one project, for the Step 1
# independent-pipeline validation gate (design §3.4 step 0). Same procedure
# as DATA/bulkRNAseq/recount3_TCGA_COAD_READ_GTEx_COLON (gencode_v26,
# transform_counts coverage -> read counts).
#
# Usage: Rscript 00_fetch_recount3_tcga.R <dataset_dir> <TCGA project, e.g. LUAD>
#   <dataset_dir>/raw/recount3_cache   BiocFileCache (immutable source files)
#   <dataset_dir>/processed/v0.1       counts RSE + coldata + sessionInfo

suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript 00_fetch_recount3_tcga.R <dataset_dir> <project>")
dataset_dir <- args[[1]]
proj <- args[[2]]
cache_dir <- file.path(dataset_dir, "raw", "recount3_cache")
out_dir <- file.path(dataset_dir, "processed", "v0.1")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(recount3_url = "https://duffel.rail.bio/recount3")
bfc <- BiocFileCache::BiocFileCache(cache_dir, ask = FALSE)

projects <- available_projects(bfc = bfc)
wanted <- subset(projects, file_source == "tcga" & project == proj)
stopifnot(nrow(wanted) == 1)

rse <- create_rse(wanted, type = "gene", annotation = "gencode_v26", bfc = bfc)
assay(rse, "counts") <- transform_counts(rse)  # coverage -> read counts
assays(rse) <- assays(rse)["counts"]

write.csv(as.data.frame(colData(rse)), file.path(out_dir, sprintf("coldata_tcga_%s.csv", proj)), row.names = TRUE)
saveRDS(rse, file.path(out_dir, sprintf("recount3_TCGA_%s_gene_counts.rds", proj)))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_00_fetch_recount3.txt"))
cat(sprintf("Saved %d genes x %d samples\n", nrow(rse), ncol(rse)))
print(table(rse$tcga.cgc_sample_sample_type, useNA = "ifany"))
