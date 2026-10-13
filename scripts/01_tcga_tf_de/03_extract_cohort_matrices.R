#!/usr/bin/env Rscript
# Design §3.4 step 0 (reconstruction QC) + §4: read the Toil expected_count and
# TPM matrices once, write one RDS per cohort to processed/v0.1:
#   counts  integer matrix, all genes (Toil versioned IDs) x cohort samples
#   tpm_tf  linear TPM, TFs x cohort samples
#   samples sample-manifest rows
# Memory-safe for a 16 GB machine (computingResource.md): counts are read for
# the manifest columns only and converted to integer column by column; TPM is
# row-filtered to TFs while streaming.
# Usage (repo root): Rscript scripts/01_tcga_tf_de/03_extract_cohort_matrices.R
source("scripts/01_tcga_tf_de/common.R")
check_md5(c("counts", "tpm"))
setDTthreads(as.integer(Sys.getenv("NSLOTS", "4")))

sm <- fread(file.path(repo_results, "sample_manifest.tsv"))
tf <- fread(file.path(repo_results, "TF_reference.tsv"))[in_expression_matrix == TRUE]
ids <- unique(sm$sample)
gunzip <- function(f) sprintf("gzip -dc %s", shQuote(f))
log_msg("Reading counts for %d samples", length(ids))

cnt <- fread(cmd = gunzip(in_path("counts")), select = c("sample", ids), showProgress = FALSE)
genes <- cnt$sample; cnt[, sample := NULL]
stopifnot(ncol(cnt) == length(ids), !anyDuplicated(genes))

# Reconstruction QC by expression bin; convert each column to integer in place.
labs <- c("0", "1-10", "10-100", "100-1e3", ">1e3")
brk <- c(-Inf, 0.5, 10, 100, 1000, Inf)
n <- sum_abs <- max_abs <- max_rel <- setNames(numeric(5), labs)
min_log <- Inf
for (j in names(cnt)) {
  x <- cnt[[j]]; min_log <- min(min_log, x)
  y <- 2^x - 1; r <- pmax(round(y), 0); e <- abs(y - r)
  b <- as.integer(cut(r, brk))
  n <- n + tabulate(b, 5)
  sum_abs <- sum_abs + vapply(1:5, function(i) sum(e[b == i]), 0)
  max_abs <- pmax(max_abs, vapply(1:5, function(i) max(c(0, e[b == i])), 0))
  max_rel <- pmax(max_rel, vapply(1:5, function(i) max(c(0, (e / pmax(r, 1))[b == i])), 0))
  set(cnt, j = j, value = as.integer(r))
}
stopifnot(min_log >= 0)
qc_dir <- file.path(repo_results, "qc", "input_validation")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
acc <- data.table(bin = labs, n = n, mean_abs = sum_abs / pmax(n, 1), max_abs = max_abs, max_rel = max_rel)
fwrite_tsv(acc, file.path(qc_dir, "toil_reconstruction_error_by_bin.tsv"))
print(acc)

libdt <- data.table(sample = names(cnt), library_size = vapply(cnt, function(v) sum(as.numeric(v)), 0),
                    zero_fraction = vapply(cnt, function(v) mean(v == 0L), 0))
libdt <- merge(sm[, .(sample, cohort, group)], libdt, by = "sample", allow.cartesian = TRUE)
libdt[, outlier_3mad := abs(library_size - median(library_size)) > 3 * mad(library_size), by = .(cohort, group)]
fwrite_tsv(libdt, file.path(qc_dir, "toil_library_sizes.tsv"))
fwrite_tsv(libdt[, .(n = .N, median_library = median(library_size), median_zero_fraction = median(zero_fraction),
                     n_outlier_3mad = sum(outlier_3mad)), by = .(cohort, group)],
           file.path(qc_dir, "toil_library_summary.tsv"))

# TPM: TF rows only (header + TF IDs), filtered while streaming.
idf <- tempfile(); writeLines(tf$gencode_v23_id, idf)
tpm <- fread(cmd = sprintf("%s | awk -F'\\t' 'NR==FNR{k[$1];next} FNR==1 || ($1 in k)' %s -", gunzip(in_path("tpm")), shQuote(idf)),
             select = c("sample", ids), showProgress = FALSE)
stopifnot(nrow(tpm) == nrow(tf))
tpm_tf <- as.matrix(tpm[, -1]); rownames(tpm_tf) <- tpm$sample
tpm_tf <- pmax(2^tpm_tf[tf$gencode_v23_id, , drop = FALSE] - 0.001, 0)
rm(tpm); gc(verbose = FALSE)

for (k in unique(sm$cohort)) {
  s <- sm[cohort == k]
  m <- as.matrix(cnt[, s$sample, with = FALSE]); rownames(m) <- genes
  saveRDS(list(counts = m, tpm_tf = tpm_tf[, s$sample, drop = FALSE], samples = s),
          file.path(processed_dir, sprintf("toil_%s.rds", k)))
  log_msg("wrote %s: %d samples", k, nrow(s))
  rm(m); gc(verbose = FALSE)
}
write_run_record(file.path(processed_dir, "run_records"), "03_extract_cohort_matrices")
