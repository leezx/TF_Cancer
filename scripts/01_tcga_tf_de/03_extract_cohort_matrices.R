#!/usr/bin/env Rscript
# Design §3.4 step 0 (reconstruction QC) + §4: read the Toil expected_count and
# TPM matrices once, write one RDS per cohort to processed/v0.1:
#   counts  integer matrix, all genes (Toil versioned IDs) x cohort samples
#   tpm_tf  linear TPM, TFs x cohort samples
#   samples sample-manifest rows
# Usage (repo root, qsub): Rscript scripts/01_tcga_tf_de/03_extract_cohort_matrices.R
source("scripts/01_tcga_tf_de/common.R")
check_md5(c("counts", "tpm"))
setDTthreads(as.integer(Sys.getenv("NSLOTS", "4")))

sm <- fread(file.path(repo_results, "sample_manifest.tsv"))
tf <- fread(file.path(repo_results, "TF_reference.tsv"))[in_expression_matrix == TRUE]
ids <- unique(sm$sample)
log_msg("Reading %d samples", length(ids))

read_cols <- function(f) {
  x <- fread(cmd = sprintf("zcat %s", shQuote(f)), select = c("sample", ids), showProgress = FALSE)
  m <- as.matrix(x[, -1]); rownames(m) <- x$sample; m
}

cnt_log <- read_cols(in_path("counts"))
log_msg("counts: %d x %d; min log value %.4f", nrow(cnt_log), ncol(cnt_log), min(cnt_log))
stopifnot(min(cnt_log) >= 0)

# Reconstruction QC by expression bin, in column chunks to bound memory.
bins <- c(-Inf, 0.5, 10, 100, 1000, Inf)
labs <- c("0", "1-10", "10-100", "100-1e3", ">1e3")
acc <- data.table(bin = labs, n = 0, sum_abs = 0, max_abs = 0, max_rel = 0)
lib <- vector("list", ncol(cnt_log))
cnt <- matrix(0L, nrow(cnt_log), ncol(cnt_log), dimnames = dimnames(cnt_log))
for (s in split(seq_len(ncol(cnt_log)), ceiling(seq_len(ncol(cnt_log)) / 500))) {
  y <- 2^cnt_log[, s, drop = FALSE] - 1
  r <- round(y); r[r < 0] <- 0
  e <- abs(y - r)
  b <- cut(r, bins, labels = labs)
  agg <- data.table(bin = as.character(b), e = as.vector(e), rel = as.vector(e / pmax(r, 1)))[
    , .(n = .N, sum_abs = sum(e), max_abs = max(e), max_rel = max(rel)), by = bin]
  acc <- merge(acc, agg, by = "bin", all.x = TRUE, suffixes = c("", ".new"))
  for (v in c("n", "sum_abs")) acc[!is.na(get(paste0(v, ".new"))), (v) := get(v) + get(paste0(v, ".new"))]
  for (v in c("max_abs", "max_rel")) acc[!is.na(get(paste0(v, ".new"))), (v) := pmax(get(v), get(paste0(v, ".new")))]
  acc[, grep("\\.new$", names(acc), value = TRUE) := NULL]
  storage.mode(r) <- "integer"
  cnt[, s] <- r
  rm(y, r, e, b); gc(verbose = FALSE)
}
acc[, mean_abs := sum_abs / n]
acc <- acc[match(labs, bin), .(bin, n, mean_abs, max_abs, max_rel)]
qc_dir <- file.path(repo_results, "qc", "input_validation")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
fwrite_tsv(acc, file.path(qc_dir, "toil_reconstruction_error_by_bin.tsv"))
print(acc)
rm(cnt_log); gc(verbose = FALSE)

libdt <- data.table(sample = colnames(cnt), library_size = colSums(cnt), zero_fraction = colMeans(cnt == 0))
libdt <- merge(sm[, .(sample, cohort, group)], libdt, by = "sample", allow.cartesian = TRUE)
libdt[, med := median(library_size), by = .(cohort, group)]
libdt[, mad := mad(library_size), by = .(cohort, group)]
libdt[, outlier_3mad := abs(library_size - med) > 3 * mad]
fwrite_tsv(libdt[, .(sample, cohort, group, library_size, zero_fraction, outlier_3mad)],
           file.path(qc_dir, "toil_library_sizes.tsv"))
fwrite_tsv(libdt[, .(n = .N, median_library = median(library_size), median_zero_fraction = median(zero_fraction),
                     n_outlier_3mad = sum(outlier_3mad)), by = .(cohort, group)],
           file.path(qc_dir, "toil_library_summary.tsv"))

tpm_log <- read_cols(in_path("tpm"))
tf_rows <- match(tf$gencode_v23_id, rownames(tpm_log))
stopifnot(!anyNA(tf_rows), identical(rownames(tpm_log), rownames(cnt)))
tpm_tf <- pmax(2^tpm_log[tf_rows, , drop = FALSE] - 0.001, 0)
rm(tpm_log); gc(verbose = FALSE)

for (k in unique(sm$cohort)) {
  s <- sm[cohort == k]
  saveRDS(list(counts = cnt[, s$sample, drop = FALSE], tpm_tf = tpm_tf[, s$sample, drop = FALSE], samples = s),
          file.path(processed_dir, sprintf("toil_%s.rds", k)))
  log_msg("wrote %s: %d samples", k, nrow(s))
}
write_run_record(file.path(processed_dir, "run_records"), "03_extract_cohort_matrices")
