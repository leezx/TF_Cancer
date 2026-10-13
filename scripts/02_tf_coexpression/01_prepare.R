#!/usr/bin/env Rscript
# Step 2 design §3.1–§3.4: per cohort, vst (tumors; normals separately for the
# ≥40-normal cohorts), the analysed TF universe, the proliferation score, MKI67,
# and ABSOLUTE purity. Writes prep/<cohort>.rds under the result root.
# Usage (repo root): Rscript scripts/02_tf_coexpression/01_prepare.R [cohort ...]
source("scripts/02_tf_coexpression/common2.R")
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(TRUE); todo <- if (length(a)) a else cohorts

tf <- fread("results/01_tcga_tf_de/TF_reference.tsv")[in_expression_matrix == TRUE]
sig <- fread(file.path(repo_out, "proliferation_signature_genes.tsv"))[status == "in_signature"]
mki <- fread(file.path(repo_out, "proliferation_signature_genes.tsv"))[symbol == "MKI67", gencode_v23_id]
stopifnot(length(mki) == 1, !any(sub("\\..*$", "", sig$gencode_v23_id) %in% tf$ensembl_id))
pur <- fread(file.path(root, "1.Databases", "TCGA_PanCanAtlas_ABSOLUTE_purity", "raw", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"))
pur <- pur[`call status` == "called", .(sample = array, purity)]
st <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(root, "1.Databases", "TCGA_PanCanAtlas_subtypes", "raw", "TCGASubtype.20170308.tsv.gz"))))
pam <- st[, .(sample = sampleID, pam50 = Subtype_mRNA)]

do_vst <- function(cnt) {
  dds <- DESeqDataSetFromMatrix(cnt, data.frame(row.names = colnames(cnt), x = rep(1, ncol(cnt))), ~1)
  assay(vst(dds, blind = TRUE, nsub = 1000, fitType = "parametric"))
}
score_of <- function(V) {  # V: genes x samples (vst); mean of per-gene z over signature genes present
  g <- intersect(sig$gencode_v23_id, rownames(V))
  Z <- t(scale(t(V[g, , drop = FALSE]))); Z <- Z[is.finite(rowSums(Z)), , drop = FALSE]
  list(score = colMeans(Z), n_genes = nrow(Z))
}

for (k in todo) {
  x <- readRDS(file.path(processed_dir, sprintf("toil_%s.rds", k)))
  s <- x$samples
  tum <- s[group == "tumor", sample]
  # §3.3 TF universe: tumor median TPM >= 1 and prevalence (TPM >= 1) >= 25%
  tp <- x$tpm_tf[tf$gencode_v23_id, tum, drop = FALSE]
  keep <- apply(tp, 1, median) >= P$tpm_min & rowMeans(tp >= P$tpm_min) >= P$prev_min
  tf_ids <- tf$gencode_v23_id[keep]
  stopifnot(all(tf_ids %in% rownames(x$counts)))
  out <- list(cohort = k, tf_ids = tf_ids, tf_symbols = tf$hgnc_symbol[keep])
  for (grp in c("tumor", "normal")) {
    ids <- s[group == grp, sample]
    if (grp == "normal" && length(ids) < P$normal_n_min) next
    V <- do_vst(x$counts[, ids, drop = FALSE])
    sc <- score_of(V)
    smp <- data.table(sample = ids, prolif = sc$score, MKI67 = V[mki, ids])
    smp <- merge(smp, pur, by = "sample", all.x = TRUE, sort = FALSE)
    if (k == "BRCA") smp <- merge(smp, pam, by = "sample", all.x = TRUE, sort = FALSE)
    smp <- smp[match(ids, sample)]
    X <- t(V[tf_ids, ids, drop = FALSE]); colnames(X) <- out$tf_symbols
    out[[grp]] <- list(X = X, samples = smp, n_signature_genes = sc$n_genes)
    rm(V); gc(verbose = FALSE)
  }
  out$n_tumor <- length(tum)
  out$n_purity <- sum(!is.na(out$tumor$samples$purity))
  out$primary_eligible <- out$n_purity / out$n_tumor >= P$purity_cov_min && out$n_purity >= P$purity_n_min
  saveRDS(out, file.path(res_dir, "prep", sprintf("%s.rds", k)))
  log_msg("%s: %d TFs analysed; tumors %d (purity %d, eligible %s); normals %s; signature genes %d",
          k, length(tf_ids), out$n_tumor, out$n_purity, out$primary_eligible,
          if (is.null(out$normal)) "-" else nrow(out$normal$X), out$tumor$n_signature_genes)
}
writeLines(capture.output(sessionInfo()), file.path(res_dir, "report", "sessionInfo_01_prepare.txt"))
writeLines(c(sprintf("vst: blind=TRUE, nsub=1000, fitType=parametric, no prior filter; DESeq2 %s", packageVersion("DESeq2"))),
           file.path(res_dir, "report", "vst_parameters.txt"))
