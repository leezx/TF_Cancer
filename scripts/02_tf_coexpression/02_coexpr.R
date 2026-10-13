#!/usr/bin/env Rscript
# Step 2 design §3.4–§3.6, §4: per cohort,
#   TF–proliferation association (purity-adjusted for primary-eligible cohorts; unadjusted for all),
#   TF–TF co-expression for three analyses (primary = purity + prolif on purity-covered tumors of
#   eligible cohorts; prolif_only and unadjusted on all tumors), edges (|rho| >= 0.4, BH q < 0.05),
#   20-permutation null (diagnostic), normals (prolif_only, cohorts with >= 40 normals), QC.
# Usage (repo root): Rscript scripts/02_tf_coexpression/02_coexpr.R [cohort ...]
source("scripts/02_tf_coexpression/common2.R")
a <- commandArgs(TRUE); todo <- if (length(a)) a else cohorts

edges_of <- function(r, df, analysis, k) {
  ut <- which(upper.tri(r), arr.ind = TRUE)
  rho <- r[ut]; ok <- is.finite(rho)
  p <- r2p(rho, df); q <- rep(NA_real_, length(p)); q[ok] <- p.adjust(p[ok], "BH")
  sel <- which(ok & abs(rho) >= P$edge_rho & q < P$fdr)
  list(edges = data.table(cohort = k, analysis = analysis, TF_a = rownames(r)[ut[sel, 1]], TF_b = colnames(r)[ut[sel, 2]],
                          rho = rho[sel], q = q[sel], sign = fifelse(rho[sel] > 0, "positive", "negative")),
       n_pairs = sum(ok))
}
perm_null <- function(E, df) {
  set.seed(P$seed)
  vapply(seq_len(P$n_perm), function(i) {
    Ep <- apply(E, 2, sample)
    nrow(edges_of(suppressWarnings(cor(Ep)), df, "perm", "")$edges)
  }, 0)
}

for (k in todo) {
  pr <- readRDS(file.path(res_dir, "prep", sprintf("%s.rds", k)))
  X <- pr$tumor$X; sm <- pr$tumor$samples
  stopifnot(identical(rownames(X), sm$sample) || is.null(rownames(X)))
  has_p <- !is.na(sm$purity)
  an <- list()
  if (pr$primary_eligible) an$primary <- list(X = X[has_p, ], Z = as.matrix(sm[has_p, .(purity, prolif)]))
  an$prolif_only <- list(X = X, Z = as.matrix(sm[, .(prolif)]))
  an$unadjusted <- list(X = X, Z = NULL)

  # ---- TF–TF co-expression ----
  mats <- list(); edges <- list(); qc <- list()
  for (nm in names(an)) {
    pc <- pcor_mat(an[[nm]]$X, an[[nm]]$Z)
    e <- edges_of(pc$r, pc$df, nm, k)
    perm <- perm_null(pc$E, pc$df)
    mats[[nm]] <- pc$r; edges[[nm]] <- e$edges
    qc[[nm]] <- data.table(cohort = k, analysis = nm, n_samples = nrow(an[[nm]]$X), n_TF = ncol(X), n_pairs = e$n_pairs,
                           n_edges_pos = sum(e$edges$sign == "positive"), n_edges_neg = sum(e$edges$sign == "negative"),
                           perm_edges_mean = mean(perm), perm_edges_max = max(perm))
    log_msg("%s %s: n=%d, edges +%d/-%d, permuted mean %.1f max %d", k, nm, nrow(an[[nm]]$X),
            qc[[nm]]$n_edges_pos, qc[[nm]]$n_edges_neg, mean(perm), max(perm))
  }
  # normals: prolif_only, same TFs (design §3.6)
  if (!is.null(pr$normal)) {
    Xn <- pr$normal$X
    mats$normal_prolif_only <- pcor_mat(Xn, as.matrix(pr$normal$samples[, .(prolif)]))$r
  }

  # ---- TF–proliferation association (§3.4) ----
  pa <- list(data.table(cohort = k, analysis = "unadjusted", TF = colnames(X), n = nrow(X),
                        rho = suppressWarnings(cor(apply(X, 2, rank), rank(sm$prolif)))[, 1]))
  if (pr$primary_eligible) {
    Zp <- cbind(1, rank(sm$purity[has_p]))
    ex <- qr.resid(qr(Zp), apply(X[has_p, ], 2, rank)); es <- qr.resid(qr(Zp), rank(sm$prolif[has_p]))
    pa[[2]] <- data.table(cohort = k, analysis = "purity_adjusted", TF = colnames(X), n = sum(has_p),
                          rho = suppressWarnings(cor(ex, es))[, 1])
  }
  pa <- rbindlist(pa)
  pa[, df := n - 2 - fifelse(analysis == "purity_adjusted", 1L, 0L)]
  pa[, p := r2p(rho, df)][, q := p.adjust(p, "BH"), by = analysis]
  pa[, coupled := rho >= P$coupled_rho & q < P$fdr]

  # ---- cohort QC (§4) ----
  cq <- data.table(cohort = k, n_tumor = pr$n_tumor, n_purity = pr$n_purity, purity_coverage = pr$n_purity / pr$n_tumor,
                   primary_eligible = pr$primary_eligible, n_TF = ncol(X), n_signature_genes = pr$tumor$n_signature_genes,
                   rho_purity_prolif = if (any(has_p)) cor(sm$purity[has_p], sm$prolif[has_p], method = "spearman") else NA_real_,
                   rho_score_MKI67 = cor(sm$prolif, sm$MKI67, method = "spearman"),
                   n_normal = if (is.null(pr$normal)) 0L else nrow(pr$normal$X),
                   rho_normal_score_MKI67 = if (is.null(pr$normal)) NA_real_ else cor(pr$normal$samples$prolif, pr$normal$samples$MKI67, method = "spearman"))
  adj <- if (all(c("E2F1", "FOXM1") %in% colnames(X)))
    data.table(cohort = k, analysis = names(mats), rho_E2F1_FOXM1 = vapply(mats, function(m) m["E2F1", "FOXM1"], 0)) else NULL

  saveRDS(list(mats = mats, tf = colnames(X), primary_eligible = pr$primary_eligible),
          file.path(res_dir, "coexpr", sprintf("%s.rds", k)))
  fwrite_tsv(rbindlist(edges), file.path(res_dir, "coexpr", sprintf("%s_edges.tsv.gz", k)))
  fwrite_tsv(pa, file.path(res_dir, "coexpr", sprintf("%s_prolif_assoc.tsv", k)))
  fwrite_tsv(rbindlist(qc), file.path(res_dir, "coexpr", sprintf("%s_analysis_qc.tsv", k)))
  fwrite_tsv(cq, file.path(res_dir, "coexpr", sprintf("%s_cohort_qc.tsv", k)))
  if (!is.null(adj)) fwrite_tsv(adj, file.path(res_dir, "coexpr", sprintf("%s_adjustment_check.tsv", k)))
  log_msg("%s: rho(score, MKI67) = %.2f; rho(purity, prolif) = %.2f; coupled TFs (unadj/purity-adj) = %d/%s",
          k, cq$rho_score_MKI67, cq$rho_purity_prolif, pa[analysis == "unadjusted", sum(coupled)],
          if (pr$primary_eligible) pa[analysis == "purity_adjusted", sum(coupled)] else "-")
}
