# Shared settings and helpers for Step 2 (docs/STEP2_TF_COEXPRESSION_DESIGN.md).
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
toil_dir <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx")
processed_dir <- file.path(toil_dir, "processed", "v0.1")
analysis_id <- "2026-10-13_tf_cancer_step2_coexpr_v0.1"
res_dir <- file.path(toil_dir, "result", analysis_id)
repo_out <- "results/02_tf_coexpression"
for (d in c(file.path(res_dir, c("prep", "coexpr", "report")), file.path(repo_out, c("qc", "figures"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Preset parameters (design §3; never changed after seeing results)
P <- list(tpm_min = 1, prev_min = 0.25, coupled_rho = 0.5, edge_rho = 0.4, fdr = 0.05,
          recur_cohort_frac = 0.5, recur_min_cohorts = 3, recur_organ_systems = 3,
          purity_cov_min = 0.70, purity_n_min = 50, normal_n_min = 40,
          n_perm = 20, n_boot = 1000, seed = 1, leiden_res = 1, leiden_iter = 10, module_min = 5,
          tumor_enriched_normal_rho_max = 0.2, mki67_check = 0.7, luma_n_min = 100)
cohorts <- c("COADREAD", "LUAD", "LUSC", "KIRC", "KIRP", "KICH", "LIHC", "LAML", "PAAD", "BRCA",
             "GBM", "LGG", "PRAD", "THCA", "STAD", "HNSC")
controls <- data.table(cohort = c("BRCA", "COADREAD", "LIHC", "LUSC"),
                       a = c("FOXA1", "CDX1", "HNF1A", "TP63"), b = c("GATA3", "CDX2", "HNF4A", "SOX2"))

log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", sprintf(...))
fwrite_tsv <- function(x, f) fwrite(x, f, sep = "\t", na = "NA")

# Rank-based partial correlation matrix (design §3.4/§3.5): rank each column of X
# (samples x genes), residualize on the ranks of covariates Z (intercept added),
# Pearson correlation of residuals. Returns list(r, df).
resid_ranks <- function(X, Z = NULL) {
  R <- apply(X, 2, rank)
  if (is.null(Z) || ncol(as.matrix(Z)) == 0) return(scale(R, scale = FALSE))
  Zr <- cbind(1, apply(as.matrix(Z), 2, rank))
  qr.resid(qr(Zr), R)
}
pcor_mat <- function(X, Z = NULL) {
  k <- if (is.null(Z)) 0 else ncol(as.matrix(Z))
  E <- resid_ranks(X, Z)
  list(r = suppressWarnings(cor(E)), df = nrow(X) - 2 - k, E = E)
}
r2p <- function(r, df) 2 * pt(-abs(r * sqrt(df / pmax(1 - r^2, 1e-12))), df)
