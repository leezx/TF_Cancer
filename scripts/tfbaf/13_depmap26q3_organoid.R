#!/usr/bin/env Rscript
# TF–BAF dependency correlation in DepMap Public 26Q3 NextGen organoid screens (ScreenType 3DO).
# Units: one value per model (screens of the same model averaged, as DepMap does for organoids);
# malignant models only (adenoma / non-cancerous excluded). CRC = Bowel, Breast, Prostate; AML and
# melanoma have no organoid screens. Spearman rho, BH within pair (min n 10), and background z:
# (rho - mean) / SD of the TF's rho with every gene in the same models.
# Usage (repo root): Rscript scripts/tfbaf/13_depmap26q3_organoid.R
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
raw <- file.path(root, "1.Databases", "DepMap_26Q3", "raw")
out <- "results/tfbaf"
baf <- fread(file.path(out, "BAF_members.tsv"))
lin <- data.table(pair = c("CRC", "Breast", "Prostate"), lineage = c("Bowel", "Breast", "Prostate"), TF = c("SOX9", "GATA3", "FOXA1"))

s <- fread(file.path(raw, "ScreenSequenceMap.csv"))
m <- fread(file.path(raw, "Model.csv"))
G <- fread(file.path(raw, "ScreenGeneEffect.csv"))
setnames(G, 1, "ScreenID"); setnames(G, sub(" \\(.*$", "", names(G)))
G <- G[, !duplicated(names(G)), with = FALSE]
u <- unique(s[ScreenType == "3DO", .(ScreenID, ModelID)])[ScreenID %in% G$ScreenID]
u <- merge(u, m[, .(ModelID, OncotreeLineage, OncotreePrimaryDisease, ModelType)], by = "ModelID")
u <- u[!OncotreePrimaryDisease %in% c("Non-Cancerous", "Tubular Adenoma of the Colon")]
u <- merge(u, lin, by.x = "OncotreeLineage", by.y = "lineage")
fwrite(u[order(pair, ModelID)], file.path(out, "qc", "depmap26q3_organoid_units.tsv"), sep = "\t")

bh <- function(p) p.adjust(p, "BH")
res <- list(); bgs <- list()
for (k in seq_len(nrow(lin))) {
  pp <- lin$pair[k]; tf <- lin$TF[k]; uu <- u[pair == pp]
  X <- G[ScreenID %in% uu$ScreenID]
  X[, ModelID := uu$ModelID[match(ScreenID, uu$ScreenID)]]
  X <- X[, lapply(.SD, mean, na.rm = TRUE), by = ModelID, .SDcols = setdiff(names(X), c("ScreenID", "ModelID"))]
  M <- as.matrix(X[, -1]); n_mod <- nrow(M)
  tfv <- M[, tf]
  Gm <- M[, colnames(M) != tf, drop = FALSE]
  ok <- colSums(is.na(Gm)) == 0 & apply(Gm, 2, function(v) var(v) > 0)
  rbg <- as.vector(cor(apply(Gm[, ok], 2, rank), rank(tfv)))
  bgs[[pp]] <- data.table(source = "DepMap26Q3_3DO", pair = pp, bg_mean = mean(rbg), bg_sd = sd(rbg), bg_n = sum(ok), n_dropped_na_or_zero_var = sum(!ok))
  r <- rbindlist(lapply(baf$symbol, function(g) {
    v <- M[, g]; okk <- !is.na(v) & !is.na(tfv); nn <- sum(okk)
    ct <- if (nn >= 10) suppressWarnings(cor.test(tfv[okk], v[okk], method = "spearman", exact = FALSE)) else NULL
    data.table(BAF = g, n = nn, rho = cor(rank(tfv[okk]), rank(v[okk])), p = if (is.null(ct)) NA_real_ else ct$p.value,
               BAF_median_effect = median(v, na.rm = TRUE))
  }))
  r[, `:=`(pair = pp, TF = tf, q = if (n_mod >= 10) bh(p) else NA_real_, underpowered_n_lt_10 = n_mod < 10,
           z = (rho - mean(rbg)) / sd(rbg), TF_median_effect = median(tfv), TF_frac_lt_m0.5 = mean(tfv < -0.5))]
  res[[pp]] <- r
}
res <- rbindlist(res); bgs <- rbindlist(bgs)
setcolorder(res, c("pair", "TF", "BAF", "n", "rho", "p", "q", "z"))
fwrite(res, file.path(out, "DepMap26Q3_organoid_TF_BAF_correlation.tsv"), sep = "\t", na = "NA")
fwrite(bgs, file.path(out, "qc", "depmap26q3_organoid_background_mean_sd.tsv"), sep = "\t")
print(u[, .(models = uniqueN(ModelID), screens = .N), by = .(pair, OncotreePrimaryDisease)])
print(bgs)
print(res[order(pair, -z)][, head(.SD, 5), by = pair][, .(pair, BAF, n, rho = round(rho, 2), q = signif(q, 2), z = round(z, 1))])
print(res[!is.na(q) & q < 0.05, .(pair, BAF, rho, q)])
