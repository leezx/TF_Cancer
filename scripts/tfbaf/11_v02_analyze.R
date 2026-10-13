#!/usr/bin/env Rscript
# TF–BAF v0.2 (docs/TFBAF_V02_DESIGN.md): §3 organoid co-dependency (CRC SOX9) with
# depth / library / primary-only sensitivities and positive control; §4 DepMap 2D
# context (SD, BAF prevalence, bootstrap CI, dependency pattern); §5 background
# percentiles (TCGA, DepMap, organoid); §6 v0.2 figures. v0.1 outputs are not touched.
# Usage (repo root): Rscript scripts/tfbaf/11_v02_analyze.R   (after 10_v02_extract.R)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
res_root <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx", "result", "2026-10-13_tf_cancer_tfbaf_v0.2")
org <- file.path(res_root, "organoid")
for (d in file.path(res_root, c("depmap_context", "background", "report"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
out <- "results/tfbaf"; qc <- file.path(out, "qc"); fig <- file.path(out, "figures")
baf <- fread(file.path(out, "BAF_members.tsv"))
cls_order <- c("ATPase", "core", "cBAF-specific", "PBAF-specific", "ncBAF-specific")
pair_order <- c("CRC", "AML", "Prostate", "Breast", "Melanoma")
tf_of <- setNames(c("SOX9", "SPI1", "FOXA1", "GATA3", "SOX10"), pair_order)
N_BOOT <- 1000; SEED <- 1; DEP_PREV <- 0.25

call_fun <- function(rho, q) fcase(is.na(rho) | is.na(q) | q >= 0.05, "none",
                                   rho >= 0.3, "strong", rho >= 0.2, "moderate",
                                   rho <= -0.3, "strong_negative", rho <= -0.2, "moderate_negative", default = "none")
spear <- function(a, b) {
  ok <- complete.cases(a, b)
  if (sum(ok) < 10) return(list(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman", exact = FALSE))
  list(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
pspear <- function(a, b, Z) {  # as v0.1 02_analyze.R
  ok <- complete.cases(a, b, Z); n <- sum(ok)
  if (n < 10) return(list(n = n, rho = NA_real_, p = NA_real_))
  Zr <- apply(as.matrix(Z)[ok, , drop = FALSE], 2, rank)
  ra <- resid(lm(rank(a[ok]) ~ Zr)); rb <- resid(lm(rank(b[ok]) ~ Zr))
  r <- cor(ra, rb); df <- n - 2 - ncol(Zr)
  list(n = n, rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df))
}
# bootstrap percentile CI of Spearman rho: same resampled rows for all members of a pair
boot_ci <- function(tfv, B) {
  ok_any <- !is.na(tfv); n <- length(tfv)
  set.seed(SEED); idx <- replicate(N_BOOT, sample.int(n, n, replace = TRUE))
  t(vapply(colnames(B), function(g) {
    r <- apply(idx, 2, function(i) suppressWarnings(cor(tfv[i], B[i, g], method = "spearman", use = "complete.obs")))
    quantile(r, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  }, numeric(2)))
}
pattern_of <- function(tf_dep, baf_dep, call) {
  pos <- call %in% c("strong", "moderate")
  fcase(tf_dep & baf_dep & pos, "both dependent, correlated",
        tf_dep & baf_dep, "both dependent, no positive co-dependency",
        tf_dep, "TF-specific dependency", baf_dep, "BAF-specific dependency", default = "neither dependent")
}
pctl <- function(v, bgv) 100 * mean(bgv <= v)

# ==== §3 organoids (CRC SOX9) ====
u <- fread(file.path(org, "crc_selected_units.tsv"))
M <- readRDS(file.path(org, "crc_units_lfc_matrix.rds"))
stopifnot(identical(rownames(M), u$screen), all(c("SOX9", baf$symbol) %in% colnames(M)))
depth <- fread(file.path(org, "screen_depth.tsv"))
a33 <- fread(file.path(org, "all_screens_tf_baf.tsv"))
u[, depth := depth$depth[match(screen, depth$screen)]]
prev_o <- a33[screen %in% u$screen, .(prev = mean(is_depleted == 1)), by = gene]
o_cor <- function(rows, label, adjust = FALSE) {
  rbindlist(lapply(baf$symbol, function(g) {
    r <- if (adjust) pspear(M[rows, "SOX9"], M[rows, g], u$depth[rows]) else spear(M[rows, "SOX9"], M[rows, g])
    data.table(BAF = g, n = r$n, rho = r$rho, p = r$p)
  }))[, `:=`(analysis = label, q = p.adjust(p, "BH"))][, call := call_fun(rho, q)]
}
all_r <- seq_len(nrow(u))
om <- o_cor(all_r, "main")
ci <- boot_ci(M[, "SOX9"], M[, baf$symbol])
om[, `:=`(ci_low = ci[BAF, 1], ci_high = ci[BAF, 2])]
od <- o_cor(all_r, "depth_adjusted", adjust = TRUE)
om <- merge(om, od[, .(BAF, rho_depth_adjusted = rho, q_depth_adjusted = q, call_depth_adjusted = call)], by = "BAF")
om[, rho_change := rho_depth_adjusted - rho]
om[, depth_flag := fcase(call != "none" & sign(rho_depth_adjusted) != sign(rho), "depth_sensitive_reversed",
                         call != "none" & abs(rho_depth_adjusted) < 0.5 * abs(rho), "depth_sensitive_attenuated_gt50pct", default = "")]
om[, `:=`(pair = "CRC", TF = "SOX9", TF_prev_depleted = prev_o[gene == "SOX9", prev],
          BAF_prev_depleted = prev_o$prev[match(BAF, prev_o$gene)],
          TF_median_LFC = median(M[, "SOX9"]), BAF_median_LFC = apply(M[, BAF, drop = FALSE], 2, median),
          TF_sd_LFC = sd(M[, "SOX9"]), BAF_sd_LFC = apply(M[, BAF, drop = FALSE], 2, sd))]
om[, pattern := pattern_of(TF_prev_depleted >= DEP_PREV, BAF_prev_depleted >= DEP_PREV, call)]
om[, negative_flag := call %in% c("strong_negative", "moderate_negative")]
om <- merge(om, baf[, .(BAF = symbol, class)], by = "BAF")
osens <- rbind(o_cor(which(u$library == "minLib"), "minLib_only"), o_cor(which(u$met == "No"), "primary_derived_only"), od)
# positive control over all 178 screens
w33 <- dcast(a33, screen ~ gene, value.var = "LFC")
octrl <- rbindlist(lapply(list(c("SMARCB1", "SMARCE1"), c("ARID1A", "SMARCB1")), function(g) {
  r <- spear(w33[[g[1]]], w33[[g[2]]]); data.table(gene_a = g[1], gene_b = g[2], n = r$n, rho = r$rho, p = r$p) }))
# library QC (N1)
scr <- fread(file.path(org, "screen_inventory.tsv"))
lib_qc <- merge(a33[screen %in% u$screen], scr[, .(screen, library)], by = "screen")[
  , .(n_units = uniqueN(screen), median_LFC = median(LFC), q25 = quantile(LFC, 0.25), q75 = quantile(LFC, 0.75)), by = .(library, gene)]
inv <- scr[, .(n_screens = .N, n_organoids = uniqueN(sample_ID), n_individuals = uniqueN(individual_ID),
               n_malignant_screens = sum(malignant), n_metastasis_derived_screens = sum(met == "Yes", na.rm = TRUE)), by = primary_tumour_type]
inv <- rbind(inv, data.table(primary_tumour_type = "CRC units selected (§3.1)", n_screens = nrow(u), n_organoids = uniqueN(u$sample_ID),
                             n_individuals = uniqueN(u$individual_ID), n_malignant_screens = nrow(u), n_metastasis_derived_screens = sum(u$met == "Yes")))
o_avail <- data.table(pair = pair_order, organoid = c("analysed", rep("unavailable", 4)),
                      reason = c(sprintf("%d CRC carcinoma individuals", nrow(u)), "no haematological organoids",
                                 "no prostate organoids", "no breast organoids", "no melanoma organoids"))

# ==== §4 DepMap 2D context ====
dd <- file.path(root, "1.Databases", "DepMap3D_202608")
ce <- fread(file.path(dd, "CRISPRGeneEffect_DepMap26Q1.csv"))
setnames(ce, 1, "ModelID"); setnames(ce, sub(" \\(.*$", "", names(ce)))
stopifnot(!anyDuplicated(names(ce)))
mod <- fread(file.path(dd, "Model_DepMap26Q1.csv"))
ce <- merge(mod[, .(ModelID, OncotreeLineage, OncotreePrimaryDisease, ModelType)], ce, by = "ModelID")
ce[, pair := fcase(OncotreeLineage == "Bowel", "CRC", OncotreePrimaryDisease == "Acute Myeloid Leukemia", "AML",
                   OncotreeLineage == "Prostate", "Prostate", OncotreeLineage == "Breast", "Breast",
                   OncotreePrimaryDisease %like% "Melanoma", "Melanoma", default = NA_character_)]
dm <- ce[!is.na(pair) & ModelType == "Cell Line"]
dep1 <- fread(file.path(out, "DepMap_TF_BAF_codependency.tsv"))   # v0.1, frozen
genes_all <- setdiff(names(ce), c("ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "ModelType", "pair"))
ctx <- rbindlist(lapply(pair_order, function(pp) {
  d <- dm[pair == pp]; tfv <- d[[tf_of[[pp]]]]
  B <- as.matrix(d[, baf$symbol, with = FALSE]); ci <- boot_ci(tfv, B)
  data.table(pair = pp, BAF = baf$symbol, TF_sd_effect = sd(tfv, na.rm = TRUE), BAF_sd_effect = apply(B, 2, sd, na.rm = TRUE),
             TF_frac_lt_m0.5 = mean(tfv < -0.5, na.rm = TRUE), BAF_frac_lt_m0.5 = colMeans(B < -0.5, na.rm = TRUE),
             ci_low = ci[, 1], ci_high = ci[, 2])
}))
ctx <- merge(dep1[, .(pair, TF, BAF, class, n, rho, q, call, underpowered, TF_not_dependency)], ctx, by = c("pair", "BAF"))
ctx[, pattern := pattern_of(TF_frac_lt_m0.5 >= DEP_PREV, BAF_frac_lt_m0.5 >= DEP_PREV, call)]
ctx[, negative_flag := call %in% c("strong_negative", "moderate_negative")]

# ==== §5 background percentiles ====
bg_tcga <- fread(file.path(res_root, "background", "tcga_tf_allgene_rho.tsv.gz"))
allrank_cor <- function(tfv, G) {   # Spearman of tfv vs each column of G; exact for columns with NA
  na_col <- colSums(is.na(G)) > 0
  r <- rep(NA_real_, ncol(G)); names(r) <- colnames(G)
  if (any(!na_col)) r[!na_col] <- as.vector(cor(apply(G[, !na_col, drop = FALSE], 2, rank), rank(tfv)))
  for (j in which(na_col)) { ok <- !is.na(G[, j]); if (sum(ok) >= 10) r[j] <- cor(rank(tfv[ok]), rank(G[ok, j])) }
  r
}
bg_dm <- rbindlist(lapply(pair_order, function(pp) {
  d <- dm[pair == pp]; tfv <- d[[tf_of[[pp]]]]; stopifnot(!anyNA(tfv))
  G <- as.matrix(d[, setdiff(genes_all, tf_of[[pp]]), with = FALSE])
  zv <- apply(G, 2, function(v) isTRUE(var(v, na.rm = TRUE) == 0)); nmiss <- sum(colSums(is.na(G)) > 0)
  r <- allrank_cor(tfv, G[, !zv, drop = FALSE])
  data.table(pair = pp, symbol = names(r), rho = r, n_models = nrow(d), n_genes_with_missing = nmiss, n_zero_var = sum(zv))
}))
Go <- M[, setdiff(colnames(M), "SOX9"), drop = FALSE]
zvo <- apply(Go, 2, function(v) isTRUE(var(v, na.rm = TRUE) == 0))
ro <- allrank_cor(M[, "SOX9"], Go[, !zvo, drop = FALSE])
bg_o <- data.table(pair = "CRC", symbol = names(ro), rho = ro, n_models = nrow(M), n_genes_with_missing = sum(colSums(is.na(Go)) > 0), n_zero_var = sum(zvo))
fwrite(bg_dm, file.path(res_root, "background", "depmap_tf_allgene_rho.tsv.gz"), sep = "\t")
fwrite(bg_o, file.path(res_root, "background", "organoid_tf_allgene_rho.tsv.gz"), sep = "\t")

tc1 <- fread(file.path(out, "TCGA_TF_BAF_correlation.tsv"))   # v0.1, frozen
perc <- function(bgd, obs, src, universe_col = NULL) {
  rbindlist(lapply(unique(obs$pair), function(pp) {
    b <- bgd[pair == pp]; uni <- if (is.null(universe_col)) b else b[get(universe_col) == TRUE]
    v <- uni$rho[is.finite(uni$rho)]
    o <- obs[pair == pp]
    o[, .(pair, BAF, source = src, rho, signed_percentile = vapply(rho, pctl, 0, bgv = v),
          abs_percentile = vapply(abs(rho), pctl, 0, bgv = abs(v)), n_background = length(v),
          n_background_dropped_na = sum(!is.finite(uni$rho)),
          BAF_in_universe = if (is.null(universe_col)) BAF %in% b$symbol else BAF %in% uni$symbol,
          background_median_rho = median(v), exploratory_only = src == "DepMap" & pp == "Prostate")]
  }))
}
# consistency: the background rho for BAF members must reproduce the v0.1 rho
chk_t <- merge(bg_tcga[symbol %in% baf$symbol, .(pair, BAF = symbol, rho_bg = rho)], tc1[, .(pair, BAF, rho)], by = c("pair", "BAF"))
chk_d <- merge(bg_dm[symbol %in% baf$symbol, .(pair, BAF = symbol, rho_bg = rho)], dep1[, .(pair, BAF, rho)], by = c("pair", "BAF"))
chk_o <- merge(bg_o[symbol %in% baf$symbol, .(BAF = symbol, rho_bg = rho)], om[, .(BAF, rho)], by = "BAF")
stopifnot(nrow(chk_t) == 140, nrow(chk_d) == 140, nrow(chk_o) == 28,
          max(abs(chk_t$rho_bg - chk_t$rho)) < 1e-8, max(abs(chk_d$rho_bg - chk_d$rho)) < 1e-8, max(abs(chk_o$rho_bg - chk_o$rho)) < 1e-8)
bgp <- rbind(perc(bg_tcga, tc1[, .(pair, BAF, rho)], "TCGA", "in_universe"),
             perc(bg_dm, dep1[, .(pair, BAF, rho)], "DepMap"),
             perc(bg_o, om[, .(pair, BAF, rho)], "Organoid"))
bgp[, n_samples := fcase(source == "TCGA", bg_tcga$n_samples[match(pair, bg_tcga$pair)],
                         source == "DepMap", bg_dm$n_models[match(pair, bg_dm$pair)], default = nrow(M))]

# ==== write tables ====
fwrite(om[order(match(class, cls_order), BAF), .(pair, TF, BAF, class, n, rho, ci_low, ci_high, p, q, call, rho_depth_adjusted, q_depth_adjusted,
                                                 call_depth_adjusted, rho_change, depth_flag, TF_median_LFC, BAF_median_LFC, TF_sd_LFC, BAF_sd_LFC,
                                                 TF_prev_depleted, BAF_prev_depleted, pattern, negative_flag)],
       file.path(out, "Organoid_TF_BAF_codependency.tsv"), sep = "\t")
fwrite(osens[order(analysis, match(BAF, baf$symbol))], file.path(out, "Organoid_TF_BAF_sensitivity.tsv"), sep = "\t")
fwrite(ctx[order(match(pair, pair_order), match(class, cls_order), BAF)], file.path(out, "DepMap_TF_BAF_context.tsv"), sep = "\t")
fwrite(bgp[order(source, match(pair, pair_order), BAF)], file.path(out, "TF_BAF_background_percentile.tsv"), sep = "\t")
fwrite(inv, file.path(qc, "organoid_inventory.tsv"), sep = "\t")
fwrite(o_avail, file.path(qc, "organoid_availability.tsv"), sep = "\t")
fwrite(octrl, file.path(qc, "organoid_positive_control.tsv"), sep = "\t")
fwrite(lib_qc[gene %in% c("SOX9", baf$symbol)][order(gene, library)], file.path(qc, "organoid_library_LFC.tsv"), sep = "\t")
fwrite(unique(bgp[, .(source, pair, n_samples, n_background, n_background_dropped_na, background_median_rho)]),
       file.path(qc, "background_universe.tsv"), sep = "\t")

# ==== §6 figures (v0.2 files; v0.1 figures untouched) ====
pl <- rbind(tc1[, .(pair, BAF, class, rho, q, call, source = "TCGA expr")],
            dep1[, .(pair, BAF, class, rho, q, call, source = "DepMap 2D CRISPR")],
            om[, .(pair, BAF, class, rho, q, call, source = "Organoid CRISPR")],
            CJ(pair = pair_order[-1], BAF = baf$symbol)[, `:=`(class = baf$class[match(BAF, baf$symbol)], rho = NA_real_, q = NA_real_,
                                                              call = "unavailable", source = "Organoid CRISPR")])
pl[, pair_lab := factor(paste0(pair, "\n", tf_of[pair]), levels = paste0(pair_order, "\n", tf_of))]
pl[, BAF := factor(BAF, levels = rev(baf[order(match(class, cls_order))]$symbol))]
pl[, source := factor(source, levels = c("TCGA expr", "DepMap 2D CRISPR", "Organoid CRISPR"))]
pl[, sig := call %in% c("strong", "moderate", "strong_negative", "moderate_negative")]
ds <- om[depth_flag != "", BAF]
cap <- sprintf(paste0("Organoid: CRC only (n = %d individuals; DepMap 3D organoid biobank LFC, analysed separately, never merged with 2D gene effects); ",
                      "grey = unavailable (no matching organoids).\nDepth-sensitive organoid calls: %s. Prostate DepMap n = 10 (underpowered). Association only; not binding or recruitment."),
               nrow(u), if (length(ds)) paste(ds, collapse = ", ") else "none")
p2 <- ggplot(pl[!is.na(rho)], aes(source, BAF)) +
  geom_point(data = pl[is.na(rho)], shape = 4, size = 0.8, colour = "grey75") +
  geom_point(aes(fill = pmax(pmin(rho, 0.6), -0.6), size = pmin(-log10(q), 10), colour = sig), shape = 21, stroke = 0.6) +
  facet_grid(class ~ pair_lab, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-0.6, 0.6), name = "Spearman ρ") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), name = "|ρ| ≥ 0.2 & q < 0.05") +
  scale_size_continuous(range = c(0.6, 4), name = "−log10 q (cap 10)") +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, title = "TF–BAF association: TCGA co-expression, DepMap 2D co-dependency, organoid co-dependency (v0.2)", caption = cap) +
  theme_bw(base_size = 7) + theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_text(angle = 0),
                                  plot.caption = element_text(hjust = 0, size = 5.5))
ggsave(file.path(fig, "TFBAF2_summary_dotplot_v02.pdf"), p2, width = 9, height = 7.5, device = cairo_pdf)

sel <- fread(file.path(out, "TFBAF1_scatter_selection.tsv"))[pair == "CRC"]
pts <- rbindlist(lapply(seq_len(nrow(sel)), function(i) data.table(BAF = sel$BAF[i], why = sel$why[i], tfv = M[, "SOX9"], bv = M[, sel$BAF[i]])))
pts[, col := factor(sprintf("%s (%s)", BAF, why), levels = sprintf("%s (%s)", sel$BAF, sel$why))]
lab <- merge(sel, om[, .(BAF, rho, q, n, call, depth_flag, rho_depth_adjusted)], by = "BAF")
lab[, col := factor(sprintf("%s (%s)", BAF, why), levels = levels(pts$col))]
lab[, label := sprintf("ρ=%.2f q=%.1e n=%d\ndepth-adj ρ=%.2f%s%s", rho, q, n, rho_depth_adjusted,
                       ifelse(call == "none", "\nnonsignificant", ""), ifelse(depth_flag != "", paste0("\n", depth_flag), ""))]
p1 <- ggplot(pts, aes(tfv, bv)) + geom_point(size = 0.7, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "firebrick") +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.1, size = 2.2, inherit.aes = FALSE) +
  facet_wrap(~col, scales = "free", nrow = 1) +
  labs(x = "SOX9 LFC (organoid CRISPR)", y = "BAF member LFC", title = "CRC organoids: SOX9 vs v0.1-selected BAF members (trend line for display only)") +
  theme_bw(base_size = 7)
ggsave(file.path(fig, "TFBAF1_scatter_CRC_organoid_v02.pdf"), p1, width = 8, height = 3, device = cairo_pdf)

pat_lv <- c("both dependent, correlated", "both dependent, no positive co-dependency", "TF-specific dependency", "BAF-specific dependency", "neither dependent")
pg <- rbind(ctx[, .(pair, BAF, class, pattern, negative_flag, source = "DepMap 2D (effect < −0.5)")],
            om[, .(pair, BAF, class, pattern, negative_flag, source = "Organoid (is_depleted)")])
pg[, pair_lab := factor(paste0(pair, "\n", tf_of[pair]), levels = paste0(pair_order, "\n", tf_of))]
pg[, BAF := factor(BAF, levels = rev(baf[order(match(class, cls_order))]$symbol))]
pg[, pattern := factor(pattern, levels = pat_lv)]
p3 <- ggplot(pg, aes(pair_lab, BAF, fill = pattern)) + geom_tile(colour = "white") +
  geom_point(data = pg[negative_flag == TRUE], shape = 45, size = 4, show.legend = FALSE) +
  facet_grid(class ~ source, scales = "free", space = "free") +
  scale_fill_manual(values = setNames(c("#B2182B", "#F4A582", "#4393C3", "#92C5DE", "grey90"), pat_lv), drop = FALSE, name = NULL) +
  labs(x = NULL, y = NULL, title = sprintf("Dependency pattern per TF–BAF pair (dependent = ≥ %d%% of models); '−' = significant negative ρ", 100 * DEP_PREV),
       caption = "Descriptive. 2D and organoid dependency definitions differ and are not directly comparable. Prostate 2D n = 10.") +
  theme_bw(base_size = 7) + theme(strip.text.y = element_text(angle = 0))
ggsave(file.path(fig, "TFBAF3_dependency_pattern_v02.pdf"), p3, width = 7.5, height = 7, device = cairo_pdf)

writeLines(capture.output(sessionInfo()), file.path(res_root, "report", "sessionInfo_11_v02_analyze.txt"))
writeLines(system("git rev-parse HEAD", intern = TRUE), file.path(res_root, "report", "code_commit_11.txt"))

cat("\n== Organoid CRC SOX9 ==\n"); print(om[order(-rho), .(BAF, n, rho = round(rho, 2), ci = sprintf("[%.2f, %.2f]", ci_low, ci_high), q = signif(q, 2), call,
                                                          adj = round(rho_depth_adjusted, 2), depth_flag, BAF_prev = round(BAF_prev_depleted, 2), pattern)])
cat("SOX9 organoid prevalence:", round(om$TF_prev_depleted[1], 2), "\n")
print(octrl); print(osens[call != "none", .(analysis, BAF, rho = round(rho, 2), q = signif(q, 2), call)])
cat("\n== DepMap patterns ==\n"); print(ctx[, .N, by = .(pair, pattern)][order(match(pair, pair_order))])
cat("\n== Background (BAF with v0.1/v0.2 call) ==\n")
calls <- rbind(tc1[, .(source = "TCGA", pair, BAF, call)], dep1[, .(source = "DepMap", pair, BAF, call)], om[, .(source = "Organoid", pair, BAF, call)])
print(merge(bgp, calls, by = c("source", "pair", "BAF"))[call != "none", .(source, pair, BAF, rho = round(rho, 2), call,
                                                                           signed_pct = round(signed_percentile, 1), abs_pct = round(abs_percentile, 1), bg_med = round(background_median_rho, 2))])
print(unique(bgp[, .(source, pair, n_background, background_median_rho = round(background_median_rho, 3))]))
