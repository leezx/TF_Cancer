#!/usr/bin/env Rscript
# TF–BAF design §3.2–§3.5 and §4: TCGA TF–BAF Spearman (main + sensitivities),
# DepMap 26Q1 CRISPR co-dependency in lineage-matched cell lines, QC, figures.
# Usage (repo root): Rscript scripts/tfbaf/02_analyze.R   (after 01_extract_tcga.R)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
res_root <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx", "result", "2026-10-13_tf_cancer_tfbaf_v0.1")
out <- "results/tfbaf"; qc <- file.path(out, "qc"); fig <- file.path(out, "figures")
for (d in c(qc, fig, file.path(res_root, "report"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
baf <- fread(file.path(out, "BAF_members.tsv"))
cls_order <- c("ATPase", "core", "cBAF-specific", "PBAF-specific", "ncBAF-specific")
pair_order <- c("CRC", "AML", "Prostate", "Breast", "Melanoma")

call_fun <- function(rho, q) fcase(is.na(rho) | is.na(q) | q >= 0.05, "none",
                                   rho >= 0.3, "strong", rho >= 0.2, "moderate",
                                   rho <= -0.3, "strong_negative", rho <= -0.2, "moderate_negative", default = "none")
spear <- function(a, b) {
  ok <- complete.cases(a, b)
  if (sum(ok) < 10) return(list(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman", exact = FALSE))
  list(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
# partial Spearman: rank, residualize on covariate ranks, Pearson of residuals (df = n - 2 - k)
pspear <- function(a, b, Z) {
  ok <- complete.cases(a, b, Z); n <- sum(ok)
  if (n < 10) return(list(n = n, rho = NA_real_, p = NA_real_))
  Zr <- apply(as.matrix(Z)[ok, , drop = FALSE], 2, rank)
  ra <- resid(lm(rank(a[ok]) ~ Zr)); rb <- resid(lm(rank(b[ok]) ~ Zr))
  r <- cor(ra, rb); k <- ncol(Zr); df <- n - 2 - k
  list(n = n, rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df))
}

# ---- TCGA ----
x <- fread(file.path(res_root, "tcga", "tcga_tf_baf_log2TPM1.tsv"))
pur <- fread(file.path(root, "1.Databases", "TCGA_PanCanAtlas_ABSOLUTE_purity", "raw", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"))
pur <- pur[`call status` == "called", .(sample = array, purity)]
x <- merge(x, pur, by = "sample", all.x = TRUE)
st <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(root, "1.Databases", "TCGA_PanCanAtlas_subtypes", "raw", "TCGASubtype.20170308.tsv.gz"))))
x <- merge(x, st[, .(sample = sampleID, pam50 = Subtype_mRNA)], by = "sample", all.x = TRUE)
lin2tpm <- function(v) 2^v - 1

tcga_cor <- function(d, set_label, adjust = NULL) {
  rbindlist(lapply(baf$symbol, function(g) {
    tfv <- d[[d$TF[1]]]
    r <- if (is.null(adjust)) spear(tfv, d[[g]]) else pspear(tfv, d[[g]], d[, ..adjust])
    data.table(BAF = g, n = r$n, rho = r$rho, p = r$p,
               TF_median_TPM = median(lin2tpm(tfv)), BAF_median_TPM = median(lin2tpm(d[[g]])))
  }))[, `:=`(analysis = set_label, q = p.adjust(p, "BH"))]
}
main <- x[set == "main", tcga_cor(.SD, "main"), by = .(pair, TF, cohort), .SDcols = names(x)]
main[, call := call_fun(rho, q)]
main[, low_expression := BAF_median_TPM < 1]
sens <- rbind(
  x[set == "main" & !is.na(purity), tcga_cor(.SD, "purity_adjusted", "purity"), by = .(pair, TF, cohort), .SDcols = names(x)],
  x[set == "main" & cohort == "BRCA" & pam50 %in% c("LumA", "LumB"), tcga_cor(.SD, "BRCA_luminal_only"), by = .(pair, TF, cohort), .SDcols = names(x)],
  x[set == "skcm_metastatic", tcga_cor(.SD, "SKCM_metastatic"), by = .(pair, TF, cohort), .SDcols = names(x)])
sens[, call := call_fun(rho, q)]
# flag purity reversal / strong attenuation (review Q4)
pa <- merge(main[, .(pair, BAF, rho_main = rho, call_main = call)], sens[analysis == "purity_adjusted", .(pair, BAF, rho_adj = rho)], by = c("pair", "BAF"))
pa[, purity_flag := fcase(call_main != "none" & sign(rho_adj) != sign(rho_main), "reversed",
                          call_main != "none" & abs(rho_adj) < 0.5 * abs(rho_main), "attenuated_gt50pct", default = "")]
main <- merge(main, pa[, .(pair, BAF, rho_purity_adjusted = rho_adj, purity_flag)], by = c("pair", "BAF"), all.x = TRUE)
main <- merge(main, baf[, .(BAF = symbol, class)], by = "BAF")
fwrite(main[order(match(pair, pair_order), match(class, cls_order), BAF)], file.path(out, "TCGA_TF_BAF_correlation.tsv"), sep = "\t")
fwrite(sens[order(analysis, match(pair, pair_order), BAF)], file.path(out, "TCGA_TF_BAF_sensitivity.tsv"), sep = "\t")

# ---- DepMap CRISPR co-dependency ----
dd <- file.path(root, "1.Databases", "DepMap3D_202608")
hdr <- names(fread(file.path(dd, "CRISPRGeneEffect_DepMap26Q1.csv"), nrows = 0))
want <- c(pair_order_tf <- c("SOX9", "SPI1", "FOXA1", "GATA3", "SOX10"), baf$symbol)
sym <- sub(" \\(.*$", "", hdr)
gcols <- hdr[sym %in% want]
stopifnot(length(gcols) == length(want))
ce <- fread(file.path(dd, "CRISPRGeneEffect_DepMap26Q1.csv"), select = c(hdr[1], gcols))
setnames(ce, c("ModelID", sub(" \\(.*$", "", gcols)))
mod <- fread(file.path(dd, "Model_DepMap26Q1.csv"))
ce <- merge(ce, mod[, .(ModelID, CellLineName, OncotreeLineage, OncotreePrimaryDisease, ModelType)], by = "ModelID")
stopifnot(!anyDuplicated(ce$ModelID))
ce[, pair := fcase(OncotreeLineage == "Bowel", "CRC", OncotreePrimaryDisease == "Acute Myeloid Leukemia", "AML",
                   OncotreeLineage == "Prostate", "Prostate", OncotreeLineage == "Breast", "Breast",
                   OncotreePrimaryDisease %like% "Melanoma", "Melanoma", default = NA_character_)]
dm <- ce[!is.na(pair) & ModelType == "Cell Line"]
tf_of <- setNames(c("SOX9", "SPI1", "FOXA1", "GATA3", "SOX10"), pair_order)
dep <- dm[, {
  tfv <- .SD[[tf_of[[.BY$pair]]]]
  rbindlist(lapply(baf$symbol, function(g) { r <- spear(tfv, .SD[[g]])
    data.table(BAF = g, n = r$n, rho = r$rho, p = r$p, TF_median_effect = median(tfv, na.rm = TRUE),
               TF_frac_lt_m0.5 = mean(tfv < -0.5, na.rm = TRUE), BAF_median_effect = median(.SD[[g]], na.rm = TRUE)) }))
}, by = pair, .SDcols = setdiff(names(dm), "pair")]
dep[, TF := tf_of[pair]]
dep[, q := p.adjust(p, "BH"), by = pair]
dep[, call := call_fun(rho, q)]
dep[, underpowered := n < 20]
dep[, TF_not_dependency := TF_median_effect > -0.3]
dep <- merge(dep, baf[, .(BAF = symbol, class)], by = "BAF")
fwrite(dep[order(match(pair, pair_order), match(class, cls_order), BAF)], file.path(out, "DepMap_TF_BAF_codependency.tsv"), sep = "\t")

# ---- summary ----
sm <- merge(main[, .(pair, TF, BAF, class, TCGA_n = n, TCGA_rho = rho, TCGA_q = q, TCGA_call = call, low_expression, purity_flag)],
            dep[, .(pair, BAF, DepMap_n = n, DepMap_rho = rho, DepMap_q = q, DepMap_call = call, underpowered, TF_not_dependency)],
            by = c("pair", "BAF"), all = TRUE)
sm[, both_positive := TCGA_call %in% c("strong", "moderate") & DepMap_call %in% c("strong", "moderate") & !underpowered]
fwrite(sm[order(match(pair, pair_order), match(class, cls_order), BAF)], file.path(out, "TF_BAF_summary.tsv"), sep = "\t")

# ---- QC (§4) ----
tfdep <- unique(dep[, .(pair, TF, n_models = n, TF_median_effect, TF_frac_lt_m0.5)], by = "pair")
tfexpr <- unique(main[, .(pair, TF, TF_median_TPM)], by = "pair")
fwrite(merge(tfexpr, tfdep, by = c("pair", "TF")), file.path(qc, "TF_expression_and_dependency.tsv"), sep = "\t")
allcl <- ce[ModelType == "Cell Line"]
ctrl <- rbindlist(lapply(list(c("SMARCB1", "SMARCE1"), c("ARID1A", "SMARCB1")), function(g) {
  r <- spear(allcl[[g[1]]], allcl[[g[2]]]); data.table(gene_a = g[1], gene_b = g[2], n = r$n, rho = r$rho, p = r$p) }))
fwrite(ctrl, file.path(qc, "BAF_codependency_positive_control.tsv"), sep = "\t")
fwrite(data.table(check = c("TCGA genes present", "CRISPR genes present"), n = c(length(intersect(c(want), names(x))), length(gcols)), expected = length(want)),
       file.path(qc, "gene_mapping.tsv"), sep = "\t")

# ---- Fig TFBAF-2: dual-evidence summary dot plot ----
pl <- rbind(main[, .(pair, BAF, class, rho, q, call, source = "TCGA expr")],
            dep[, .(pair, BAF, class, rho, q, call, source = "DepMap CRISPR")])
pl[, pair_lab := factor(paste0(pair, "\n", fcase(pair == "CRC", "SOX9", pair == "AML", "SPI1", pair == "Prostate", "FOXA1",
                                                    pair == "Breast", "GATA3", default = "SOX10")),
                        levels = paste0(pair_order, "\n", c("SOX9", "SPI1", "FOXA1", "GATA3", "SOX10")))]
pl[, BAF := factor(BAF, levels = rev(baf[order(match(class, cls_order))]$symbol))]
pl[, source := factor(source, levels = c("TCGA expr", "DepMap CRISPR"))]
pl[, sig := call %in% c("strong", "moderate", "strong_negative", "moderate_negative")]
p2 <- ggplot(pl, aes(source, BAF)) +
  geom_point(aes(fill = pmax(pmin(rho, 0.6), -0.6), size = pmin(-log10(q), 10), colour = sig), shape = 21, stroke = 0.6) +
  facet_grid(class ~ pair_lab, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-0.6, 0.6), name = "Spearman ρ") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), name = "|ρ| ≥ 0.2 & q < 0.05") +
  scale_size_continuous(range = c(0.6, 4), name = "−log10 q (cap 10)") +
  labs(x = NULL, y = NULL, title = "TF–BAF association: TCGA tumor co-expression and DepMap co-dependency",
       caption = "Prostate DepMap n = 10 (underpowered). Association only; not binding or recruitment.") +
  theme_bw(base_size = 7) + theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_text(angle = 0))
ggsave(file.path(fig, "TFBAF2_summary_dotplot.pdf"), p2, width = 8, height = 7, device = cairo_pdf)

# ---- Fig TFBAF-1: scatter panels, preset selection (SMARCA4 + top TCGA + top DepMap) ----
sel <- rbindlist(lapply(pair_order, function(pp) {
  a <- main[pair == pp][order(-rho)]$BAF; b <- dep[pair == pp][order(-rho)]$BAF
  s1 <- "SMARCA4"; s2 <- setdiff(a, s1)[1]; s3 <- setdiff(b, c(s1, s2))[1]
  data.table(pair = pp, BAF = c(s1, s2, s3), why = c("fixed reference", "top TCGA rho", "top DepMap rho"))
}))
fwrite(sel, file.path(out, "TFBAF1_scatter_selection.tsv"), sep = "\t")
pts <- rbindlist(lapply(seq_len(nrow(sel)), function(i) {
  pp <- sel$pair[i]; g <- sel$BAF[i]; tfg <- tf_of[[pp]]
  rbind(x[set == "main" & pair == pp, .(pair = pp, BAF = g, slot = i, source = "TCGA log2(TPM+1)", tfv = get(tfg), bv = get(g))],
        dm[pair == pp, .(pair = pp, BAF = g, slot = i, source = "DepMap gene effect", tfv = get(tfg), bv = get(g))])
}))
lab <- merge(sel[, .(pair, BAF, slot = seq_len(.N), why)],
             rbind(main[, .(pair, BAF, source = "TCGA log2(TPM+1)", rho, q, n, call)], dep[, .(pair, BAF, source = "DepMap gene effect", rho, q, n, call)]),
             by = c("pair", "BAF"))
lab[, label := sprintf("%s\nρ=%.2f q=%.1e n=%d%s", BAF, rho, q, n, ifelse(call == "none", "\nnonsignificant", ""))]
pts[, panel := factor(paste(pair, slot)), ]
cairo_pdf(file.path(fig, "TFBAF1_scatter.pdf"), width = 9, height = 4.2, onefile = TRUE)
for (pp in pair_order) {
  d <- pts[pair == pp]; l <- lab[pair == pp]
  d[, col := factor(paste0(BAF, " (", sel[pair == pp]$why[match(slot, which(sel$pair == pp))], ")"))]
  l[, col := factor(paste0(BAF, " (", why, ")"))]
  print(ggplot(d, aes(tfv, bv)) + geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "firebrick") +
    geom_text(data = l, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.1, size = 2.2, inherit.aes = FALSE) +
    facet_wrap(source ~ col, scales = "free", nrow = 2) +
    labs(x = paste(tf_of[[pp]], "(TF)"), y = "BAF member", title = sprintf("%s: %s vs BAF members (trend line for display only)", pp, tf_of[[pp]])) +
    theme_bw(base_size = 7))
}
dev.off()
writeLines(capture.output(sessionInfo()), file.path(res_root, "report", "sessionInfo_02_analyze.txt"))
writeLines(system("git rev-parse HEAD", intern = TRUE), file.path(res_root, "report", "code_commit.txt"))

cat("\n== TCGA main calls ==\n"); print(dcast(main, BAF ~ pair, value.var = "rho", fun.aggregate = function(v) round(v, 2))[order(match(BAF, baf$symbol))])
cat("\n== DepMap rho ==\n"); print(dcast(dep, BAF ~ pair, value.var = "rho", fun.aggregate = function(v) round(v, 2))[order(match(BAF, baf$symbol))])
print(merge(tfexpr, tfdep, by = c("pair", "TF"))); print(ctrl)
cat("both positive:\n"); print(sm[both_positive == TRUE, .(pair, TF, BAF, TCGA_rho, DepMap_rho)])
