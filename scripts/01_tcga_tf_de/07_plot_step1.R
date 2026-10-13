#!/usr/bin/env Rscript
# Design §1.6 / §3.6: draft Figs 1A, 1B (primary C1 + supplementary exploratory)
# and Fig 2 (UpSet + dot plot). Labeled "TF expression discovery" (review N9).
# Usage (repo root): Rscript scripts/01_tcga_tf_de/07_plot_step1.R
source("scripts/01_tcga_tf_de/common.R")
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(ggplot2); library(UpSetR) })
fig <- file.path(repo_results, "figures"); dir.create(fig, showWarnings = FALSE)
de <- fread(file.path(repo_results, "TCGA_TF_DE_all.tsv"))
rec <- fread(file.path(repo_results, "TF_recurrent_upregulation.tsv"))
tfr <- fread(file.path(repo_results, "TF_reference.tsv"))
cm <- fread(file.path(repo_results, "cohort_manifest.tsv"))
col_order <- c("COADREAD", "STAD", "LIHC", "PAAD", "LUAD", "LUSC", "HNSC", "KIRC", "KIRP", "KICH",
               "BRCA", "PRAD", "THCA", "GBM", "LGG")
lfc_col <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))

heat <- function(tfs, ct, file, title) {
  x <- de[contrast == ct & TF %in% tfs]
  cols <- intersect(col_order, unique(x$cohort))
  m <- as.matrix(dcast(x, TF ~ cohort, value.var = "log2FC_shrunk")[, c("TF", cols), with = FALSE], rownames = "TF")
  t <- as.matrix(dcast(x, TF ~ cohort, value.var = "tier")[, c("TF", cols), with = FALSE], rownames = "TF")
  dbd <- tfr$DBD[match(rownames(m), tfr$hgnc_symbol)]
  top <- names(sort(table(dbd), decreasing = TRUE))[1:min(8, length(unique(dbd)))]
  dbd[!dbd %in% top] <- "Other"
  pal <- setNames(c(RColorBrewer::brewer.pal(max(3, length(top)), "Set2")[seq_along(top)], "grey80"), c(top, "Other"))
  ra <- rowAnnotation(DBD = dbd, col = list(DBD = pal))
  h <- Heatmap(m, name = "shrunk\nlog2FC", col = lfc_col, na_col = "grey92", right_annotation = ra,
               cluster_columns = FALSE, show_row_names = nrow(m) <= 120, row_names_gp = gpar(fontsize = 5),
               column_title = title, column_title_gp = gpar(fontsize = 9),
               cell_fun = function(j, i, x0, y0, w, h, f) {
                 if (!is.na(t[i, j])) grid.points(x0, y0, pch = if (t[i, j] == "Tier1") 16 else 1, size = unit(1.2, "mm"))
               })
  pdf(file, width = 6.5, height = max(4, min(14, 1.5 + nrow(m) * 0.07))); draw(h); dev.off()
  nrow(m)
}

up_c1 <- de[contrast == "C1" & !is.na(tier), unique(TF)]
n1 <- heat(up_c1, "C1", file.path(fig, "Fig1A_all_up_TFs_heatmap.pdf"),
           "TF expression discovery: TFs up (Tier 1 ● / Tier 2 ○) in ≥1 TCGA tumor-vs-normal (C1) cohort")
n2 <- heat(rec[fig1B_primary == TRUE, TF], "C1", file.path(fig, "Fig1B_recurrent_TFs_heatmap.pdf"),
           sprintf("TF expression discovery: C1-up in ≥%d organ systems", cfg$thresholds$fig1b_min_organ_systems))
n3 <- heat(rec[fig1B_supplementary == TRUE, TF], "C2", file.path(fig, "Fig1B_supp_exploratory_TFs_heatmap_C2.pdf"),
           "Exploratory (cross-source, TCGA vs GTEx): up via C1 or C2 in ≥3 organ systems; shown C2 LFC")

# Fig 2A: UpSet of C1 up sets
sets <- lapply(split(de[contrast == "C1" & !is.na(tier), TF], de[contrast == "C1" & !is.na(tier), cohort]), unique)
pdf(file.path(fig, "Fig2_TF_UpSet.pdf"), width = 9, height = 5, onefile = FALSE)
print(upset(fromList(sets), nsets = length(sets), nintersects = 40, order.by = "freq",
            mainbar.y.label = "TFs up (C1, Tier 1/2)", sets.x.label = "TFs up per cohort"))
dev.off()

# Fig 2B: dot plot of top recurrent TFs (C1)
top <- head(rec[n_organ_systems_C1_up > 0], 40)$TF
x <- de[contrast == "C1" & TF %in% top]
x[, TF := factor(TF, levels = rev(top))]
x[, cohort := factor(cohort, levels = intersect(col_order, unique(cohort)))]
p <- ggplot(x, aes(cohort, TF)) +
  geom_point(aes(size = tumor_prevalence, fill = pmax(pmin(log2FC_shrunk, 3), -3),
                 colour = !is.na(padj) & padj < 0.05), shape = 21, stroke = 0.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "shrunk log2FC") +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"), name = "padj < 0.05") +
  scale_size_continuous(range = c(0.5, 3.5), name = "tumor prevalence\n(TPM ≥ 1)") +
  labs(x = NULL, y = NULL, title = "TF expression discovery: top 40 C1-recurrent TFs") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig, "Fig2_TF_dotplot.pdf"), p, width = 6, height = 7)
log_msg("Fig1A %d TFs; Fig1B primary %d; Fig1B supp %d; UpSet sets %d", n1, n2, n3, length(sets))
