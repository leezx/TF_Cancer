#!/usr/bin/env Rscript
# Step 2 design §5 figures: Fig3A per-cohort TF co-expression heatmaps, Fig3B recurrent
# TF network, FigS proliferation-coupling heatmap, focused Fig 1B view.
# Usage (repo root): Rscript scripts/02_tf_coexpression/04_plot.R
source("scripts/02_tf_coexpression/common2.R")
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(igraph) })
fig <- file.path(repo_out, "figures")
nodes <- fread(file.path(repo_out, "TF_network_modules.tsv"))
pa <- fread(file.path(repo_out, "TF_proliferation_association.tsv"))
ps <- fread(file.path(repo_out, "TF_proliferation_summary.tsv"))
foc <- fread(file.path(repo_out, "TF_focused_fig1B.tsv"))
de <- fread("results/01_tcga_tf_de/TCGA_TF_DE_all.tsv")
col_order <- c("COADREAD", "STAD", "LIHC", "PAAD", "LUAD", "LUSC", "HNSC", "KIRC", "KIRP", "KICH",
               "BRCA", "PRAD", "THCA", "GBM", "LGG", "LAML")
rho_col <- colorRamp2(c(-0.8, 0, 0.8), c("#2166AC", "white", "#B2182B"))

# ---- Fig3A: per cohort, top 40 network TFs present in the cohort ----
top_nodes <- nodes[order(-weighted_degree)]
for (k in cohorts) {
  m <- readRDS(file.path(res_dir, "coexpr", sprintf("%s.rds", k)))
  an <- if (m$primary_eligible) "primary" else "prolif_only"
  tfs <- head(top_nodes[TF %in% m$tf, TF], 40)
  if (length(tfs) < 3) next
  r <- m$mats[[an]][tfs, tfs]; diag(r) <- NA
  mod <- factor(nodes$module[match(tfs, nodes$TF)])
  h <- Heatmap(r, name = "partial\nSpearman ρ", col = rho_col, na_col = "grey90",
               row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
               top_annotation = HeatmapAnnotation(module = mod, show_legend = FALSE, simple_anno_size = unit(2, "mm")),
               column_title = sprintf("%s: TF co-expression across tumors (%s), top recurrent-network TFs", k,
                                      if (an == "primary") "purity + proliferation adjusted" else "proliferation-only; not primary-eligible"),
               column_title_gp = gpar(fontsize = 7))
  cairo_pdf(file.path(fig, sprintf("Fig3A_%s_TF_coexpression_heatmap.pdf", k)), width = 6, height = 5.6); draw(h); dev.off()
}

# ---- Fig3B: recurrent-edge network (primary) ----
rec <- fread(file.path(res_dir, "report", "TF_pair_recurrence.tsv.gz"))[analysis == "primary"]
g <- graph_from_data_frame(rec[, .(TF_a, TF_b, weight = n_organ_systems_positive)], directed = FALSE,
                           vertices = nodes[, .(TF, module, weighted_degree, module_reported)])
keep <- V(g)$module_reported
g2 <- induced_subgraph(g, V(g)[keep])
mods <- sort(unique(V(g2)$module))
pal <- setNames(grDevices::hcl.colors(max(3, length(mods)), "Dark 3")[seq_along(mods)], mods)
lab_nodes <- nodes[module_reported == TRUE][order(-weighted_degree), head(TF, 4), by = module]$V1
set.seed(P$seed)
lay <- layout_with_fr(g2, weights = E(g2)$weight)
cairo_pdf(file.path(fig, "Fig3B_recurrent_TF_network.pdf"), width = 8, height = 8)
plot(g2, layout = lay, vertex.size = 1.5 + 4 * sqrt(V(g2)$weighted_degree / max(V(g2)$weighted_degree)),
     vertex.color = pal[as.character(V(g2)$module)], vertex.frame.color = NA,
     vertex.label = ifelse(V(g2)$name %in% lab_nodes, V(g2)$name, NA), vertex.label.cex = 0.45, vertex.label.color = "black",
     edge.width = 0.2 + 0.3 * (E(g2)$weight - 2), edge.color = adjustcolor("grey50", 0.25))
title(sprintf("Recurrent TF co-expression network (purity + proliferation adjusted; ≥%d organ systems)\n%d TFs, %d edges, %d modules ≥ %d TFs; co-expression, not binding",
              P$recur_organ_systems, vcount(g2), ecount(g2), length(mods), P$module_min), cex.main = 0.7)
dev.off()

# ---- FigS: proliferation coupling of Step 1 Fig 1B TFs ----
tfs <- foc$TF
w <- dcast(pa[analysis == "purity_adjusted" & TF %in% tfs], TF ~ cohort, value.var = "rho")
cols <- intersect(col_order, names(w))
M <- as.matrix(w[, ..cols]); rownames(M) <- w$TF
M <- M[intersect(foc[order(-n_organ_systems_C1_up)]$TF, rownames(M)), , drop = FALSE]
# rows ordered by mean rho (TFs absent from some cohorts leave NA cells, so no hclust)
M <- M[order(-rowMeans(M, na.rm = TRUE)), , drop = FALSE]
rc <- foc$recurrently_proliferation_coupled[match(rownames(M), foc$TF)]
h <- Heatmap(M, name = "ρ(TF, prolif)\npurity-adj.", col = colorRamp2(c(-0.8, 0, 0.8), c("#2166AC", "white", "#B2182B")), na_col = "grey90",
             cluster_columns = FALSE, show_row_names = FALSE, cluster_rows = FALSE,
             left_annotation = rowAnnotation(recurrently_coupled = rc, col = list(recurrently_coupled = c(`TRUE` = "black", `FALSE` = "grey85"))),
             column_title = sprintf("Proliferation association of the %d Step 1 Fig 1B TFs (purity-adjusted; KIRC/LAML not eligible)", nrow(M)),
             column_title_gp = gpar(fontsize = 7))
cairo_pdf(file.path(fig, "FigS_TF_proliferation_coupling_heatmap.pdf"), width = 6, height = 8); draw(h); dev.off()

# ---- focused Fig 1B: Step 1 C1 log2FC of Fig 1B TFs not recurrently proliferation-coupled ----
ft <- foc[in_focused_view == TRUE, TF]
x <- de[contrast == "C1" & TF %in% ft]
cc <- intersect(col_order, unique(x$cohort))
L <- as.matrix(dcast(x, TF ~ cohort, value.var = "log2FC_shrunk")[, c("TF", cc), with = FALSE], rownames = "TF")
Ti <- as.matrix(dcast(x, TF ~ cohort, value.var = "tier")[, c("TF", cc), with = FALSE], rownames = "TF")
h <- Heatmap(L, name = "shrunk\nlog2FC (C1)", col = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")), na_col = "grey92",
             cluster_columns = FALSE, show_row_names = nrow(L) <= 150, row_names_gp = gpar(fontsize = 4),
             column_title = sprintf("Focused Fig 1B (secondary): %d C1-recurrent TFs not recurrently proliferation-coupled (Tier 1 ● / Tier 2 ○)", nrow(L)),
             column_title_gp = gpar(fontsize = 7),
             cell_fun = function(j, i, x0, y0, w0, h0, f) if (!is.na(Ti[i, j])) grid.points(x0, y0, pch = if (Ti[i, j] == "Tier1") 16 else 1, size = unit(0.9, "mm")))
cairo_pdf(file.path(fig, "Fig1B_focused_nonproliferative_TFs.pdf"), width = 6.5, height = max(4, min(14, 1.5 + nrow(L) * 0.07))); draw(h); dev.off()
log_msg("figures written: Fig3A x%d, Fig3B (%d TFs), FigS (%d TFs), focused Fig1B (%d TFs)",
        length(list.files(fig, "^Fig3A_")), vcount(g2), nrow(M), nrow(L))
