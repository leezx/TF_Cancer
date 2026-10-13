#!/usr/bin/env Rscript
# Manhattan-style TF–BAF correlation plot: one point per BAF member, one column per
# TF–cancer pair (CRC, AML, Prostate, Breast, Melanoma). Three stacked panels:
# TCGA expression, DepMap 2D CRISPR dependency, organoid CRISPR dependency.
# Top significant (q < 0.05) positive members are labelled; negatives are not.
# Input: results/tfbaf/TF_BAF_core_correlation.tsv
# Usage (repo root): Rscript scripts/tfbaf/12_manhattan_plot.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
out <- "results/tfbaf"; fig <- file.path(out, "figures")
TOP_N <- 3   # labels per column: top significant positive members by rho
x <- fread(file.path(out, "TF_BAF_core_correlation.tsv"))
pair_order <- c("CRC", "AML", "Prostate", "Breast", "Melanoma")
tf_of <- setNames(c("SOX9", "SPI1", "FOXA1", "GATA3", "SOX10"), pair_order)
cls_order <- c("ATPase", "core", "cBAF-specific", "PBAF-specific", "ncBAF-specific")
cls_col <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00"), cls_order)
lab_x <- paste0(pair_order, "\n", tf_of)

long <- rbindlist(lapply(c("TCGA", "DepMap", "Organoid"), function(s)
  x[, .(pair, BAF, class, source = s, n = get(paste0(s, "_n")), rho = get(paste0(s, "_rho")), q = get(paste0(s, "_q")))]))
long[, col := match(pair, pair_order)]
# deterministic jitter: members spread across the column in BAF-list order
long[, xpos := col + (match(BAF, unique(x$BAF)) - 14.5) / 14.5 * 0.32]
long[, sig := !is.na(q) & q < 0.05]
long[, class := factor(class, levels = cls_order)]
long[, label := ""]
long[sig & rho > 0, label := fifelse(frank(-rho, ties.method = "first") <= TOP_N, BAF, ""), by = .(source, pair)]

titles <- c(TCGA = "TCGA tumors: TF–BAF expression correlation (log2 TPM+1)",
            DepMap = "DepMap 2D cell lines: TF–BAF dependency correlation (CRISPR gene effect)",
            Organoid = "DepMap 3D organoids: TF–BAF dependency correlation (CRISPR LFC)")
panel <- function(s) {
  d <- long[source == s]
  nn <- d[, .(n = n[1]), by = .(pair, col)]
  nn[, txt := fifelse(is.na(n), "not available", paste0("n = ", n))]
  ggplot(d[!is.na(rho)], aes(xpos, rho)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 1.5:4.5, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(colour = class, shape = sig), size = 1.9, stroke = 0.6) +
    geom_text_repel(aes(label = label, colour = class), size = 2.4, min.segment.length = 0, segment.size = 0.25,
                    box.padding = 0.3, max.overlaps = Inf, seed = 1, show.legend = FALSE) +
    geom_text(data = nn, aes(x = col, y = -0.85, label = txt), size = 2.2, colour = "grey35", inherit.aes = FALSE) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."), name = NULL, drop = FALSE) +
    scale_colour_manual(values = cls_col, name = "BAF class", drop = FALSE) +
    scale_x_continuous(breaks = 1:5, labels = lab_x, limits = c(0.5, 5.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.9, 0.9), breaks = seq(-0.8, 0.8, 0.4)) +
    labs(x = NULL, y = "Spearman ρ", title = titles[[s]]) +
    theme_bw(base_size = 8) + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                                    plot.title = element_text(size = 8, face = "bold"))
}
p <- (panel("TCGA") / panel("DepMap") / panel("Organoid")) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = sprintf(paste0("Each point is one of 28 BAF members. Labels: top %d significant positive members per column (BH within pair). ",
                                           "Prostate DepMap n = 10 (underpowered). Association only; not binding or recruitment."), TOP_N),
                  theme = theme(plot.caption = element_text(size = 6, hjust = 0)))
ggsave(file.path(fig, "TFBAF_manhattan_rho.pdf"), p, width = 7.5, height = 9, device = cairo_pdf)
ggsave(file.path(fig, "TFBAF_manhattan_rho.png"), p, width = 7.5, height = 9, dpi = 200)
print(long[label != "", .(source, pair, BAF, rho = round(rho, 2), q = signif(q, 2))][order(source, match(pair, pair_order), -rho)])
