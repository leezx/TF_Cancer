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
# background z (dataset-comparable): z = (rho - mean_bg) / sd_bg, where bg = the TF's rho with every
# gene in the same dataset and pair (TF excluded; TCGA: protein-coding, median TPM >= 1), from 11_v02_analyze.R
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
bgd <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx", "result", "2026-10-13_tf_cancer_tfbaf_v0.2", "background")
bg <- rbind(fread(file.path(bgd, "tcga_tf_allgene_rho.tsv.gz"))[in_universe == TRUE, .(source = "TCGA", pair, rho)],
            fread(file.path(bgd, "depmap_tf_allgene_rho.tsv.gz"))[, .(source = "DepMap", pair, rho)],
            fread(file.path(bgd, "organoid_tf_allgene_rho.tsv.gz"))[, .(source = "Organoid", pair, rho)])
bgs <- bg[is.finite(rho), .(bg_mean = mean(rho), bg_sd = sd(rho), bg_n = .N), by = .(source, pair)]
for (s in c("TCGA", "DepMap", "Organoid")) {
  b <- bgs[source == s]; i <- match(x$pair, b$pair)
  set(x, j = paste0(s, "_bg_mean"), value = b$bg_mean[i]); set(x, j = paste0(s, "_bg_sd"), value = b$bg_sd[i])
  set(x, j = paste0(s, "_z"), value = (x[[paste0(s, "_rho")]] - b$bg_mean[i]) / b$bg_sd[i])
}
fwrite(x, file.path(out, "TF_BAF_core_correlation.tsv"), sep = "\t", na = "NA")
fwrite(bgs, file.path(out, "qc", "background_mean_sd.tsv"), sep = "\t")
YVAR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "z"   # "z" (default) or "rho"
# optional 2nd arg "internal": organoid panel from DepMap internal 25Q2 3DO screens (not public; output to DATA only)
INTERNAL <- length(commandArgs(TRUE)) > 1 && commandArgs(TRUE)[2] == "internal"
if (INTERNAL) {
  idir <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx", "result", "2026-10-13_tf_cancer_tfbaf_v0.3_internal3D")
  io <- fread(file.path(idir, "Organoid_internal25Q2_TF_BAF_correlation.tsv"))
  ib <- fread(file.path(idir, "organoid_internal25Q2_background_mean_sd.tsv"))
  x[, c("Organoid_n", "Organoid_rho", "Organoid_q", "Organoid_z") := NA]
  i <- match(paste(x$pair, x$BAF), paste(io$pair, io$BAF))
  x[, `:=`(Organoid_n = io$n[i], Organoid_rho = io$rho[i], Organoid_q = io$q[i], Organoid_z = io$z[i])]
  fig <- idir
}
pair_order <- c("CRC", "AML", "Prostate", "Breast", "Melanoma")
tf_of <- setNames(c("SOX9", "SPI1", "FOXA1", "GATA3", "SOX10"), pair_order)
cls_order <- c("ATPase", "core", "cBAF-specific", "PBAF-specific", "ncBAF-specific")
cls_col <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00"), cls_order)
lab_x <- paste0(pair_order, "\n", tf_of)

long <- rbindlist(lapply(c("TCGA", "DepMap", "Organoid"), function(s)
  x[, .(pair, BAF, class, source = s, n = get(paste0(s, "_n")), rho = get(paste0(s, "_rho")), q = get(paste0(s, "_q")),
      z = get(paste0(s, "_z")))]))
long[, y := if (YVAR == "z") z else rho]
long[, col := match(pair, pair_order)]
# deterministic jitter: members spread across the column in BAF-list order
long[, xpos := col + (match(BAF, unique(x$BAF)) - 14.5) / 14.5 * 0.32]
long[, sig := factor(!is.na(q) & q < 0.05, levels = c(TRUE, FALSE))]
long[, class := factor(class, levels = cls_order)]
long[, label := ""]
long[sig == "TRUE" & rho > 0, label := fifelse(frank(-y, ties.method = "first") <= TOP_N, BAF, ""), by = .(source, pair)]

titles <- c(TCGA = "TCGA tumors: TF–BAF expression correlation (log2 TPM+1)",
            DepMap = "DepMap 2D cell lines: TF–BAF dependency correlation (CRISPR gene effect)",
            Organoid = if (INTERNAL) "DepMap 3D organoids (internal 25Q2, 3DO screens): TF–BAF dependency correlation (Chronos gene effect)" else "Sanger organoid biobank: TF–BAF dependency correlation (CRISPR LFC)")
ylim <- if (YVAR == "z") { m <- ceiling(max(abs(long$z), na.rm = TRUE)); c(-m - 0.5, m + 0.5) } else c(-0.9, 0.9)
panel <- function(s) {
  d <- long[source == s]
  nn <- d[, .(n = n[1]), by = .(pair, col)]
  nn[, txt := fifelse(is.na(n), "not available", paste0("n = ", n, fifelse(n < 10, " (underpowered)", "")))]
  ggplot(d[!is.na(y)], aes(xpos, y)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 1.5:4.5, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(colour = class, shape = sig), size = 1.9, stroke = 0.6) +
    geom_text_repel(aes(label = label, colour = class), size = 2.4, min.segment.length = 0, segment.size = 0.25,
                    box.padding = 0.3, max.overlaps = Inf, seed = 1, show.legend = FALSE) +
    geom_text(data = nn, aes(x = col, y = ylim[1] * 0.94, label = txt), size = 2.2, colour = "grey35", inherit.aes = FALSE) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."), name = NULL, drop = FALSE, limits = c("TRUE", "FALSE")) +
    scale_colour_manual(values = cls_col, name = "BAF class", drop = FALSE) +
    scale_x_continuous(breaks = 1:5, labels = lab_x, limits = c(0.5, 5.5), expand = c(0, 0)) +
    scale_y_continuous(limits = ylim) +
    labs(x = NULL, y = if (YVAR == "z") "background z of ρ" else "Spearman ρ", title = titles[[s]]) +
    theme_bw(base_size = 8) + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                                    plot.title = element_text(size = 8, face = "bold"))
}
p <- (panel("TCGA") / panel("DepMap") / panel("Organoid")) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = sprintf(paste0(if (YVAR == "z") "y = (ρ − mean) / SD of the TF's ρ with all genes in the same dataset (comparable across pairs and datasets).\n" else "",
                                           "Each point is one of 28 BAF members. Labels: top %d significant positive members per column (q < 0.05, BH within pair).\n",
                                           "Prostate DepMap n = 10 (underpowered). Association only; not binding or recruitment."), TOP_N),
                  theme = theme(plot.caption = element_text(size = 6, hjust = 0)))
fn <- paste0(if (YVAR == "z") "TFBAF_manhattan_z" else "TFBAF_manhattan_rho", if (INTERNAL) "_internal3D" else "")
ggsave(file.path(fig, paste0(fn, ".pdf")), p, width = 7.5, height = 9, device = cairo_pdf)
ggsave(file.path(fig, paste0(fn, ".png")), p, width = 7.5, height = 9, dpi = 200)
print(bgs)
print(long[label != "", .(source, pair, BAF, rho = round(rho, 2), z = round(z, 1), q = signif(q, 2))][order(source, match(pair, pair_order), -z)])
