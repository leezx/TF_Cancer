#!/usr/bin/env Rscript
# Design §3.3–3.5: one DESeq2 contrast for one cohort.
# Usage (repo root):
#   Rscript scripts/01_tcga_tf_de/04_run_deseq2.R <cohort> <C1|C2|C3> [toil|recount3] [production|matched|full|read_only]
# production = Toil, design sample set (default)
# matched    = Toil or recount3 restricted to the same TCGA samples and the shared
#              bare-Ensembl gene universe (validation, primary; C1 only)
# full       = recount3 with all of its eligible samples (validation, secondary; C1 only)
# read_only  = Toil COADREAD C1 restricted to READ (design §5.7 sensitivity)
source("scripts/01_tcga_tf_de/common.R")
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(BiocParallel) })

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("Usage: 04_run_deseq2.R <cohort> <contrast> [toil|recount3] [production|matched|full]")
k <- a[[1]]; ctr <- a[[2]]
src <- if (length(a) >= 3) a[[3]] else "toil"
set <- if (length(a) >= 4) a[[4]] else "production"
stopifnot(ctr %in% c("C1", "C2", "C3"), src %in% c("toil", "recount3"),
          set %in% c("production", "matched", "full", "read_only"),
          set != "read_only" || (k == "COADREAD" && src == "toil"),
          !(src == "recount3" && set %in% c("production", "read_only")), set == "production" || ctr == "C1")
th <- cfg$thresholds
nslots <- as.integer(Sys.getenv("NSLOTS", "1"))
bp <- if (nslots > 1) MulticoreParam(nslots) else SerialParam()
run_id <- sprintf("%s_%s", src, set)
out_dir <- file.path(result_dir, "de", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_msg("cohort=%s contrast=%s source=%s set=%s slots=%d", k, ctr, src, set, nslots)

toil <- readRDS(file.path(processed_dir, sprintf("toil_%s.rds", k)))
tf <- fread(file.path(repo_results, "TF_reference.tsv"))[in_expression_matrix == TRUE]
s <- copy(toil$samples)

# ---- recount3 loader (validation only) ----
load_recount3 <- function() {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  f <- dpath(cfg$inputs$recount3[[k]])
  rse <- readRDS(f)
  # The combined COADREAD+GTEx object kept only columns shared with GTEx; the TCGA
  # sample type and barcode live in the per-project coldata CSVs next to it.
  if (!"tcga.cgc_sample_sample_type" %in% colnames(colData(rse))) {
    extra <- rbindlist(lapply(Sys.glob(file.path(dirname(f), "coldata_tcga_*.csv")), function(x)
      fread(x, select = c("V1", "tcga.cgc_sample_sample_type", "tcga.tcga_barcode"))), use.names = TRUE)
    i <- match(colnames(rse), extra$V1)
    colData(rse)$tcga.cgc_sample_sample_type <- extra$tcga.cgc_sample_sample_type[i]
    colData(rse)$tcga.tcga_barcode <- extra$tcga.tcga_barcode[i]
  }
  cd <- as.data.frame(colData(rse))
  keep <- cd$study %in% names(cfg$cohorts[[k]]$projects) &
    cd$tcga.cgc_sample_sample_type %in% c("Primary Tumor", cfg$normal_sample_type)
  rse <- rse[!grepl("_PAR_Y$", rownames(rse)), keep]
  cd <- as.data.frame(colData(rse))
  m <- assay(rse, "counts")
  storage.mode(m) <- "integer"
  rownames(m) <- bare_id(rownames(m))
  smp <- data.table(col = colnames(m), sample = substr(cd$tcga.tcga_barcode, 1, 15),
                    patient = substr(cd$tcga.tcga_barcode, 1, 12), project = cd$study,
                    group = fifelse(cd$tcga.cgc_sample_sample_type == cfg$normal_sample_type, "normal", "tumor"),
                    library_size = colSums(m))
  # Several aliquots per sample/patient can exist in recount3: keep the largest library
  # (implementation choice for the validation runs; logged).
  setorder(smp, patient, group, -library_size)
  n0 <- nrow(smp)
  smp <- smp[!duplicated(smp[, .(patient, group)])]
  log_msg("recount3 %s: %d columns, %d after one-per-patient-per-group", k, n0, nrow(smp))
  m <- m[, smp$col, drop = FALSE]; colnames(m) <- smp$sample
  list(counts = m, samples = smp[, .(sample, patient, project, group, source = "TCGA")])
}

# ---- choose samples and count matrix ----
if (ctr == "C1") s <- s[group %in% c("tumor", "normal")]
if (ctr == "C2") s <- s[group %in% c("tumor", "gtex")]
if (ctr == "C3") s <- s[group %in% c("tumor", "normal") & in_pair == TRUE]
toil_counts <- toil$counts
rownames(toil_counts) <- bare_id(rownames(toil_counts))

if (set == "read_only") s <- s[project == "READ"]
if (set %in% c("production", "read_only")) {
  counts <- toil_counts[, s$sample, drop = FALSE]
} else {
  r3 <- load_recount3()
  if (set == "matched") {
    shared <- intersect(s$sample, r3$samples$sample)
    genes <- intersect(rownames(toil_counts), rownames(r3$counts))
    s <- s[sample %in% shared]
    counts <- if (src == "toil") toil_counts[genes, s$sample, drop = FALSE] else r3$counts[genes, s$sample, drop = FALSE]
    log_msg("matched: %d samples (of %d Toil), %d shared genes", length(shared), nrow(toil$samples[group %in% c("tumor","normal")]), length(genes))
  } else {
    s <- r3$samples
    counts <- r3$counts[, s$sample, drop = FALSE]
  }
}

s[, condition := factor(fifelse(group == "tumor", "tumor", "ref"), levels = c("ref", "tumor"))]
if (anyDuplicated(s[condition == "tumor", patient]) || anyDuplicated(s[condition == "ref", patient]))
  stop("duplicate patient within a condition")
form <- if (ctr == "C3") ~ patient + condition else if (k == "COADREAD" && ctr == "C1" && set != "read_only") ~ project + condition else ~ condition
s[, patient := factor(patient)]; s[, project := factor(project)]
n_t <- sum(s$condition == "tumor"); n_r <- sum(s$condition == "ref")

keep <- rowSums(counts >= th$min_count) >= min(n_t, n_r)
log_msg("samples tumor=%d ref=%d; genes kept %d / %d", n_t, n_r, sum(keep), nrow(counts))
cd <- as.data.frame(s); rownames(cd) <- s$sample
stopifnot(identical(colnames(counts), rownames(cd)))
dds <- DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], cd, design = form)
dds <- DESeq(dds, parallel = nslots > 1, BPPARAM = bp, quiet = TRUE)
coef <- "condition_tumor_vs_ref"
res <- results(dds, name = coef, alpha = th$padj)
shr <- lfcShrink(dds, coef = coef, type = "apeglm", parallel = nslots > 1, BPPARAM = bp, quiet = TRUE)

all <- data.table(gene = rownames(res), baseMean = res$baseMean, log2FC_MLE = res$log2FoldChange,
                  lfcSE = res$lfcSE, log2FC_shrunk = shr$log2FoldChange, pvalue = res$pvalue, padj = res$padj)
fwrite(all, file.path(out_dir, sprintf("%s_%s_allgenes.tsv.gz", k, ctr)), sep = "\t")

# ---- TF table (all TFs, untested ones flagged) ----
ref_group <- if (ctr == "C2") "gtex" else "normal"
tpm <- toil$tpm_tf
rownames(tpm) <- bare_id(tf$gencode_v23_id[match(rownames(tpm), tf$gencode_v23_id)])
tum_ids <- intersect(s[condition == "tumor", sample], colnames(tpm))
ref_ids <- intersect(s[condition == "ref", sample], colnames(tpm))
# recount3-full samples absent from Toil have no TPM: the TPM gate uses the Toil cohort
# (production C1 samples) so both pipelines are filtered identically.
if (set == "full") {
  tum_ids <- toil$samples[group == "tumor", sample]; ref_ids <- toil$samples[group == ref_group, sample]
}
med <- function(ids) if (length(ids)) apply(tpm[, ids, drop = FALSE], 1, median) else NA_real_
tft <- data.table(ensembl_id = tf$ensembl_id, TF = tf$hgnc_symbol)
tft[, tumor_median_TPM := med(tum_ids)[ensembl_id]]
tft[, ref_median_TPM := med(ref_ids)[ensembl_id]]
tft[, tumor_prevalence := rowMeans(tpm[ensembl_id, tum_ids, drop = FALSE] >= th$prevalence_tpm)]
tft <- merge(tft, all, by.x = "ensembl_id", by.y = "gene", all.x = TRUE, sort = FALSE)
tft[, `:=`(cohort = k, contrast = ctr, run = run_id, n_tumor = n_t, n_ref = n_r,
           ref_source = fifelse(ctr == "C2", "GTEx", "TCGA"), n_tested_genes = sum(!is.na(res$padj)),
           shrinkage_method = "apeglm", tested = !is.na(padj),
           in_gene_universe = ensembl_id %in% rownames(counts))]
tft[, up := tested & padj < th$padj & tumor_median_TPM >= th$tumor_tpm_min &
           tumor_prevalence >= th$tumor_prevalence_min & log2FC_shrunk >= th$lfc_tier2]
tft[, tier := fifelse(up & log2FC_shrunk >= th$lfc_tier1, "Tier1", fifelse(up, "Tier2", NA_character_))]
fwrite_tsv(tft, file.path(out_dir, sprintf("%s_%s_TF.tsv", k, ctr)))
saveRDS(dds, file.path(out_dir, sprintf("%s_%s_dds.rds", k, ctr)))
log_msg("done: TF tested %d, Tier1 %d, Tier2 %d", sum(tft$tested), sum(tft$tier == "Tier1", na.rm = TRUE), sum(tft$tier == "Tier2", na.rm = TRUE))
write_run_record(file.path(out_dir, "run_records"), sprintf("04_%s_%s", k, ctr))
