# recount3 loader for the validation runs (design §3.4.0). Sourced by
# 04_run_deseq2.R; also writes the aliquot-selection audit.
load_recount3 <- function(k) {
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
  smp[, kept := !duplicated(smp[, .(patient, group)])]
  # Audit of retained/excluded aliquots (review round 3, N1).
  qd <- file.path(repo_results, "qc", "input_validation"); dir.create(qd, recursive = TRUE, showWarnings = FALSE)
  fwrite_tsv(smp[, .(cohort = k, recount3_column = col, sample, patient, project, group, library_size, kept)],
             file.path(qd, sprintf("recount3_aliquot_selection_%s.tsv", k)))
  smp <- smp[kept == TRUE]
  log_msg("recount3 %s: %d columns, %d after one-per-patient-per-group", k, n0, nrow(smp))
  m <- m[, smp$col, drop = FALSE]; colnames(m) <- smp$sample
  list(counts = m, samples = smp[, .(sample, patient, project, group, source = "TCGA")])
}

