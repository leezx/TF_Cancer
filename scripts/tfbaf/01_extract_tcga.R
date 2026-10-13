#!/usr/bin/env Rscript
# TF–BAF design §2/§3.2: extract log2(TPM + 1) for the 5 TFs + 28 BAF members
# for the TCGA tumors of the 5 pairs, straight from the Toil TPM matrix
# (row-filtered while streaming; the TPM file is log2(TPM + 0.001)).
# COADREAD / LAML / PRAD / BRCA tumors come from the Step 1 sample manifest;
# SKCM (not in Step 1) comes from the Toil phenotype file: Primary Tumor (main)
# and Metastatic (sensitivity), one sample per patient within each set.
# Usage (repo root): Rscript scripts/tfbaf/01_extract_tcga.R
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
raw <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx", "raw")
out <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx", "result", "2026-10-13_tf_cancer_tfbaf_v0.1", "tcga")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pairs <- data.table(pair = c("CRC", "AML", "Prostate", "Breast", "Melanoma"),
                    TF = c("SOX9", "SPI1", "FOXA1", "GATA3", "SOX10"),
                    cohort = c("COADREAD", "LAML", "PRAD", "BRCA", "SKCM"))

tf <- fread("results/01_tcga_tf_de/TF_reference.tsv")
baf <- fread("results/tfbaf/BAF_members.tsv")
genes <- rbind(data.table(symbol = pairs$TF, gencode_v23_id = tf$gencode_v23_id[match(pairs$TF, tf$hgnc_symbol)]),
               baf[, .(symbol, gencode_v23_id)])
stopifnot(!anyNA(genes$gencode_v23_id), !anyDuplicated(genes$gencode_v23_id))

# ---- samples ----
sm <- fread("results/01_tcga_tf_de/sample_manifest.tsv")[group == "tumor" & cohort %in% pairs$cohort, .(sample, cohort, patient)]
sm[, set := "main"]
ph <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(raw, "TcgaTargetGTEX_phenotype.txt.gz"))))
sk <- ph[`_study` == "TCGA" & detailed_category == "Skin Cutaneous Melanoma" & `_sample_type` %in% c("Primary Tumor", "Metastatic"),
         .(sample, cohort = "SKCM", patient = substr(sample, 1, 12),
           set = fifelse(`_sample_type` == "Primary Tumor", "main", "skcm_metastatic"))]
setorder(sk, set, patient, sample)
n_dup <- sk[, sum(duplicated(patient)), by = set]
sk <- sk[!duplicated(sk[, .(set, patient)])]
s <- rbind(sm, sk)
fwrite(merge(n_dup, sk[, .N, by = set], by = "set"), file.path(out, "skcm_samples_dedup.tsv"), sep = "\t")

# ---- TPM rows (33 genes) for the selected samples ----
idf <- tempfile(); writeLines(genes$gencode_v23_id, idf)
cmd <- sprintf("gzip -dc %s | awk -F'\\t' 'NR==FNR{k[$1];next} FNR==1 || ($1 in k)' %s -",
               shQuote(file.path(raw, "TcgaTargetGtex_rsem_gene_tpm.gz")), shQuote(idf))
hdr <- names(fread(cmd = sprintf("gzip -dc %s | head -1", shQuote(file.path(raw, "TcgaTargetGtex_rsem_gene_tpm.gz")))))
miss <- setdiff(s$sample, hdr)
if (length(miss)) stop("samples missing from TPM matrix: ", paste(head(miss), collapse = ","))
x <- fread(cmd = cmd, select = c("sample", unique(s$sample)))
stopifnot(nrow(x) == nrow(genes))
m <- as.matrix(x[, -1]); rownames(m) <- genes$symbol[match(x$sample, genes$gencode_v23_id)]
m <- log2(pmax(2^m - 0.001, 0) + 1)

long <- data.table(gene = rep(rownames(m), ncol(m)), sample = rep(colnames(m), each = nrow(m)), log2TPM1 = as.vector(m))
w <- dcast(long, sample ~ gene, value.var = "log2TPM1")
res <- merge(s, w, by = "sample", allow.cartesian = TRUE)
res <- merge(res, pairs, by = "cohort")
fwrite(res, file.path(out, "tcga_tf_baf_log2TPM1.tsv"), sep = "\t")
print(res[, .N, by = .(pair, cohort, set)])
