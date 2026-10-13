#!/usr/bin/env Rscript
# TF–BAF v0.2 (docs/TFBAF_V02_DESIGN.md) inputs:
#   §2/§3.1 organoid inventory and CRC unit selection, from the DepMap 3D organoid
#     biobank supplementary tables (zip streamed; nothing written next to it);
#   §3.3 common-essential reference set and per-screen depth;
#   §5 TCGA background: Spearman of each pair's TF against every protein-coding gene
#     with median TPM >= 1, over the v0.1 main samples.
# Usage (repo root): Rscript scripts/tfbaf/10_v02_extract.R
suppressPackageStartupMessages({ library(data.table); library(readxl) })
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
toil <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx")
res_root <- file.path(toil, "result", "2026-10-13_tf_cancer_tfbaf_v0.2")
for (d in file.path(res_root, c("organoid", "background", "report"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
zip <- file.path(root, "1.Databases", "DepMap3D_202608", "01_tumour_derived_organoid_biobank", "MOESM3_Supplementary_Tables.zip")
baf <- fread("results/tfbaf/BAF_members.tsv")
tfs <- c(CRC = "SOX9", AML = "SPI1", Prostate = "FOXA1", Breast = "GATA3", Melanoma = "SOX10")
log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", sprintf(...))
unz <- function(f) sprintf("unzip -p %s %s", shQuote(zip), f)

# ---- organoids: clinical, screens, LFC, depletion ----
xd <- tempfile(); dir.create(xd)
system2("unzip", c("-q", "-o", shQuote(zip), "supplementary_table_2_revision.xlsx", "-d", shQuote(xd)))
cl <- as.data.table(read_excel(file.path(xd, "supplementary_table_2_revision.xlsx"), "clinical_data"))
cl <- cl[, .(sample_ID, individual_ID, primary_tumour_type, icd = `ICD-O-3 Code and Morphology`,
             met = organoid_derived_from_metastasis)]
lfc <- fread(cmd = unz("supplementary_table_6_revision.csv"), sep = " ", select = c("sample_ID", "gene", "library", "LFC"))
dep <- fread(cmd = unz("supplementary_table_5_revision.csv"), sep = " ", select = c("sample_ID", "gene", "library", "is_depleted"))
log_msg("LFC rows %d, depletion rows %d", nrow(lfc), nrow(dep))
lfc[, screen := paste(sample_ID, library, sep = "|")]; dep[, screen := paste(sample_ID, library, sep = "|")]
stopifnot(!anyDuplicated(lfc[, .(screen, gene)]), !anyDuplicated(dep[, .(screen, gene)]))

scr <- lfc[, .(sample_ID = sample_ID[1], library = library[1], n_genes = .N), by = screen]
scr <- merge(scr, cl, by = "sample_ID", all.x = TRUE)
stopifnot(!anyNA(scr$primary_tumour_type))
scr[, malignant := grepl("/3", icd, fixed = TRUE)]
fwrite(scr[order(primary_tumour_type, sample_ID, library)], file.path(res_root, "organoid", "screen_inventory.tsv"), sep = "\t")

# §3.1 unit selection: CRC, malignant -> minLib screen if two -> primary-derived first -> first sample_ID
u <- scr[primary_tumour_type == "Colorectal" & malignant]
u[, lib_rank := fifelse(library == "minLib", 1L, 2L)]
setorder(u, sample_ID, lib_rank); u <- u[!duplicated(sample_ID)]
u[, met_rank := fifelse(met == "No", 1L, 2L)]
setorder(u, individual_ID, met_rank, sample_ID); u <- u[!duplicated(individual_ID)]
fwrite(u[, .(screen, sample_ID, individual_ID, library, met, icd)], file.path(res_root, "organoid", "crc_selected_units.tsv"), sep = "\t")
log_msg("CRC units: %d individuals (minLib %d, v1_1 %d; primary-derived %d)", nrow(u),
        sum(u$library == "minLib"), sum(u$library == "v1_1"), sum(u$met == "No"))

# CRC unit LFC matrix (units x genes), all genes
M <- dcast(lfc[screen %in% u$screen], screen ~ gene, value.var = "LFC")
M <- as.matrix(M, rownames = "screen")[u$screen, , drop = FALSE]
saveRDS(M, file.path(res_root, "organoid", "crc_units_lfc_matrix.rds"))

# §3.3 reference set: is_depleted = 1 in >= 90% of all 178 screens, excluding the 33 TF/BAF genes
n_scr <- uniqueN(dep$screen)
ref <- dep[, .(n_depleted = sum(is_depleted == 1)), by = gene][n_depleted >= 0.9 * n_scr & !gene %in% c(tfs, baf$symbol)]
fwrite(ref[order(gene)], file.path(res_root, "organoid", "depth_reference_genes.tsv"), sep = "\t")
depth <- lfc[gene %in% ref$gene, .(depth = median(LFC), n_ref = .N), by = screen]
fwrite(depth, file.path(res_root, "organoid", "screen_depth.tsv"), sep = "\t")
log_msg("depth reference: %d genes (of %d screens)", nrow(ref), n_scr)

# 33 TF/BAF genes over all screens (positive control, prevalence)
g33 <- c(tfs, baf$symbol)
fwrite(merge(lfc[gene %in% g33, .(screen, gene, LFC)], dep[gene %in% g33, .(screen, gene, is_depleted)], by = c("screen", "gene")),
       file.path(res_root, "organoid", "all_screens_tf_baf.tsv"), sep = "\t")
rm(lfc, dep); gc(verbose = FALSE)

# ---- TCGA background (§5) ----
x1 <- fread(file.path(toil, "result", "2026-10-13_tf_cancer_tfbaf_v0.1", "tcga", "tcga_tf_baf_log2TPM1.tsv"))[set == "main"]
hg <- fread(file.path(root, "1.Databases", "HGNC_gene_id_mapping", "raw", "hgnc_custom_download.tsv"))
pc_ens <- hg[Status == "Approved" & `UniProt ID(supplied by UniProt)` != "" & `Ensembl gene ID` != "", unique(`Ensembl gene ID`)]
tpm_f <- file.path(toil, "raw", "TcgaTargetGtex_rsem_gene_tpm.gz")
X <- fread(cmd = sprintf("gzip -dc %s", shQuote(tpm_f)), select = c("sample", unique(x1$sample)))
ids <- X$sample; X <- as.matrix(X[, -1]); rownames(X) <- ids
X <- log2(pmax(2^X - 0.001, 0) + 1)
pm <- fread(file.path(toil, "raw", "gencode.v23.annotation.gene.probemap"))
sym <- pm$gene[match(ids, pm$id)]
is_pc <- sub("\\..*$", "", ids) %in% pc_ens
log_msg("TPM matrix %d genes x %d samples; protein-coding (HGNC+UniProt) %d", nrow(X), ncol(X), sum(is_pc))
bg <- rbindlist(lapply(names(tfs), function(pp) {
  s <- x1[pair == pp, sample]; Y <- X[, s, drop = FALSE]
  tf_row <- which(sym == tfs[[pp]]); stopifnot(length(tf_row) == 1)
  med_ok <- apply(Y, 1, median) >= log2(2)
  zero_var <- apply(Y, 1, var) == 0
  keep <- (is_pc & med_ok & !zero_var) | sym %in% baf$v23_symbol   # BAF members always scored (v23 names); flagged below
  keep[tf_row] <- FALSE
  R <- t(apply(Y[keep, , drop = FALSE], 1, rank)); rt <- rank(Y[tf_row, ])
  rho <- as.vector(cor(t(R), rt))
  log_msg("TCGA background %s: %d genes scored", pp, sum(keep))
  sk <- sym[keep]; i23 <- match(sk, baf$v23_symbol); sk[!is.na(i23)] <- baf$symbol[i23[!is.na(i23)]]  # v23 -> current BAF names
  data.table(pair = pp, TF = tfs[[pp]], gene_id = ids[keep], symbol = sk, rho = rho,
             in_universe = (is_pc & med_ok & !zero_var)[keep],
             n_samples = length(s), n_pc_genes = sum(is_pc), n_pc_median_ge1 = sum(is_pc & med_ok),
             n_pc_zero_var = sum(is_pc & zero_var), n_missing = sum(!is.finite(Y)))
}))
fwrite(bg, file.path(res_root, "background", "tcga_tf_allgene_rho.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(res_root, "report", "sessionInfo_10_v02_extract.txt"))
writeLines(system("git rev-parse HEAD", intern = TRUE), file.path(res_root, "report", "code_commit_10.txt"))
log_msg("done")
