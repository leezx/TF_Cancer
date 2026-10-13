#!/usr/bin/env Rscript
# Design §3.1: TF universe = Lambert v1.01 `Is TF? == Yes`, matched by bare
# Ensembl ID to the GENCODE v23 probemap used by the Toil matrices.
# Usage (repo root): Rscript scripts/01_tcga_tf_de/01_build_tf_reference.R
source("scripts/01_tcga_tf_de/common.R")
check_md5(c("lambert", "probemap"))

lam <- fread(in_path("lambert"))
setnames(lam, make.names(names(lam)))
tf <- lam[Is.TF. == "Yes", .(ensembl_id = Ensembl.ID, hgnc_symbol = HGNC.symbol,
                             DBD, binding_mode = Binding.mode, motif_status = Motif.status)]
stopifnot(nrow(tf) == 1639, !anyDuplicated(tf$ensembl_id))

pm <- fread(in_path("probemap"))
pm[, bare := bare_id(id)]
stopifnot(!anyDuplicated(pm$bare))

tf <- merge(tf, pm[, .(ensembl_id = bare, gencode_v23_id = id, gencode_v23_symbol = gene)],
            by = "ensembl_id", all.x = TRUE, sort = FALSE)
tf[, in_expression_matrix := !is.na(gencode_v23_id)]
tf[, symbol_mismatch := in_expression_matrix & hgnc_symbol != gencode_v23_symbol]
setorder(tf, hgnc_symbol)

log_msg("TFs: %d; in matrix: %d; missing: %s; symbol mismatches: %d",
        nrow(tf), sum(tf$in_expression_matrix),
        paste(tf[!(in_expression_matrix), hgnc_symbol], collapse = ","), sum(tf$symbol_mismatch))
stopifnot(sum(tf$in_expression_matrix) == 1637)
fwrite_tsv(tf, file.path(repo_results, "TF_reference.tsv"))
