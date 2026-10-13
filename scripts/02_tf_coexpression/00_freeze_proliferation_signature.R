#!/usr/bin/env Rscript
# Step 2 design §3.4 (review round 4, B3): freeze the proliferation signature by
# GENCODE v23 gene ID before any expression is analysed. Annotation-only.
# HALLMARK_G2M_CHECKPOINT ∪ HALLMARK_E2F_TARGETS, minus Lambert TFs (matched by
# Ensembl ID, not symbol) and minus MKI67 (held-out check).
# Usage (repo root): Rscript scripts/02_tf_coexpression/00_freeze_proliferation_signature.R
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
gmt <- file.path(root, "genesets", "h.all.v2026.1.Hs.symbols.gmt")
pm  <- file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx", "raw", "gencode.v23.annotation.gene.probemap")
hg  <- file.path(root, "1.Databases", "HGNC_gene_id_mapping", "processed", "v0.1", "hgnc_symbol_ensembl_map.tsv")
out <- "results/02_tf_coexpression"

g <- strsplit(readLines(gmt), "\t"); names(g) <- vapply(g, `[`, "", 1)
sets <- c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS")
sig <- unique(rbindlist(lapply(sets, function(s) data.table(set = s, symbol = g[[s]][-(1:2)]))), by = c("set", "symbol"))
sig <- sig[, .(sets = paste(sort(set), collapse = ";")), by = symbol]

# Symbol -> Ensembl via HGNC approved symbols (MSigDB 2026 uses current names;
# GENCODE v23 has older ones, e.g. H2AFX for H2AX), then Ensembl -> v23 ID.
# Fallback: v23 probemap symbol when HGNC has no Ensembl ID.
p <- fread(pm)[!grepl("_PAR_Y$", id), .(gencode_v23_id = id, v23_symbol = gene, chrom,
                                         ensembl_id = sub("\\..*$", "", id))]
h <- fread(hg)[ensembl_id != "", .(symbol, ensembl_id)]
m <- merge(sig, h, by = "symbol", all.x = TRUE)
m[, map_source := fifelse(is.na(ensembl_id), NA_character_, "HGNC")]
fb <- m[is.na(ensembl_id), symbol]
m[symbol %in% fb, `:=`(ensembl_id = p$ensembl_id[match(symbol, p$v23_symbol)], map_source = "v23_symbol")]
m <- merge(m, p, by = "ensembl_id", all.x = TRUE)
m[is.na(gencode_v23_id), map_source := NA_character_]
m[, n_ids := sum(!is.na(gencode_v23_id)), by = symbol]

tf <- fread("results/01_tcga_tf_de/TF_reference.tsv")
m[, status := fcase(is.na(gencode_v23_id), "unmapped",
                    n_ids > 1, "ambiguous_symbol",
                    ensembl_id %in% tf$ensembl_id, "excluded_lambert_TF",
                    symbol == "MKI67", "excluded_heldout_MKI67",
                    default = "in_signature")]
setorder(m, status, symbol)
fwrite(m[, .(symbol, v23_symbol, sets, gencode_v23_id, ensembl_id, chrom, map_source, status)], file.path(out, "proliferation_signature_genes.tsv"), sep = "\t")
print(m[, .N, by = status])
