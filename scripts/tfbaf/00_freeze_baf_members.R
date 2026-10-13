#!/usr/bin/env Rscript
# TF–BAF module design §3.1: freeze the mSWI/SNF (BAF) member list by Ensembl ID
# before any data are analysed. Annotation-only. Membership after Mashtalir et al.
# 2018 (Cell) / Centore et al. 2020 (Trends Genet): cBAF, PBAF, ncBAF (GBAF).
# ACTB is excluded (ubiquitous, common-essential); BCL11A/B excluded (Lambert TFs).
# Usage (repo root): Rscript scripts/tfbaf/00_freeze_baf_members.R
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("TF_CANCER_DATA_ROOT", file.path(Sys.getenv("BIOWORKSPACE_ROOT", "/Volumes/Stelligen_SSD/Stelligen"), "DATA"))
hg <- fread(file.path(root, "1.Databases", "HGNC_gene_id_mapping", "processed", "v0.1", "hgnc_symbol_ensembl_map.tsv"))
pm <- fread(file.path(root, "bulkRNAseq", "UCSC_Toil_TCGA_GTEx", "raw", "gencode.v23.annotation.gene.probemap"))[!grepl("_PAR_Y$", id)]
pm[, ensembl_id := sub("\\..*$", "", id)]

baf <- fread(text = "
symbol,class,cBAF,PBAF,ncBAF
SMARCA4,ATPase,1,1,1
SMARCA2,ATPase,1,1,1
SMARCC1,core,1,1,1
SMARCC2,core,1,1,0
SMARCD1,core,1,1,1
SMARCD2,core,1,1,0
SMARCD3,core,1,1,0
ACTL6A,core,1,1,1
ACTL6B,core,1,1,1
BCL7A,core,1,1,1
BCL7B,core,1,1,1
BCL7C,core,1,1,1
SMARCB1,core,1,1,0
SMARCE1,core,1,1,0
SS18,core,1,0,1
SS18L1,core,1,0,1
ARID1A,cBAF-specific,1,0,0
ARID1B,cBAF-specific,1,0,0
DPF1,cBAF-specific,1,0,0
DPF2,cBAF-specific,1,0,0
DPF3,cBAF-specific,1,0,0
PBRM1,PBAF-specific,0,1,0
ARID2,PBAF-specific,0,1,0
BRD7,PBAF-specific,0,1,0
PHF10,PBAF-specific,0,1,0
BRD9,ncBAF-specific,0,0,1
BICRA,ncBAF-specific,0,0,1
BICRAL,ncBAF-specific,0,0,1")
baf <- merge(baf, hg[, .(symbol, ensembl_id)], by = "symbol", all.x = TRUE, sort = FALSE)
baf <- merge(baf, pm[, .(ensembl_id, gencode_v23_id = id, v23_symbol = gene)], by = "ensembl_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(baf) == 28, !anyNA(baf$ensembl_id), !anyNA(baf$gencode_v23_id), !anyDuplicated(baf$symbol))
tf <- fread("results/01_tcga_tf_de/TF_reference.tsv")
# ARID2, DPF1, DPF3 are also in the Lambert list; recorded, not excluded.
baf[, in_lambert_TF := ensembl_id %in% tf$ensembl_id]
setcolorder(baf, c("symbol", "v23_symbol", "ensembl_id", "gencode_v23_id", "class", "cBAF", "PBAF", "ncBAF", "in_lambert_TF"))
fwrite(baf[order(match(class, c("ATPase", "core", "cBAF-specific", "PBAF-specific", "ncBAF-specific")))],
       "results/tfbaf/BAF_members.tsv", sep = "\t")
print(baf[symbol != v23_symbol, .(symbol, v23_symbol)])
cat(nrow(baf), "BAF members frozen\n")
