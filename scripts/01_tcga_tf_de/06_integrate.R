#!/usr/bin/env Rscript
# Design §3.6 + §5: integrate production DESeq2 TF tables into the long-format
# atlas, TF x cohort matrices, recurrence metrics, and QC tables.
# Usage (repo root): Rscript scripts/01_tcga_tf_de/06_integrate.R
source("scripts/01_tcga_tf_de/common.R")
th <- cfg$thresholds
qc <- file.path(repo_results, "qc")
cm <- fread(file.path(repo_results, "cohort_manifest.tsv"))
tasks <- fread(file.path(repo_results, "deseq2_tasks.tsv"))
prod_dir <- file.path(result_dir, "de", "toil_production")

f <- file.path(prod_dir, sprintf("%s_%s_TF.tsv", tasks$cohort, tasks$contrast))
missing <- f[!file.exists(f)]
if (length(missing)) stop("Missing DESeq2 outputs:\n", paste(missing, collapse = "\n"))
de <- rbindlist(lapply(f, fread), use.names = TRUE, fill = TRUE)

# LAML: expression only, DE_status = no_valid_normal (design §3.2).
laml <- readRDS(file.path(processed_dir, "toil_LAML.rds"))
tf <- fread(file.path(repo_results, "TF_reference.tsv"))[in_expression_matrix == TRUE]
tp <- laml$tpm_tf[tf$gencode_v23_id, laml$samples[group == "tumor", sample], drop = FALSE]
de <- rbind(de, data.table(ensembl_id = tf$ensembl_id, TF = tf$hgnc_symbol, cohort = "LAML", contrast = "none",
  run = "toil_production", tumor_median_TPM = apply(tp, 1, median),
  tumor_prevalence = rowMeans(tp >= th$prevalence_tpm), n_tumor = ncol(tp), n_ref = 0L,
  tested = FALSE, up = FALSE), fill = TRUE)

de <- merge(de, cm[, .(cohort, organ_system, normal_reference_quality, composition_risk,
                       eligible_for_high_confidence_recurrence, DE_status)], by = "cohort", sort = FALSE)
de[, evidence_class := fcase(contrast == "C1", "high_confidence_discovery",
                             contrast == "C2", "exploratory_cross_source",
                             contrast == "C3", "supporting_sensitivity",
                             default = "expression_only")]
de[, cross_source_only := contrast == "C2"]
de[, primary_for_ranking := contrast == "C1"]
cols <- c("cohort", "organ_system", "contrast", "evidence_class", "cross_source_only", "primary_for_ranking",
          "TF", "ensembl_id", "baseMean", "log2FC_MLE", "lfcSE", "log2FC_shrunk", "shrinkage_method",
          "pvalue", "padj", "n_tested_genes", "tumor_median_TPM", "ref_median_TPM", "tumor_prevalence",
          "n_tumor", "n_ref", "ref_source", "normal_reference_quality", "composition_risk", "tested", "up", "tier", "DE_status")
de <- de[, ..cols]
setorder(de, cohort, contrast, TF)
fwrite_tsv(de, file.path(repo_results, "TCGA_TF_DE_all.tsv"))
fwrite_tsv(de[!is.na(tier)], file.path(repo_results, "TCGA_TF_DE_significant.tsv"))

# ---- matrices (one set per evidence class) ----
wide <- function(dt, v) dcast(dt, TF + ensembl_id ~ cohort, value.var = v)
for (ct in c("C1", "C2")) {
  x <- de[contrast == ct]
  fwrite_tsv(wide(x, "log2FC_shrunk"), file.path(repo_results, sprintf("TF_x_cohort_log2FC_%s.tsv", ct)))
  fwrite_tsv(wide(x, "padj"), file.path(repo_results, sprintf("TF_x_cohort_padj_%s.tsv", ct)))
}
tum <- unique(de[contrast %in% c("C1", "C2", "none")][order(match(contrast, c("C1", "C2", "none")))], by = c("TF", "cohort"))
fwrite_tsv(wide(tum, "tumor_median_TPM"), file.path(repo_results, "TF_x_cohort_tumorTPM.tsv"))

# ---- recurrence (design §3.6, review B4) ----
d <- de[organ_system != "blood" & contrast %in% c("C1", "C2", "C3")]
w <- dcast(d, TF + ensembl_id + cohort + organ_system ~ contrast, value.var = c("up", "tier", "log2FC_shrunk", "tested"))
for (v in grep("^(up|tested)_", names(w), value = TRUE)) w[is.na(get(v)), (v) := FALSE]
w[, has_C1 := cohort %in% cm[c1_eligible == TRUE, cohort]]
w[, has_C3 := cohort %in% cm[c3_eligible == TRUE, cohort]]
org <- w[, .(
  C1_up = any(up_C1), C1_tier1 = any(tier_C1 %in% "Tier1"),
  C1_C3 = any(up_C1 & has_C3 & tested_C3 & log2FC_shrunk_C3 > 0),
  C2_up = any(up_C2), exploratory = any(up_C1 | up_C2),
  discordant = any(has_C1 & tested_C1 & tested_C2 &
                     ((up_C1 & log2FC_shrunk_C2 < 0) | (up_C2 & log2FC_shrunk_C1 < 0))),
  elig_C1 = any(has_C1)), by = .(TF, ensembl_id, organ_system)]
rec <- org[, .(
  n_organ_systems_C1_up = sum(C1_up), n_organ_systems_C1_tier1 = sum(C1_tier1),
  n_organ_systems_C1_C3_supported = sum(C1_C3), n_organ_systems_C2_only_up = sum(C2_up & !C1_up),
  n_organ_systems_exploratory_up = sum(exploratory), n_organ_systems_discordant = sum(discordant),
  n_organ_systems_eligible_C1 = sum(elig_C1), n_organ_systems_eligible_exploratory = .N,
  organ_systems_C1_up = paste(sort(organ_system[C1_up]), collapse = ";"),
  organ_systems_C2_only_up = paste(sort(organ_system[C2_up & !C1_up]), collapse = ";")), by = .(TF, ensembl_id)]
cc <- de[contrast == "C1", .(n_cohorts_C1_tier1 = sum(tier %in% "Tier1"), n_cohorts_C1_up = sum(up, na.rm = TRUE)), by = TF]
rec <- merge(rec, cc, by = "TF", all.x = TRUE)
rec[, fig1B_primary := n_organ_systems_C1_up >= th$fig1b_min_organ_systems]
rec[, fig1B_supplementary := n_organ_systems_exploratory_up >= th$fig1b_min_organ_systems]
# recount3-discordant Toil Tier 1 TFs (review round 3, N3): carry the flag downstream.
disc <- fread(file.path(qc, "input_validation", "recount3_discordant_tier1_TFs.tsv"))
rec[, pipeline_sensitive := TF %in% disc$TF]
setorder(rec, -n_organ_systems_C1_up, -n_organ_systems_C1_tier1, -n_organ_systems_exploratory_up, TF)
fwrite_tsv(rec, file.path(repo_results, "TF_recurrent_upregulation.tsv"))

# ---- QC: contrast concordance (§5.6) ----
conc <- rbindlist(lapply(unique(de[contrast == "C1", cohort]), function(k) rbindlist(lapply(c("C2", "C3"), function(b) {
  x <- merge(de[cohort == k & contrast == "C1", .(TF, a = log2FC_shrunk, ta = tested, ua = up)],
             de[cohort == k & contrast == b, .(TF, b = log2FC_shrunk, tb = tested, ub = up)], by = "TF")
  if (!nrow(x)) return(NULL)
  x <- x[ta & tb]; big <- x[abs(a) >= 1 | abs(b) >= 1]
  data.table(cohort = k, comparison = paste0("C1_vs_", b), n_TF = nrow(x),
             spearman_rho = cor(x$a, x$b, method = "spearman"),
             direction_agreement_absLFCge1 = mean(sign(big$a) == sign(big$b)),
             up_jaccard = sum(x$ua & x$ub) / max(1, sum(x$ua | x$ub)))
}))))
conc[, flag_low_rho := spearman_rho < th$concordance_rho_flag]
fwrite_tsv(conc, file.path(qc, "contrast_concordance.tsv"))

# READ-only sensitivity (§5.7)
ro <- file.path(result_dir, "de", "toil_read_only", "COADREAD_C1_TF.tsv")
if (file.exists(ro)) {
  r <- merge(de[cohort == "COADREAD" & contrast == "C1", .(TF, a = log2FC_shrunk, ta = tested, ua = up)],
             fread(ro)[, .(TF, b = log2FC_shrunk, tb = tested, ub = up)], by = "TF")[ta & tb]
  fwrite_tsv(data.table(comparison = "COADREAD_C1_vs_READonly_C1", n_TF = nrow(r),
                        spearman_rho = cor(r$a, r$b, method = "spearman"),
                        up_jaccard = sum(r$ua & r$ub) / max(1, sum(r$ua | r$ub))),
             file.path(qc, "read_only_sensitivity.tsv"))
}

# ---- QC: controls (§5.4, §5.5) ----
idc <- rbindlist(lapply(names(cfg$identity_controls), function(k)
  data.table(cohort = k, TF = unlist(cfg$identity_controls[[k]]))))
idc <- merge(idc, tum[, .(cohort, TF, tumor_median_TPM)], by = c("cohort", "TF"), all.x = TRUE)
idc[, pass := tumor_median_TPM >= 1]
fwrite_tsv(idc, file.path(qc, "identity_controls.tsv"))
pro <- de[contrast == "C1" & TF %in% unlist(cfg$proliferation_controls),
          .(cohort, TF, log2FC_shrunk, padj, up, tier)]
fwrite_tsv(pro, file.path(qc, "proliferation_controls.tsv"))

# ---- QC: TF coverage and filter sensitivity (§5.2, N2) ----
cov <- de[contrast != "none", .(n_TF_tested = sum(tested), n_TF_untested = sum(!tested),
                                n_untested_with_tumorTPMge1 = sum(!tested & tumor_median_TPM >= 1, na.rm = TRUE),
                                n_tested_genes = max(n_tested_genes, na.rm = TRUE)), by = .(cohort, contrast)]
fwrite_tsv(cov, file.path(qc, "tf_coverage_by_contrast.tsv"))
never <- de[contrast != "none", .(ever = any(tested)), by = TF][ever == FALSE, TF]
writeLines(sort(never), file.path(qc, "tf_never_tested.txt"))

log_msg("integrated %d rows; C1 Tier1/2 calls %d; Fig1B primary %d TFs, supplementary %d TFs; identity controls pass %d/%d",
        nrow(de), de[contrast == "C1" & !is.na(tier), .N], sum(rec$fig1B_primary), sum(rec$fig1B_supplementary),
        sum(idc$pass, na.rm = TRUE), nrow(idc))
write_run_record(file.path(result_dir, "report", "run_records"), "06_integrate")
