#!/usr/bin/env Rscript
# Design §3.4 step 0 + review round 2 §8: independent-pipeline validation gate.
# Compares Toil vs recount3 TF results for COADREAD and LUAD C1:
#   matched = same samples + shared gene universe (primary; gated)
#   full    = Toil production vs recount3 all eligible samples (secondary; reported)
# Thresholds are preset in config.yaml and must not be changed after seeing results.
# Usage (repo root): Rscript scripts/01_tcga_tf_de/05_validate_recount3.R
source("scripts/01_tcga_tf_de/common.R")
suppressPackageStartupMessages(library(SummarizedExperiment))
th <- cfg$thresholds
qc_dir <- file.path(repo_results, "qc", "input_validation")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(run, k) fread(file.path(result_dir, "de", run, sprintf("%s_C1_TF.tsv", k)))

compare <- function(k, a_run, b_run, label) {
  a <- rd(a_run, k); b <- rd(b_run, k)
  x <- merge(a[, .(ensembl_id, TF, lfc_a = log2FC_shrunk, tested_a = tested, tier_a = tier)],
             b[, .(ensembl_id, lfc_b = log2FC_shrunk, tested_b = tested, tier_b = tier)], by = "ensembl_id")
  both <- x[tested_a & tested_b]
  big <- both[abs(lfc_a) > th$recount3_direction_lfc | abs(lfc_b) > th$recount3_direction_lfc]
  t1 <- x[tier_a %in% "Tier1"]
  t1_eval <- t1[tested_b == TRUE]
  up_a <- x[!is.na(tier_a), ensembl_id]; up_b <- x[!is.na(tier_b), ensembl_id]
  disc <- t1[!(tier_b %in% c("Tier1", "Tier2")), .(cohort = k, comparison = label, TF, lfc_toil = lfc_a,
                                                   lfc_recount3 = lfc_b, recount3_tested = tested_b, recount3_tier = tier_b)]
  m <- data.table(cohort = k, comparison = label,
    n_TF_tested_both = nrow(both),
    spearman_rho = suppressWarnings(cor(both$lfc_a, both$lfc_b, method = "spearman")),
    n_big_effect = nrow(big),
    direction_agreement = mean(sign(big$lfc_a) == sign(big$lfc_b)),
    n_toil_tier1 = nrow(t1), n_toil_tier1_evaluable = nrow(t1_eval),
    tier1_recovery_evaluable = mean(t1_eval$tier_b %in% c("Tier1", "Tier2")),
    tier1_recovery_all = mean(t1$tier_b %in% c("Tier1", "Tier2")),
    n_up_toil = length(up_a), n_up_recount3 = length(up_b),
    up_jaccard = length(intersect(up_a, up_b)) / max(1, length(union(up_a, up_b))))
  m[, pass_direction := direction_agreement >= th$recount3_direction_min]
  m[, pass_tier1 := tier1_recovery_evaluable >= th$recount3_tier1_recovery_min]
  list(metrics = m, discordant = disc)
}

res <- list()
for (k in names(cfg$inputs$recount3)) {
  res[[length(res) + 1]] <- compare(k, "toil_matched", "recount3_matched", "matched")
  res[[length(res) + 1]] <- compare(k, "toil_production", "recount3_full", "full")
}
met <- rbindlist(lapply(res, `[[`, "metrics"))
met[, gated := comparison == "matched"]
disc <- rbindlist(lapply(res, `[[`, "discordant"))
fwrite_tsv(met, file.path(qc_dir, "recount3_concordance.tsv"))
fwrite_tsv(disc, file.path(qc_dir, "recount3_discordant_tier1_TFs.tsv"))
print(met)

# Q1 audit: recount3 eligible samples absent from Toil (COADREAD, LUAD).
sm <- fread(file.path(repo_results, "sample_manifest.tsv"))
audit <- rbindlist(lapply(names(cfg$inputs$recount3), function(k) {
  f <- dpath(cfg$inputs$recount3[[k]])
  cd <- rbindlist(lapply(Sys.glob(file.path(dirname(f), "coldata_tcga_*.csv")), fread), fill = TRUE)
  cd <- cd[tcga.cgc_sample_sample_type %in% c("Primary Tumor", cfg$normal_sample_type)]
  stage_col <- grep("tumor_stage|pathologic_stage", names(cd), value = TRUE)[1]
  cd[, `:=`(sample = substr(tcga.tcga_barcode, 1, 15), tss = substr(tcga.tcga_barcode, 6, 7),
            stage = if (is.na(stage_col)) NA_character_ else as.character(get(stage_col)))]
  cd <- unique(cd, by = "sample")
  cd[, in_toil := sample %in% sm[cohort == k, sample]]
  cd[, .(cohort = k, project = study, sample_type = tcga.cgc_sample_sample_type, tss, stage, in_toil)]
}))
fwrite_tsv(audit[, .N, by = .(cohort, project, sample_type, in_toil)], file.path(qc_dir, "toil_absent_audit_counts.tsv"))
fwrite_tsv(audit[, .N, by = .(cohort, project, in_toil, stage)][order(cohort, project, stage)], file.path(qc_dir, "toil_absent_audit_stage.tsv"))
tss <- audit[sample_type == "Primary Tumor", .N, by = .(cohort, tss, in_toil)]
tss <- dcast(tss, cohort + tss ~ in_toil, value.var = "N", fill = 0)
setnames(tss, c("FALSE", "TRUE"), c("n_absent_from_toil", "n_in_toil"), skip_absent = TRUE)
fwrite_tsv(tss[order(cohort, -n_absent_from_toil)], file.path(qc_dir, "toil_absent_audit_tss.tsv"))
print(audit[, .N, by = .(cohort, project, sample_type, in_toil)])

gate_pass <- all(met[gated == TRUE, pass_direction & pass_tier1])
writeLines(if (gate_pass) "PASS" else "FAIL", file.path(qc_dir, "recount3_gate_status.txt"))
log_msg("recount3 validation gate (matched comparisons): %s", if (gate_pass) "PASS" else "FAIL")
write_run_record(file.path(result_dir, "report", "run_records"), "05_validate_recount3")
