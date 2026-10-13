#!/usr/bin/env Rscript
# Design §3.2–3.3: sample + cohort manifests, contrast task list, and the
# pre-registered sample-count gate (§5.1). Stops on any gate mismatch.
# Usage (repo root): Rscript scripts/01_tcga_tf_de/02_build_manifests.R
source("scripts/01_tcga_tf_de/common.R")
check_md5(c("pheno"))
th <- cfg$thresholds

header_of <- function(f) strsplit(readLines(gzcon(file(f, "rb")), n = 1), "\t")[[1]][-1]
count_cols <- header_of(in_path("counts"))
tpm_cols <- header_of(in_path("tpm"))
in_both <- intersect(count_cols, tpm_cols)

ph <- fread(in_path("pheno"), encoding = "Latin-1")
setnames(ph, c("sample", "detailed_category", "primary_disease", "primary_site",
               "sample_type", "gender", "study"))
ph <- ph[sample %in% in_both]

# TCGA tumor / normal (study filter first: TARGET AML shares the AML label).
rows <- list()
for (k in names(cfg$cohorts)) {
  pj <- cfg$cohorts[[k]]$projects
  for (p in names(pj)) {
    x <- ph[study == "TCGA" & tolower(trimws(detailed_category)) == pj[[p]]]
    x <- x[sample_type %in% c(cfg$tumor_sample_types, cfg$normal_sample_type)]
    if (nrow(x)) rows[[length(rows) + 1]] <- x[, .(sample, cohort = k, project = p,
      group = fifelse(sample_type == cfg$normal_sample_type, "normal", "tumor"),
      source = "TCGA", tissue = detailed_category, sample_type)]
  }
  gt <- unlist(cfg$cohorts[[k]]$gtex)
  if (length(gt)) {
    x <- ph[study == "GTEX" & detailed_category %in% gt]
    rows[[length(rows) + 1]] <- x[, .(sample, cohort = k, project = "GTEX", group = "gtex",
                                      source = "GTEX", tissue = detailed_category, sample_type)]
  }
}
sm <- rbindlist(rows)
sm[, patient := fifelse(source == "TCGA", substr(sample, 1, 12), sample)]

# One sample per patient per group; a duplicate stops the run (design §3.2 step 5).
dup <- sm[source == "TCGA", .N, by = .(cohort, group, patient)][N > 1]
if (nrow(dup)) { print(dup); stop("Duplicate TCGA patients within a group: investigate before proceeding") }

paired_pts <- sm[source == "TCGA", .(paired = all(c("tumor", "normal") %in% group)), by = .(cohort, patient)][paired == TRUE]
sm[, in_pair := source == "TCGA" & paste(cohort, patient) %in% paired_pts[, paste(cohort, patient)]]

cm <- sm[, .(n_tumor = sum(group == "tumor"), n_normal = sum(group == "normal"),
             n_pairs = sum(in_pair & group == "tumor"), n_gtex = sum(group == "gtex")), by = cohort]
cm <- merge(cohort_meta, cm, by = "cohort", all.x = TRUE, sort = FALSE)
for (v in c("n_tumor", "n_normal", "n_pairs", "n_gtex")) cm[is.na(get(v)), (v) := 0L]

# §5.1 gate
gate <- rbindlist(lapply(names(cfg$expected_n), function(k) {
  e <- cfg$expected_n[[k]]; r <- cm[cohort == k]
  data.table(cohort = k, metric = c("tumor", "normal", "pairs", "gtex"),
             expected = c(e$tumor, e$normal, e$pairs, e$gtex),
             realized = c(r$n_tumor, r$n_normal, r$n_pairs, r$n_gtex))
}))
gate[, pass := expected == realized]
fwrite_tsv(gate, file.path(repo_results, "qc", "sample_count_gate.tsv"))
if (!all(gate$pass)) { print(gate[pass == FALSE]); stop("Sample-count gate FAILED (design §5.1)") }
log_msg("Sample-count gate passed: %d checks", nrow(gate))

cm[, `:=`(
  c1_eligible = n_normal >= th$min_normal_c1,
  c2_eligible = n_gtex > 0,
  c3_eligible = n_pairs >= th$min_pairs_c3)]
cm[, normal_reference_quality := fifelse(c1_eligible, "TCGA_matched", paste0("GTEx_", c2_quality))]
cm[, lineage_match_status := fifelse(cohort == "LAML", "no_valid_normal",
                              fifelse(c1_eligible, "same_patient_population", "cross_source_proxy"))]
cm[, eligible_for_high_confidence_recurrence := hc_eligible & c1_eligible]
cm[, DE_status := fifelse(cohort == "LAML", "no_valid_normal", "run")]
fwrite_tsv(cm, file.path(repo_results, "cohort_manifest.tsv"))
fwrite_tsv(sm, file.path(repo_results, "sample_manifest.tsv"))

tasks <- rbindlist(list(
  cm[c1_eligible == TRUE, .(cohort, contrast = "C1")],
  cm[c2_eligible == TRUE, .(cohort, contrast = "C2")],
  cm[c3_eligible == TRUE, .(cohort, contrast = "C3")]))
setorder(tasks, cohort, contrast)
fwrite_tsv(tasks, file.path(repo_results, "deseq2_tasks.tsv"))
log_msg("Tasks: C1=%d C2=%d C3=%d", sum(tasks$contrast == "C1"), sum(tasks$contrast == "C2"), sum(tasks$contrast == "C3"))
write_run_record(file.path(result_dir, "report", "run_records"), "02_build_manifests")
