# Shared helpers for Step 1 scripts. Source from the repo root:
#   source("scripts/01_tcga_tf_de/common.R")
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

step1_dir <- "scripts/01_tcga_tf_de"
cfg <- yaml::read_yaml(file.path(step1_dir, "config.yaml"))

data_root <- Sys.getenv("TF_CANCER_DATA_ROOT",
                        unset = file.path(Sys.getenv("BIOWORKSPACE_ROOT"), "DATA"))
if (!dir.exists(data_root)) stop("TF_CANCER_DATA_ROOT not found: ", data_root)

dpath <- function(rel) file.path(data_root, rel)
in_path <- function(key) dpath(cfg$inputs[[key]])
processed_dir <- dpath(cfg$outputs$processed_dir)
result_dir <- dpath(cfg$outputs$result_dir)
repo_results <- cfg$outputs$repo_results
for (d in c(processed_dir, result_dir, repo_results, file.path(repo_results, "qc")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

check_md5 <- function(keys) {
  for (k in keys) {
    got <- unname(tools::md5sum(in_path(k)))
    if (!identical(got, cfg$md5[[k]])) stop(sprintf("md5 mismatch for %s: %s", k, got))
  }
  invisible(TRUE)
}

bare_id <- function(x) sub("\\..*$", "", x)

fwrite_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")

log_msg <- function(...) cat(format(Sys.time(), "%F %T"), "|", sprintf(...), "\n")

# Cohort -> organ system lookup (LAML excluded from recurrence by design).
cohort_meta <- rbindlist(lapply(names(cfg$cohorts), function(k) {
  x <- cfg$cohorts[[k]]
  data.table(cohort = k, organ_system = x$organ_system,
             gtex_tissues = paste(unlist(x$gtex), collapse = "; "),
             c2_quality = x$c2_quality, composition_risk = x$composition_risk,
             hc_eligible = isTRUE(x$hc_eligible))
}))

write_run_record <- function(dir, script) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA)
  writeLines(commit, file.path(dir, "code_commit.txt"))
  file.copy(file.path(step1_dir, "config.yaml"), file.path(dir, "config.yaml"), overwrite = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(dir, sprintf("sessionInfo_%s.txt", script)))
  writeLines(paste(commandArgs(), collapse = " "), file.path(dir, sprintf("command_%s.txt", script)))
}
