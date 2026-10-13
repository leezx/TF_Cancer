#!/usr/bin/env bash
# Run Step 1 on the Mac mini (computingResource.md). Run from the repo root:
#   scripts/01_tcga_tf_de/run_step1_local.sh <prep|validate|smoke|full>
# Stages follow the approved order (design §8 round 2): prep -> validate (gate) -> smoke -> full.
set -Eeuo pipefail
export TF_CANCER_DATA_ROOT="${TF_CANCER_DATA_ROOT:-${BIOWORKSPACE_ROOT:-/Volumes/Stelligen_SSD/Stelligen}/DATA}"
export NSLOTS="${NSLOTS:-4}"
S=scripts/01_tcga_tf_de
mkdir -p logs
run() { echo "== $(date +%T) $*"; Rscript "$@"; }
de() { run "$S/04_run_deseq2.R" "$@"; }

case "${1:-}" in
  prep)
    run "$S/01_build_tf_reference.R"
    run "$S/02_build_manifests.R"
    run "$S/03_extract_cohort_matrices.R" ;;
  validate)
    for k in COADREAD LUAD; do
      de "$k" C1 toil matched; de "$k" C1 recount3 matched
      de "$k" C1 toil production; de "$k" C1 recount3 full
    done
    run "$S/05_validate_recount3.R" ;;
  smoke)
    for c in C1 C2 C3; do de KICH "$c" toil production; done ;;
  full)
    [[ "$(cat results/01_tcga_tf_de/qc/input_validation/recount3_gate_status.txt)" == PASS ]] || { echo "recount3 gate not PASS; stop"; exit 1; }
    tail -n +2 results/01_tcga_tf_de/deseq2_tasks.tsv | while read -r k c; do
      f="$TF_CANCER_DATA_ROOT/bulkRNAseq/UCSC_Toil_TCGA_GTEx/result/2026-10-13_tf_cancer_step1_tcga_de_v0.1/de/toil_production/${k}_${c}_TF.tsv"
      [[ -s "$f" ]] && { echo "skip $k $c (exists)"; continue; }
      de "$k" "$c" toil production
    done
    de COADREAD C1 toil read_only
    run "$S/06_integrate.R"
    run "$S/07_plot_step1.R" ;;
  *) echo "usage: $0 <prep|validate|smoke|full>"; exit 2 ;;
esac
