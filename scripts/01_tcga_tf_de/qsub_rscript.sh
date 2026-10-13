#!/usr/bin/env bash
# SGE wrapper for Step 1 R scripts (Argos). Run from the repo root.
#   Single script:  qsub -N s1_extract -pe smp 8 scripts/01_tcga_tf_de/qsub_rscript.sh scripts/01_tcga_tf_de/03_extract_cohort_matrices.R
#   DESeq2 array:   TASK_FILE=<tsv: cohort contrast source set> qsub -t 1-N -pe smp 8 scripts/01_tcga_tf_de/qsub_rscript.sh scripts/01_tcga_tf_de/04_run_deseq2.R
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/
set -Eeuo pipefail
export TF_CANCER_DATA_ROOT="${TF_CANCER_DATA_ROOT:-/home/zz950/DATA}"
export R_LIBS_USER="${R_LIBS_USER:-/home/zz950/softwares/R_lib_4}"
export TMPDIR="${TMPDIR_OVERRIDE:-/home/zz950/tmp}"
RSCRIPT="${RSCRIPT:-/home/zz950/softwares/miniconda3/envs/r4p3/bin/Rscript}"
script="$1"; shift
args=("$@")
if [[ -n "${TASK_FILE:-}" && "${SGE_TASK_ID:-undefined}" != "undefined" ]]; then
  read -r -a args <<< "$(sed -n "${SGE_TASK_ID}p" "$TASK_FILE")"
fi
echo "host=$(hostname) slots=${NSLOTS:-1} task=${SGE_TASK_ID:-none} script=$script args=${args[*]:-}"
"$RSCRIPT" "$script" ${args[@]+"${args[@]}"}
