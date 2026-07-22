#!/bin/bash
#SBATCH -J dwls_decon
#SBATCH -A r01604
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH -o /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration/results/Script_results/dwls_%j.out
#SBATCH -e /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration/results/Script_results/dwls_%j.err

set -eo pipefail

# Some module and Conda initialization scripts expect PS1 to exist.
export PS1=""

PROJECT_DIR="/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration"
SCRIPT_RESULTS_DIR="${PROJECT_DIR}/results/Script_results"
ANALYSIS_RESULTS_DIR="${PROJECT_DIR}/results/deconvolution_comparison"

mkdir -p "${SCRIPT_RESULTS_DIR}"
mkdir -p "${ANALYSIS_RESULTS_DIR}"

cd "${PROJECT_DIR}"

echo "============================================================"
echo "DWLS batch job"
echo "============================================================"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Working directory: $(pwd)"
echo "Script results directory: ${SCRIPT_RESULTS_DIR}"
echo "Analysis results directory: ${ANALYSIS_RESULTS_DIR}"
echo "Requested CPUs: ${SLURM_CPUS_PER_TASK:-1}"
echo "Start time: $(date)"
echo

# ------------------------------------------------------------
# Load software environment
# ------------------------------------------------------------

module purge
module load gnu/12.2.0
module load sqlite/3.35.5
module load rstudio
module swap r/4.3.1
module load conda

source /N/soft/rhel8/conda/26.3.2/etc/profile.d/conda.sh

conda activate \
  /N/u/echimal/Quartz/.conda/envs/integration_env

# ------------------------------------------------------------
# Prevent nested BLAS threading
# ------------------------------------------------------------

# The R adapter parallelizes across ROIs. Keep each worker single-threaded
# internally so 8 workers do not each spawn additional BLAS threads.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Make the requested output directory available to the R script.
export DWLS_OUTPUT_DIR="${ANALYSIS_RESULTS_DIR}"

# ------------------------------------------------------------
# Environment checks
# ------------------------------------------------------------

echo "R executable: $(command -v R || true)"
echo "Rscript executable: $(command -v Rscript || true)"
echo "Conda environment: ${CONDA_PREFIX:-not activated}"
echo

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not available after module setup." >&2
  echo "Loaded modules:" >&2
  module list 2>&1
  exit 1
fi

Rscript --version

echo
echo "Checking required R packages..."

Rscript - <<'RS'
required_packages <- c(
  "DWLS",
  "MAST",
  "Matrix"
)

status <- data.frame(
  package = required_packages,
  installed = vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  ),
  stringsAsFactors = FALSE
)

print(status, row.names = FALSE)

if (!all(status$installed)) {
  missing <- status$package[!status$installed]

  stop(
    "Missing required R package(s): ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}
RS

echo
echo "Starting DWLS analysis..."
echo

# ------------------------------------------------------------
# Run DWLS
# ------------------------------------------------------------

Rscript \
  "${PROJECT_DIR}/deconvolution_comparison/scripts/run_dwls_standalone.R"

echo
echo "============================================================"
echo "DWLS batch job completed"
echo "End time: $(date)"
echo "Expected output directory: ${ANALYSIS_RESULTS_DIR}"
echo "============================================================"
