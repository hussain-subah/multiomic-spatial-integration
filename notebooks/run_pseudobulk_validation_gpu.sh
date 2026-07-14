#!/bin/bash
#SBATCH -J pseudobulk_val
#SBATCH -p gpu
#SBATCH -A r01604
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o pseudobulk_validation_%j.out
#SBATCH -e pseudobulk_validation_%j.err

set -euo pipefail

echo "====================================="
echo "Job started on $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID}"
date
echo "====================================="

PROJECT_DIR="/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration"
OUTPUT_DIR="${PROJECT_DIR}/results/pseudobulk_validation"

cd "${PROJECT_DIR}"

module purge
module load gnu/12.2.0
module load sqlite/3.35.5
module load python/gpu/3.10.10
module list

source /N/soft/rhel8/conda/26.3.2/etc/profile.d/conda.sh
conda activate /N/u/echimal/Quartz/.conda/envs/integration_env

export PYTHONPATH="${PROJECT_DIR}"

export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16

echo "Python:"
which python
python --version

echo "GPU check:"
python - <<'PY'
import torch

print("CUDA available:", torch.cuda.is_available())

if not torch.cuda.is_available():
    raise RuntimeError("GPU requested, but CUDA is unavailable.")

print("GPU:", torch.cuda.get_device_name(0))
print("CUDA devices:", torch.cuda.device_count())
PY

echo "Checking Notebook 6 syntax..."
python -m py_compile notebooks/06_pseudobulk_validation.py

echo "Archiving previous pseudobulk results..."

STAMP=$(date +%Y%m%d_%H%M%S)

if [[ -d "${OUTPUT_DIR}" ]] && \
   find "${OUTPUT_DIR}" -mindepth 1 -print -quit | grep -q .
then
    ARCHIVE_DIR="${PROJECT_DIR}/results/archive/pseudobulk_validation_${STAMP}"
    mkdir -p "${ARCHIVE_DIR}"
    cp -a "${OUTPUT_DIR}/." "${ARCHIVE_DIR}/"
    rm -rf "${OUTPUT_DIR}"
    echo "Archived previous results to: ${ARCHIVE_DIR}"
fi

mkdir -p "${OUTPUT_DIR}"

echo "Starting pseudobulk validation..."
python notebooks/06_pseudobulk_validation.py

echo "Checking expected outputs..."

for condition in "AD+CAA" "Control"; do
    test -f "${OUTPUT_DIR}/${condition}/recovery_metrics.csv"
    test -f "${OUTPUT_DIR}/${condition}/ground_truth_proportions.csv"
    test -f "${OUTPUT_DIR}/${condition}/inferred_proportions.csv"
    test -f "${OUTPUT_DIR}/${condition}/recovery_scatter.pdf"
    test -f "${OUTPUT_DIR}/${condition}/recovery_bias.pdf"
done

echo "All expected pseudobulk validation outputs were found."

echo "====================================="
echo "Job finished successfully"
date
echo "====================================="
