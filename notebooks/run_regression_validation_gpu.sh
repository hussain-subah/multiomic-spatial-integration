#!/bin/bash
#SBATCH -J regression_val
#SBATCH -p gpu
#SBATCH -A r01604
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o regression_validation_%j.out
#SBATCH -e regression_validation_%j.err

set -euo pipefail

echo "====================================="
echo "Job started on $(hostname)"
date
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "====================================="

cd /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration

module purge
module load gnu/12.2.0
module load sqlite/3.35.5
module load python/gpu/3.10.10

source /N/soft/rhel8/conda/26.3.2/etc/profile.d/conda.sh
conda activate /N/u/echimal/Quartz/.conda/envs/integration_env

export PYTHONPATH="$(pwd)"

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

if torch.cuda.is_available():
    print("GPU:", torch.cuda.get_device_name(0))
    print("CUDA device count:", torch.cuda.device_count())
else:
    raise RuntimeError("GPU was requested, but CUDA is not available.")
PY

echo "Backing up current inferred signature files..."

STAMP=$(date +%Y%m%d_%H%M%S)
BACKUP_DIR="/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/Regression-model/archive_${STAMP}"

mkdir -p "${BACKUP_DIR}"

for f in \
  /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/Regression-model/AD+CAA_inferred_signatures.csv \
  /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/Regression-model/Control_inferred_signatures.csv
do
  if [[ -f "${f}" ]]; then
    cp "${f}" "${BACKUP_DIR}/"
  fi
done

echo "Backup directory: ${BACKUP_DIR}"

echo "Starting regression training and signature-order validation..."

python notebooks/02_regression_signatures.py

echo "Checking expected validation outputs..."

test -f /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/Regression-model/AD+CAA_inferred_signatures.csv
test -f /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/Regression-model/Control_inferred_signatures.csv
test -f /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/Regression-model/signature_same_position_validation.csv
test -f /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/Regression-model/signature_cross_condition_best_matches_corrected.csv

echo "All expected output files were found."

echo "====================================="
echo "Job finished successfully"
date
echo "====================================="
