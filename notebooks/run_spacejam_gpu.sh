#!/bin/bash
#SBATCH -J spacejam_fit
#SBATCH -p gpu
#SBATCH -A r01604
#SBATCH --gpus-per-node=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o spacejam_%j.out
#SBATCH -e spacejam_%j.err

set -euo pipefail

echo "====================================="
echo "Job started on $(hostname)"
date
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "====================================="

PROJECT_DIR="/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration"
SPACEJAM_DIR="${PROJECT_DIR}/results/spacejam"
REGRESSION_DIR="${PROJECT_DIR}/results/regression_model"

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
    raise RuntimeError(
        "GPU was requested, but CUDA is not available."
    )

print("GPU:", torch.cuda.get_device_name(0))
print("CUDA device count:", torch.cuda.device_count())
PY

echo "Checking Notebook 3 syntax..."
python -m py_compile notebooks/03_spacejam_cell2location.py

echo "Checking corrected regression signatures..."

test -f "${REGRESSION_DIR}/AD+CAA_inferred_signatures.csv"
test -f "${REGRESSION_DIR}/Control_inferred_signatures.csv"

python - <<'PY'
from pathlib import Path
import pandas as pd

project = Path(
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/multiomic-spatial-integration"
)

signature_dir = project / "results" / "regression_model"

ad = pd.read_csv(
    signature_dir / "AD+CAA_inferred_signatures.csv",
    index_col=0
)

ctrl = pd.read_csv(
    signature_dir / "Control_inferred_signatures.csv",
    index_col=0
)

print("AD+CAA signature shape:", ad.shape)
print("Control signature shape:", ctrl.shape)
print("Same factor set:", set(ad.columns) == set(ctrl.columns))
print("Same factor order:", ad.columns.tolist() == ctrl.columns.tolist())

if ad.columns.tolist() != ctrl.columns.tolist():
    raise RuntimeError(
        "AD+CAA and Control signature columns are not in the same order."
    )

if len(ad.columns) != 46:
    raise RuntimeError(
        f"Expected 46 factors, found {len(ad.columns)}."
    )

print("Signature validation passed.")
PY

echo "Backing up existing SpaceJam results..."

STAMP=$(date +%Y%m%d_%H%M%S)

if [[ -d "${SPACEJAM_DIR}" ]] && \
   find "${SPACEJAM_DIR}" -mindepth 1 -print -quit | grep -q .
then
    BACKUP_DIR="${PROJECT_DIR}/results/archive/spacejam_${STAMP}"
    mkdir -p "${BACKUP_DIR}"
    cp -a "${SPACEJAM_DIR}/." "${BACKUP_DIR}/"
    echo "Existing SpaceJam outputs copied to: ${BACKUP_DIR}"
fi

mkdir -p "${SPACEJAM_DIR}"

echo "Starting SpaceJam training..."
python notebooks/03_spacejam_cell2location.py

echo "Checking expected SpaceJam outputs..."

expected_files=(
    "${SPACEJAM_DIR}/ADCAA_spot_factors_abs.pt"
    "${SPACEJAM_DIR}/ADCAA_spot_factors_rel.pt"
    "${SPACEJAM_DIR}/CTRL_spot_factors_abs.pt"
    "${SPACEJAM_DIR}/CTRL_spot_factors_rel.pt"
    "${SPACEJAM_DIR}/ADCAA_param_store.pt"
    "${SPACEJAM_DIR}/CTRL_param_store.pt"
    "${SPACEJAM_DIR}/ADCAA_training_loss.csv"
    "${SPACEJAM_DIR}/CTRL_training_loss.csv"
    "${SPACEJAM_DIR}/ADCAA_manifest_rois.csv"
    "${SPACEJAM_DIR}/ADCAA_manifest_factors.csv"
    "${SPACEJAM_DIR}/ADCAA_manifest_experiments.csv"
    "${SPACEJAM_DIR}/CTRL_manifest_rois.csv"
    "${SPACEJAM_DIR}/CTRL_manifest_factors.csv"
    "${SPACEJAM_DIR}/CTRL_manifest_experiments.csv"
    "${SPACEJAM_DIR}/ADCAA_run_metadata.json"
    "${SPACEJAM_DIR}/CTRL_run_metadata.json"
)

for f in "${expected_files[@]}"; do
    if [[ ! -f "${f}" ]]; then
        echo "Missing expected output: ${f}" >&2
        exit 1
    fi
done

echo "All expected SpaceJam outputs were found."

echo "====================================="
echo "Job finished successfully"
date
echo "====================================="
