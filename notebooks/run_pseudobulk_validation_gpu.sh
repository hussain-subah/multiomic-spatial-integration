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

echo "====================================="
echo "Job started on $(hostname)"
date
echo "====================================="

cd /N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration

module purge
module load python/gpu/3.10.10
source /N/soft/rhel8/conda/26.3.2/etc/profile.d/conda.sh
conda activate /N/u/echimal/Quartz/.conda/envs/integration_env

export PYTHONPATH=$(pwd)

export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16

echo "Python:"
which python
python --version

echo "GPU check:"
python -c "import torch; print('CUDA available:', torch.cuda.is_available()); print('GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')"

echo "Starting pseudobulk validation..."
python notebooks/06_pseudobulk_validation.py

echo "====================================="
echo "Job finished"
date
echo "====================================="
