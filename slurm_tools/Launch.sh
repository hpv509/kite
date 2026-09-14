#!/bin/bash
#SBATCH -p nodes
#SBATCH --time=01:00:00
#SBATCH --job-name="Test"
#SBATCH --mem-per-cpu=4GB

module load HDF5/1.14.5-gompi-2024a
module load Eigen/3.4.0-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0
module load FFTW/3.3.10-GCC-13.3.0
module load GCC/13.3.0


output=$(python ../examples/rhomb_graphene_swave.py "$id" "$ud" "$mu" 2>&1 >/dev/null)
OMP_NUM_THREADS="$ths" ./../build/KITEx "$output"
