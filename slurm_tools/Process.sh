#!/bin/bash
#SBATCH -p test
#SBATCH --time=00:30:00
#SBATCH --job-name="Cond"
#SBATCH --mem-per-cpu=4GB

module load HDF5/1.14.5-gompi-2024a
module load Python/3.12.3-GCCcore-13.3.0

python ../scripts/process_swave.py
