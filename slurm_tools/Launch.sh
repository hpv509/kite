#!/bin/bash
#SBATCH -p nodes,preempt
#SBATCH --time=00:30:00
#SBATCH --account=pet-kite-2023
#SBATCH --job-name="Local"
#SBATCH --mem-per-cpu=4G

module load HDF5/1.14.5-gompi-2024a
module load Eigen/3.4.0-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0
module load GCC/13.3.0

output=$(python ../examples/ldos.py "$seed" 2>&1 >/dev/null)
./../build/KITEx "$output"
