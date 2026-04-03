#!/bin/bash
#SBATCH -p test
#SBATCH --time=00:10:00
#SBATCH --job-name="Local"

module load HDF5/1.14.5-gompi-2024a
module load Python/3.12.3-GCCcore-13.3.0

python ../scripts/process_mag.py

cd ../scripts/Code/build
make -j4
cd ../../../slurm_tools
./../scripts/Code/build/main
