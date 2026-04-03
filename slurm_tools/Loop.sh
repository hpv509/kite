#!/bin/bash
THS=4
Seed=1

sbatch --cpus-per-task="$THS" --export=ALL,seed="$Seed" Launch.sh
