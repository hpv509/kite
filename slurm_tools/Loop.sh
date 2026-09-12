#!/bin/bash
THS=4

for ID in {1..400}; do
	sbatch --cpus-per-task="$THS" --export=ALL,ths="$THS",id="$ID" Launch.sh
done
