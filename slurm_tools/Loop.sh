#!/bin/bash
THS=4
MU=0.45
UD=0.65

for ID in {1..100}; do
	sbatch --cpus-per-task="$THS" --export=ALL,ths="$THS",id="$ID",ud="$UD",mu="$MU" Launch.sh
done
