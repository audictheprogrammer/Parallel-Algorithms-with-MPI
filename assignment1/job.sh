#!/bin/bash
#SBATCH --nodes=1 --time=06:00:00
#SBATCH --partition=singlenode
#SBATCH --job-name=TEST01
#SBATCH --export=NONE

# run job using sbatch job.sh
# check with squeue
module load intel
./exe-ICC