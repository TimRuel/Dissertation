#!/bin/bash
#SBATCH -A p32397
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mail-user=timothyruel2024@u.northwestern.edu
#SBATCH --output=out-slurm.txt
#SBATCH --job-name testjob
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --array=0-4
#SBATCH --mem=1G  # Total memory in GB needed for a job. Also see --mem-per-cpu

module purge all

module load R/4.4.0

export MC_CORES=$(($SLURM_NTASKS-1))
export TASK_ID=$(($SLURM_ARRAY_TASK_ID))

R --max-connections=256 --vanilla < test.R

# submit this script with sbatch from the command line:
# sbatch job_submit.sh