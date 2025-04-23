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

# unload modules that may have been loaded when job was submitted
module purge all

module load R/4.4.0

module load gdal/3.1.3-R-4.1.1
module load proj/7.1.1
module load geos/3.8.1
module load hdf5/1.10.6-openmpi-3.1.3-gcc-8.4.0

# set environment variable for R for parallelization;
# MC_CORES should be 1 less that cores requested above (which is 4),
# the value of $SLURM_NPROCS is set to the same value as --ntasks-per-node above
# since it is additional processes on top of the main one
export MC_CORES=$(($SLURM_NPROCS-1))
export TASK_ID=$(($SLURM_ARRAY_TASK_ID))

Rscript test.R 

# submit this script with sbatch from the command line:
# sbatch job_submit.sh