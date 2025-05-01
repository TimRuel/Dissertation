#!/bin/bash
#SBATCH -A p32397
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothyruel2024@u.northwestern.edu
#SBATCH --output=sample_job.%A_%a.out
#SBATCH ---job-name="sample_job_\${SLURM_ARRAY_TASK_ID}"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --array=0-999
#SBATCH --mem=1G  # Total memory in GB needed for a job. Also see --mem-per-cpu

module purge all

module load R/4.4.0

module load hdf5/1.14.1-2-gcc-12.3.0 
module load gsl/2.7.1-gcc-12.3.0 
module load fftw/3.3.10-gcc-12.3.0 
module load gdal/3.7.0-gcc-12.3.0

export MC_CORES=$(($SLURM_NTASKS-1))
export TASK_ID=$(($SLURM_ARRAY_TASK_ID))

R --max-connections=256 --vanilla < test.R

# submit this script with sbatch from the command line:
# sbatch job_submit.sh