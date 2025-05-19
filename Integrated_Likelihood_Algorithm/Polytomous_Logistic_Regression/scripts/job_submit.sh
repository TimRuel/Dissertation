#!/bin/bash
#SBATCH -A p32397
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothyruel2024@u.northwestern.edu
#SBATCH --output=logs/sample_job.%A_%a.out
#SBATCH --job-name="sample_job_%a"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --array=0-999
#SBATCH --mem=1G

# Load modules
module purge all
module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0 
module load gsl/2.7.1-gcc-12.3.0 
module load fftw/3.3.10-gcc-12.3.0 
module load gdal/3.7.0-gcc-12.3.0

# Set parallelism environment variable if needed
export MC_CORES=$((SLURM_NTASKS - 1))

# Set experiment ID and iteration
EXPERIMENT_ID="experiment_A"
ITERATION_ID=$SLURM_ARRAY_TASK_ID
REQUESTED_CORES=$MC_CORES

# Run the R script using Rscript (better for command-line args)
Rscript --max-connections=256 scripts/main.R "$EXPERIMENT_ID" "$ITERATION_ID" "$REQUESTED_CORES"
