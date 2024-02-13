#!/bin/bash
#SBATCH --job-name=seg_sites
#SBATCH -A standby
#SBATCH --mem=20G
#SBATCH --array=1-100
#SBATCH -t 4:00:00

module purge
module load r/4.2

outfile=mcmc_res

Rscript seg_sites_rarefaction.R $SLURM_ARRAY_TASK_ID $outfile
