#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=cladeLevel-%j.out
#SBATCH --time=672:00:00
#SBATCH --ntasks=100
#SBATCH --mem-per-cpu=5000
#SBATCH -p scavenge

module load R

Rscript ./Fig4_cladeLevelDR_4resp_26pred_5to50Ma_forHPC.R

