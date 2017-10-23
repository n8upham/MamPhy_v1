#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=cladeLevel-%j.out
#SBATCH --time=672:00:00
#SBATCH --ntasks=100
#SBATCH --mem-per-cpu=6400
#SBATCH -p general

module load R

Rscript ./Fig4_cladeLevelDR_4resp_uni10pred_bySlice_forHPC.R

