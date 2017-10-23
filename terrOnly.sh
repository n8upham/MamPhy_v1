#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=speciesLevelDR-%j.out
#SBATCH --time=168:00:00
#SBATCH --ntasks=25
#SBATCH --mem-per-cpu=5000

module load R

Rscript ./Fig4_spLevelDR_terrOnly_8vars_forHPC.R

