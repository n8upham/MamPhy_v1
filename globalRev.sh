#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=speciesLevelDR-%j.out
#SBATCH --time=672:00:00
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=5000
#SBATCH -p pi_jetz

module load R

Rscript ./Fig4_spLevelDR_globalRev_8vars_forHPC.R

