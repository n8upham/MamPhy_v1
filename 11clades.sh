#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=speciesLevelDR-%j.out
#SBATCH --time=168:00:00
#SBATCH --ntasks=25
#SBATCH --mem-per-cpu=5G

module load R

Rscript ./Fig4_spLevelDR_11clades_8vars_forHPC.R

