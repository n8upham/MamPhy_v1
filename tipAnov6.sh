#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=speciesAnova-%j.out
#SBATCH --time=168:00:00
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=46000
#SBATCH -p bigmem

module load R

Rscript ./Fig4_tipLevelDR_1resp_multiANOVA_forHPC_6var.R

