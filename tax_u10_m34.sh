#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=taxonLevel-%j.out
#SBATCH --time=672:00:00
#SBATCH --ntasks=100
#SBATCH --mem-per-cpu=5000
#SBATCH -p general

module load R

Rscript ./Fig3_TaxonLevelRich_uni10var_forHPC.R

