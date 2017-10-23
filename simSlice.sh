#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=sliceMaking-%j.out
#SBATCH --time=672:00:00
#SBATCH --ntasks=100
#SBATCH --mem-per-cpu=5000
#SBATCH -p scavenge

module load R

Rscript ./sliceMakingOnSIMS.R

