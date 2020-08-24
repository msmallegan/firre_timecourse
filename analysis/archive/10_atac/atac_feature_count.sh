#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=fc_atac
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=50gb
#SBATCH --time=2:00:00
#SBATCH --output=fc.out
#SBATCH --error=fc.err

date

# Let's run this r script
./atac_feature_count.R

date