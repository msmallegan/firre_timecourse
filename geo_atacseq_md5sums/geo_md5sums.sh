#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=md5sums_geo
#SBATCH --mail-type=NONE
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-343
#SBATCH --mem=6gb
#SBATCH --time=0:20:00
#SBATCH --output=md5sums_geo_%A_%a.out
#SBATCH --error=md5sums_geo_%A_%a.err

hostname
date

FILES=(../geo_upload/geo_submission_atacseq_may4/*.fastq.gz)    
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo ${FILE}
OUTFILE=$(basename ${FILE} .fastq.gz).md5.txt
md5sum $FILE > $OUTFILE

date
