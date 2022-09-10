#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=bw2bg
#SBATCH --mail-type=NONE
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-8
#SBATCH --mem=10gb
#SBATCH --time=2:00:00
#SBATCH --output=logs/bw2bg_%A_%a.out
#SBATCH --error=logs/bw2bg_%A_%a.err

hostname; date

module load samtools

declare -a DBPS PEAK_FILES BATCH SLURM_ARRAY
while IFS=$'\t' read -r first second;do
BW_FILE+=("$first")
TARGET_FILE+=("$second")
done < /scratch/Shares/rinn/Michael/firre_timecourse/analysis/15_coverage_tracks/results/proseq_bg_files.tsv

bgzip < ${BW_FILE[$SLURM_ARRAY_TASK_ID]} > ../results/bedgraphs/${TARGET_FILE[$SLURM_ARRAY_TASK_ID]}.bedGraph.bgz
tabix -p bed ../results/bedgraphs/${TARGET_FILE[$SLURM_ARRAY_TASK_ID]}.bedGraph.bgz

date