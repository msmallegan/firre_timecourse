#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=rnaseq2bw
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-227
#SBATCH --mem=8gb
#SBATCH --time=2:00:00
#SBATCH --output=rnaseq2bw_%a.out
#SBATCH --error=rnaseq2bw_%a.err

hostname; date

declare -a BAM BW
while IFS=$'\t' read -r first second;do
BAM+=("$first")
BW+=("$second")
done < /scratch/Shares/rinn/Michael/firre_timecourse/analysis/15_coverage_tracks/data/rnaseq2bw.txt

source ~/anaconda3/etc/profile.d/conda.sh
conda activate deeptools

module load samtools
samtools index ${BAM[$SLURM_ARRAY_TASK_ID]}

bamCoverage -b ${BAM[$SLURM_ARRAY_TASK_ID]} -o /scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/${BW[$SLURM_ARRAY_TASK_ID]}

date

