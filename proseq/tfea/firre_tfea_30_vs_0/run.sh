#!/bin/bash

###Name the job
#SBATCH --job-name=TFEA

###Specify the queue
#SBATCH --partition=short

###Specify WallTime
#SBATCH --time=24:00:00

### Specify the number of nodes/cores
#SBATCH --ntasks=10

### Allocate the amount of memory needed
#SBATCH --mem=50gb

### Setting to mail when the job is complete
#SBATCH --error /scratch/Shares/rinn/Michael/firre_timecourse/proseq/tfea/firre_tfea_30_vs_0/e_and_o/%x.err
#SBATCH --output /scratch/Shares/rinn/Michael/firre_timecourse/proseq/tfea/firre_tfea_30_vs_0/e_and_o/%x.out

### Set your email address
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mism6893@colorado.edu

### Load required modules
module purge

module load python/3.6.3
# module load R/3.6.1
module load python/3.6.3/matplotlib/1.5.1
module load python/3.6.3/scipy/0.17.1
module load python/3.6.3/numpy/1.14.1
module load python/3.6.3/htseq/0.9.1
module load python/3.6.3/pybedtools/0.7.10

module load bedtools/2.25.0
module load meme/5.0.3
module load samtools/1.3.1
module load gcc/7.1.0

### now call your program

source ~/anaconda3/etc/profile.d/conda.sh
conda activate tfea

### python ${src} --config ${config} --sbatch SUBMITTED
TFEA --output firre_tfea_30_vs_0 --combined_file /scratch/Shares/rinn/Michael/firre_timecourse/proseq/tfea/tfea_mumerge_step/combined_file.mumerge_MUMERGE.bed \
--bam1 /scratch/Shares/rinn/Michael/firre_timecourse/proseq/results/mapped/bams/JR3113.sorted.bam /scratch/Shares/rinn/Michael/firre_timecourse/proseq/results/mapped/bams/JR3116.sorted.bam /scratch/Shares/rinn/Michael/firre_timecourse/proseq/results/mapped/bams/JR3119.sorted.bam --bam2 /scratch/Shares/rinn/Michael/firre_timecourse/proseq/results/mapped/bams/JR3115.sorted.bam /scratch/Shares/rinn/Michael/firre_timecourse/proseq/results/mapped/bams/JR3118.sorted.bam /scratch/Shares/rinn/Michael/firre_timecourse/proseq/results/mapped/bams/JR3121.sorted.bam --label1 0min --label2 30min --genomefasta /scratch/Shares/rinn/genomes/Mus_musculus/hisat2/GRCm38.fa --fimo_motifs /scratch/Shares/rinn/Michael/firre_timecourse/proseq/bin/motif_databases/MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme --mdd
