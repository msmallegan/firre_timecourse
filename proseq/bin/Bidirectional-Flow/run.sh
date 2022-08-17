#!/bin/bash
#SBATCH --job-name=firre_bidir # Job name
#SBATCH -p long
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mism6893@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=8gb # Memory limit
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=nextflow.%j.out # Standard output
#SBATCH --error=nextflow.%j.err # Standard error log

module load samtools/1.8
module load bedtools/2.28.0
module load openmpi/1.6.4
module load gcc/7.1.0
module load python/3.6.3


nextflow run main.nf -profile mm10 \
--bams "/scratch/Shares/rinn/Michael/firre_timecourse/proseq/results/mapped/bams/*.bam" \
--workdir /scratch/Shares/rinn/Michael/firre_timecourse/proseq/bidir_work/ \
--outdir /scratch/Shares/rinn/Michael/firre_timecourse/proseq/bidir_results \
--tfit \
--gene_count \
--savebidirs \
-resume
