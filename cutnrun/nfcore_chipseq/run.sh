#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=cutnrun_chip_nfcore
#SBATCH --mail-type=NONE
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=40:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date

module load singularity/3.1.1

nextflow run nf-core/chipseq -r 1.2.1 \
-profile singularity \
--input design.csv \
--fasta /scratch/Shares/rinn/genomes/Mus_musculus/Gencode/M25/GRCm38.primary_assembly.genome.fa \
--gtf /scratch/Shares/rinn/genomes/Mus_musculus/Gencode/M25/gencode.vM25.primary_assembly.annotation.gtf \
--macs_gsize 1.87e9 \
--blacklist mm10-blacklist.v2.bed \
--email michael.smallegan@colorado.edu \
-resume \
-c nextflow.config

date
