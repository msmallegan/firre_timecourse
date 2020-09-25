#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=ftc-atacseq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=50:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core."

module load singularity/3.1.1

nextflow run nf-core/atacseq -r 1.2.0 \
-resume \
-profile singularity \
--input atacseq_design.csv \
--fasta ../../../genomes/Mus_musculus/Gencode/M25/GRCm38.p6.genome.fa \
--gtf ../../../genomes/Mus_musculus/Gencode/M25/gencode.vM25.annotation.gtf \
--blacklist mm10-blacklist.v2.bed \
--macs_gsize 1.87e9 \
--email michael.smallegan@colorado.edu \
-c nextflow.config

date
