#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=ftc-rnaseq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=100:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core."

module load singularity/3.1.1

nextflow run nf-core/rnaseq -r 1.4.2 \
-resume \
-profile singularity \
--reads 'fastq/*{_read1,_read2}.fastq.gz' \
--fasta ../../../genomes/Mus_musculus/Gencode/M25/GRCm38.p6.genome.fa \
--gtf ../../../genomes/Mus_musculus/Gencode/M25/gencode.vM25.annotation.gtf \
--pseudo_aligner salmon \
--gencode \
--skipQC \
--email michael.smallegan@colorado.edu \
-c nextflow.config

date
