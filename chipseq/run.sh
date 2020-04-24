#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=firre_chip
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=12:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "Nice -- you've requested $SLURM_CPUS_ON_NODE core."

module load singularity/3.1.1

nextflow run nf-core/chipseq -r 1.1.0 \
-profile singularity \
--input design.csv \
--fasta ../util/GRCm38.p6.genome.fa \
--gtf ../../genomes/Mus_musculus/Gencode/M24/gencode.vM24.annotation.gtf \
--macs_gsize 2.6e9 \
--blacklist mm10-blacklist.v2.bed \
--save_reference \
--email michael.smallegan@colorado.edu \
-resume \
-c nextflow.config

date
