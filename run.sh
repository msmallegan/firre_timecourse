#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=ftc-rna-seq-nextflow
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb                 # Memory limit
#SBATCH --time=100:00:00            # Time limit hrs:min:sec
#SBATCH --output=./log/nextflow.out
#SBATCH --error=./log/nextflow.err



pwd; hostname; date

echo "You've requested $SLURM_CPUS_ON_NODE core."

export PATH=/opt/singularity/3.1.1/bin:$PATH

nextflow run nf-core/rnaseq -r 1.3 \
-resume \
--reads 'fastq/*{_read1,_read2}.fastq.gz' \
--skip_genebody_coverage \
--saveAlignedIntermediates \
--genome GRCm38 \
--email michael.smallegan@colorado.edu \
--igenomes_base '/scratch/Shares/rinn/Michael/genomes/references/' \
-c nextflow.config

date
