#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=tfea
#SBATCH --mail-type=NONE
#SBATCH --mail-user=michael.smallegan@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=3:00:00
#SBATCH --output=tfea.out
#SBATCH --error=tfea.err

module load singularity/3.1.1



singularity exec \
--bind /scratch/Shares/rinn/Michael/firre_timecourse/analysis/04_tfea/results/firre_induction_vs_control_rankings:/mnt/rankings \
--bind /scratch/Shares/rinn/Michael/firre_timecourse/analysis/04_tfea/results/tfea_results/ESC_KO_tp30:/mnt/output \
--bind /scratch/Shares/rinn/Michael/firre_timecourse/analysis/04_tfea/data/motif_databases/MOUSE/:/mnt/mouse_meme \
--bind /scratch/Shares/rinn/genomes/Mus_musculus/Gencode/M25:/mnt/genome \
	bin/tfea.sif TFEA \
--output /mnt/output \
--ranked_file /mnt/rankings/ESC_KO_tp30.bed \
--label1 firre_induced --label2 dox_control \
--genomefasta /mnt/genome/GRCm38.p6.genome.fa \
--fimo_motifs /mnt/mouse_meme/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme
