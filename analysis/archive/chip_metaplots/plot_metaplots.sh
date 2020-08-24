# Plot metaplots
module load samtools
samtools merge H3K27me3_ES.bam H3K27me3_ES_ENCFF006QMD.bam H3K27me3_ES_ENCFF041GIV.bam
samtools merge H3K4me3_ES.bam H3K4me3_ES_ENCFF017QQB.bam H3K4me3_ES_ENCFF246IWX.bam


ngs.plot.r -G mm10 -R tss -C frik4.config.txt -O H3K4me3_fr -D ensembl
ngs.plot.r -G mm10 -R tss -C fri.config.txt -O H3K27me3_fr -D ensembl