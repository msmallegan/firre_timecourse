#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
library(tidyverse)
library(Rsubread)

eclip_bams <- list.files("data/eclip_bams", pattern = ".bam", full.names = TRUE)
eclip_firre_counts <- featureCounts(eclip_bams, 
                                    annot.ext = "/scratch/Shares/rinn/genomes/Homo_sapiens/Gencode/v38/firre.gtf",
                                    isGTFAnnotationFile = TRUE, isPairedEnd = TRUE)
write_rds(eclip_firre_counts, "results/eclip_firre_counts.rds")