#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
library(tidyverse)
library(Rsubread)

eclip_bams <- list.files("data/eclip_bams", pattern = ".bam", full.names = TRUE)
eclip_counts <- featureCounts(eclip_bams, 
                                    annot.ext = "/scratch/Shares/rinn/genomes/Homo_sapiens/Gencode/v38/gencode.v38.annotation.gtf",
                                    isGTFAnnotationFile = TRUE, isPairedEnd = TRUE)
write_rds(eclip_counts, "results/eclip_counts.rds")
