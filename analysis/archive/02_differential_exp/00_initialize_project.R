options(stringsAsFactors = FALSE)
library(tidyverse)
setwd("/scratch/Shares/rinn/Michael/firre_timecourse1/analysis/")
data_dir <- file.path("data", Sys.Date())
results_dir <- file.path("results", Sys.Date())
figures_dir <- file.path("figures", Sys.Date())
invisible(mapply(
  FUN = dir.create,
  path = c(data_dir, results_dir, figures_dir),
  MoreArgs = list(showWarnings = FALSE, recursive = TRUE)
))

### Import gene annotations

gtf <- rtracklayer::import("../../genomes/references/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf")

# Harmonize the naming scheme to match what bcbBioRNAseq package expects.
names(gtf@elementMetadata@listData)[
  which(names(gtf@elementMetadata@listData) == "gene_id")] <- "geneID"
names(gtf@elementMetadata@listData)[
  which(names(gtf@elementMetadata@listData) == "gene_name")] <- "geneName"
names(gtf@elementMetadata@listData)[
  which(names(gtf@elementMetadata@listData) == "gene_type")] <-"geneBiotype"
gtf@elementMetadata$id <- gtf@elementMetadata$geneID

g2s <- as.data.frame(gtf@elementMetadata@listData) %>% 
  select(geneID, geneName) %>% 
  distinct()

gtf <- gtf[which(gtf$type == "gene")]
