options(stringsAsFactors = FALSE)
fl <- list.files("fastq/", full.names = T)
library(tidyverse)
samples <- read.csv("rnaseq/samplesheet.csv")
samples <- samples %>% filter(timepoint %in% c("0m", "0h"))
names(fl) <- sapply(fl, function(x) {
  unlist(strsplit(unlist(strsplit(x, "//"))[[2]], "_"))[[1]]
})
names(fl)
files <- data.frame("sample_id" = names(fl),
                    "file_path" = paste0("/scratch/Shares/rinn/Michael/firre_timecourse/", fl))
samples <- samples %>% merge(files)
write_csv(samples, "firre_timecouse_zero_timepoint_samples.csv")
