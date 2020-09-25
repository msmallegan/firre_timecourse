library(tidyverse)
# Create a script to move files to scratch
if(!dir.exists("../../atacseq/fastq")) {
  cat("Creating fastq directory") 
  dir.create("../../atacseq/fastq")
}

samples <- read.csv("../../atacseq/atacseq_samplesheet.csv")
samples$dest_read1 <- paste0(file.path("/scratch/Shares/rinn/Michael",
                                       "firre_timecourse/atacseq/fastq/"),
                             samples$sample_id, "_read1.fastq.gz")
samples$dest_read2 <- paste0(file.path("/scratch/Shares/rinn/Michael",
                                       "firre_timecourse/atacseq/fastq/"),
                             samples$sample_id, "_read2.fastq.gz")

cp_read1 <- samples %>% dplyr::select(raw_data_read1,
                                      dest_read1)
names(cp_read1) <- c("src", "dest")
cp_read2 <- samples %>% dplyr::select(raw_data_read2,
                                      dest_read2)
names(cp_read2) <- c("src", "dest")
cp_df <- bind_rows(cp_read1, cp_read2)
cp_df$cmd <- "cp"
cp_df <- cp_df %>% dplyr::select(cmd, src, dest)

# Write script to move FASTQs to scratch
cat("#!/bin/bash \n", file = "copy_atacseq_fastqs_to_scratch.sh")
write.table(cp_df, "copy_atacseq_fastqs_to_scratch.sh", append = TRUE,
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Write a file to check MD5 sums
md5_read1 <- samples %>% 
  dplyr::select(md5_read1, dest_read1)
names(md5_read1) <- c("md5sum", "file")
md5_read2 <- samples %>% 
  dplyr::select(md5_read2, dest_read2)
names(md5_read2) <- c("md5sum", "file")
md5_df <- bind_rows(md5_read1, md5_read2)
write.table(md5_df, "atacseq_md5_list.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")


# Check that all files currently exist
d1 <- samples$dest_read1
d2 <- samples$dest_read2
dest <- c(d1, d2)

for(i in 1:length(dest)) {
  if(!file.exists(dest[i])) {
    cat("file doesn't exist: ", dest[i])
  }
}


# Retrieve blacklist
system(paste0("cd /scratch/Shares/rinn/Michael/firre_timecourse/atacseq/; ",
              "wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz; ",
              "gunzip mm10-blacklist.v2.bed.gz"))



