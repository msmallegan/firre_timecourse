options(stringsAsFactors = FALSE)

of <- read.table("old_files.txt")

nfss <- read.csv("JR2156_firre_overexpression.csv")

nf <- read.table("new_files.txt") %>% 
  separate(V1, into = c("sample_name", "read", "suffix"), sep = "_", remove = FALSE) %>%
  filter(!is.na(suffix))
nf <- merge(nf, nfss)
nf$suffix <- ".fastq.gz"
nf$read <- gsub("R", "_read", nf$read)
nf$lane <- "_lane1"
nf <- nf %>% unite(new_name, id, lane, read, suffix, sep = "")
nf <- nf %>% select(V1, new_name)
names(nf)[1] <- "old_name"
nf$old_path <- "/rinnlab/rinngrp/mism6893/christian_timecourse/30-239079954/"
nf$new_path <- "/scratch/Shares/rinn/Michael/firre_timecourse1/fastq/"
nf <- nf %>% unite(op, old_path, old_name, sep = "")
nf <- nf %>% unite(np, new_path, new_name, sep = "")



## old files
of <- of %>% separate(V1, into = c("id", "lane", "read", "suffix"), remove = F)
of$lane <- gsub("L00", "_lane", of$lane)
of$read <- gsub("R", "_read", of$read)
of$suffix <- ".fastq.gz"
of <- of %>% unite(new_name, id, lane, read, suffix, sep = "")
of$new_path <- "/scratch/Shares/rinn/Michael/firre_timecourse1/fastq/"
of$old_path <- "/scratch/Shares/rinn/Michael/firre_timecourse/input/fastq/"
of <- of %>% unite(op, old_path, V1, sep = "")
of <- of %>% unite(np, new_path, new_name, sep = "")


copy_map <- bind_rows(of, nf)

# hmm, let's remove the lane1 for simplicity since they're all lane1
copy_map$np <- gsub("_lane1", "", copy_map$np)
write.table(copy_map, "copy_map.csv", row.names = F, col.names = F, quote = F, sep = ",")
