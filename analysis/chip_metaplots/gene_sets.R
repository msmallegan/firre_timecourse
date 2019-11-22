options(stringsAsFactors = FALSE)
library(tidyverse)

dat <- read.csv("../02_differential_exp/results/combined_significance.csv")
gl <- dat %>% filter(padj < 0.05) %>%
  select(gene_id, gene_set)
write.table(file = "fr_gene_names.txt",
  gl[which(gl$gene_set == "firre_responders"),"gene_name"],
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
  )

write.table(file = "frall_gene_names.txt",
            gl[which(gl$gene_set == "all_genes_changing_in_time"),"gene_name"],
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
)

# for ngsplot I'm going to need the decimal taken off the end 
# of the gene id
gl$gene_id <- sapply(gl$gene_id, function(x) {
  unlist(strsplit(x,".", fixed = T))[[1]]
})
gene_sets <- unique(gl$gene_set)

for(i in 1:length(gene_sets)) {
  write.table(file = paste0(gene_sets[[i]],".Ensembl.txt"),
              gl[which(gl$gene_set == gene_sets[[i]]),"gene_id"],
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}



# create a background gene set
hist(dat$baseMean, xlim = c(0,10))
bg <- dat %>% filter(baseMean > 10) %>%
  select(gene_id) %>%
  distinct()

bg$gene_id <- sapply(bg$gene_id, function(x) {
  unlist(strsplit(x,".", fixed = T))[[1]]
})
write.table(file = paste0("background_baseMean_gt10",".Ensembl.txt"),
            bg[,"gene_id"],
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
