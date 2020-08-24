options(stringsAsFactors = FALSE)
library(tidyverse)

# let's read in the results tables from the differential expresion analysis
fl <- list.files("results/differential_expression/", full.names = TRUE)

dat <- read.csv(fl[1])
dat$gene_set <- fl[1]
for(i in 2:length(fl)) {
  tdat <- read.csv(fl[i])
  tdat$gene_set <- fl[i]
  dat <- bind_rows(dat, tdat)
}

dat$gene_set <- gsub("results/differential_expression//", "", dat$gene_set)
dat$gene_set <- gsub(".csv", "", dat$gene_set)

dat <- dat %>% select(-X)
write.csv(dat, "results/combined_significance.csv", row.names = FALSE)
dat <- dat %>% select(gene_id, gene_name, padj, gene_set)
# dat <- dat %>% spread(key = gene_set, value = padj)

# let's use a cutoff of 0.05
sigdat <- dat %>% filter(padj < 0.05)
sigdat$sig <- 1
sigdat <- sigdat %>% select(-padj)
sigdat <- sigdat %>% spread(key = gene_set, value = sig, fill = 0)

# let's look at the stemness genes
stemness_genes <- c("Nanog", "Pou5f1", "Sox2",
                    "Klf4", "Zfp42", "Myc")
stem <- sigdat[which(sigdat$gene_name %in% stemness_genes), ]

stem <- stem %>% select(-gene_id)
stem <- stem %>% gather(key = "gene_set", value = "sig", 2:8)

stem$sig <- factor(stem$sig)
stem <- stem %>% filter(gene_set %in% c("firre_responders","firre_responders_same_in_both",
                                        "all_genes_changing_in_time", 
                                        "genes_changing_in_ko"))
g <- ggplot(stem, aes(x = gene_name, y = gene_set, color = sig))
g + geom_point() + scale_color_manual(values = c("white", "#A81E2B")) + 
  theme(axis.text.x = element_text(angle = -35, hjust =-.1, vjust = 1))
ggsave("stem_markers_setplot.pdf")
