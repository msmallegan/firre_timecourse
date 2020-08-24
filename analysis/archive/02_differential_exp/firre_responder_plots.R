options(stringsAsFactors = FALSE)
library(tidyverse)
source("_setup.R")
source("helper_functions.R")

samples <- read.csv("../../samplesheet.csv") %>%
  filter(cell_type == "mESCs",
         !grepl("A", sample_name)) %>%
  mutate(timepoint = factor(timepoint, 
                            levels = as.character(seq(from = 0, to = 330, by = 30))))
# Change from a boolean to a string
samples[which(samples$firre_induced), "firre_induced"] <- "firre_induced"
samples[which(samples$firre_induced == FALSE), "firre_induced"] <- "firre_uninduced"
samples$firre_induced <- factor(samples$firre_induced, levels = c("firre_uninduced", "firre_induced"))

# Factorfy the other design variable
samples$firre_ko <- factor(samples$firre_ko, levels = c("WT", "KO"))



rownames(samples) <- samples$id
counts_combined <- read.table("../../results/featureCounts/merged_gene_counts.txt", header = T)
g2s <- counts_combined[,c(1,2)]
names(g2s) <- c("gene_id", "gene_name")

# Read in the featureCounts
names(counts_combined) <- sapply(names(counts_combined), function(x) {unlist(strsplit(x, "_"))[[1]]})
rownames(counts_combined) <- counts_combined$Geneid
counts <- as.matrix(counts_combined[,3:ncol(counts_combined)])

# Ensure that the ordering of the columns in the counts matrix 
# is the same as in the sample sheet.
# let's subset to just the samples which we would like to use
counts <- counts[,samples$id]
stopifnot(all(rownames(samples) == colnames(counts)))



fcounts <- fCountReader("../../results/featureCounts/gene_counts/",
                       samples$id[1], 
                       "_read1Aligned.sortedByCoord.out_gene.featureCounts.txt")

countConverter<-function(fCount,return="TPM") {
  if (return=="FPKM") {
    out=fCount$count/(fCount$annot$length/1000)
    out=t(t(out)/(colSums(fCount$count)/10^6))
  } else if (return=="TPM") {
    out=fCount$count/(fCount$annot$length/1000)
    out=t(t(out)/(colSums(out)/10^6))
  }
  return(out)
}

tpm <- counts/(fcounts$annot$length/1000)
tpm <- t(t(tpm)/(colSums(tpm)/10^6))

row.names(tpm) <- row.names(counts)
colnames(tpm)

# write.csv(tpm, "mesc_firre_tc_tpm.csv")
fr_sb_overlap_96 <- read.csv("../../../ftc_long/analysis/fr_sb_overlap_96.csv")
# fr <- read.csv("results/differential_expression/firre_responders_same_in_both.csv")

# fr <- fr %>%
#   filter(padj < 0.05)



ex_counts <- tpm[fr_sb_overlap_96$gene_id,] %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key = "id", value = "tpm", 2:ncol(.)) %>%
  merge(g2s) %>%
  merge(samples)

ex_counts <- merge(ex_counts, fr_sb_overlap_96 %>% select(-gene_name))
ex_counts <- ex_counts %>% arrange(padj)
ex_counts$gene_name <- factor(ex_counts$gene_name, levels = unique(ex_counts$gene_name))

g <- ggplot(ex_counts, aes(x = timepoint, y = tpm, color = firre_induced, group = firre_induced))
g + geom_point(alpha = 0.6) + stat_summary(fun.y=mean, geom="line") +
  # scale_y_log10() + 
  facet_wrap(~ gene_name, scales = "free_y") + 
  # scale_color_manual(values = c("#414142", "#A83F4B"))
scale_color_manual(values = c("#414142", "#A81E2B"))


# paste(fr$gene_name, collapse = ", ")

# Okay, now let's plot the ES marker genes
# Nanog, Pou5f1 (Oct4), Sox2, Klf4, Zfp42 (Rex1), Myc
stemness_genes <- c("Nanog", "Pou5f1", "Sox2",
                    "Klf4", "Zfp42", "Myc")

stemness_gene_ids <- g2s[which(g2s$gene_name %in% stemness_genes),"gene_id"]
ex_counts <- tpm[stemness_gene_ids,] %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key = "id", value = "tpm", 2:ncol(.)) %>%
  merge(g2s) %>%
  merge(samples)

# ex_counts <- merge(ex_counts, fr %>% select(-gene_name))
# ex_counts <- ex_counts %>% arrange(padj)
ex_counts$gene_name <- factor(ex_counts$gene_name, levels = unique(ex_counts$gene_name))

g <- ggplot(ex_counts, aes(x = timepoint, y = tpm, color = firre_induced, group = firre_induced))
g + geom_point(alpha = 0.6) + stat_summary(fun.y=mean, geom="line") + facet_wrap(~ gene_name, scales = "free_y") + 
  # scale_color_manual(values = c("#414142", "#A83F4B"))
  scale_color_manual(values = c("#414142", "#A81E2B")) + 
  theme(axis.text.x = element_text(angle = -35, hjust =-.1, vjust = 1))
# ggsave("stem_markers_profile_plots.pdf", height = 7, width = 12)


# Now let's look at how Firre compares to other highly expressed genes
# to consider the RNA effect
which.max(tpm)
totals <- rowSums(tpm) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  merge(g2s)



# let's choose a few genes to plot
gtp <- c("Firre", "Nanog", "Dppa5a", "Npm1", "Gapdh")
gtpids <- g2s[which(g2s$gene_name %in% gtp),"gene_id"]
ex_counts <- tpm[gtpids,] %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key = "id", value = "tpm", 2:ncol(.)) %>%
  merge(g2s) %>%
  merge(samples)

g <- ggplot(ex_counts, aes(x = gene_name, y = tpm, group = gene_name)) 
g + 
  geom_jitter(width = 0.1, color = "#A81E2B", alpha = 0.3) +
  theme(axis.text.x = element_text(angle = -35, hjust =-.1, vjust = 1))
ggsave("highly_expressed_genes.pdf")
# g <- ggplot(ex_counts, aes(x = timepoint, y = tpm, color = firre_induced, group = firre_induced))
# g + geom_point(alpha = 0.6) + stat_summary(fun.y=mean, geom="line") + facet_wrap(~ gene_name, scales = "free_y") + 
#   # scale_color_manual(values = c("#414142", "#A83F4B"))
#   scale_color_manual(values = c("#414142", "#A81E2B")) + 
#   theme(axis.text.x = element_text(angle = -35, hjust =-.1, vjust = 1))
