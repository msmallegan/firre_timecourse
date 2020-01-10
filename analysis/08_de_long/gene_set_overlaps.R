options(stringsAsFactors = FALSE)


dat <- read.csv("../../firre_timecourse/analysis/02_differential_exp/results/combined_significance.csv")
ftc96 <- read.csv("results/firre_responders_96.csv")

unique(dat$gene_set)
fr_tc <- dat %>% filter(gene_set == "firre_responders") %>%
  filter(padj <= 0.05)
fr_96 <- ftc96 %>% filter(padj <= 0.05)


length(which(fr_96$gene_id %in% fr_tc$gene_id))
length(which(fr_tc$gene_id %in% fr_96$gene_id))

library(VennDiagram)

v <- venn.diagram(list("fr_tc" = fr_tc$gene_id, "fr_96" = fr_96$gene_id), 
                  filename = NULL)
grid.newpage()
grid.draw(v)
# intersect(fr_96$gene_id,fr_tc$gene_id)


fr_sb <- dat %>% filter(gene_set == "firre_responders_same_in_both") %>%
  filter(padj <= 0.05)
v <- venn.diagram(list("fr_sb" = fr_sb$gene_id, "fr_96" = fr_96$gene_id), 
                  filename = NULL)
grid.newpage()
grid.draw(v)

fr_sb_overlap_96 <- fr_sb[which(fr_sb$gene_id %in% intersect(fr_96$gene_id,fr_sb$gene_id)),]
fr_sb_overlap_96$gene_set <- "fr_sb_overlap_96"
write.csv(fr_sb_overlap_96, "fr_sb_overlap_96.csv")

# q = the number of white balls drawn from the urn (without replacement)
# 
# m = the number of white balls in the urn
# 
# n = the number of black balls in the urn
# 
# k = the number of balls drawn from the urn (sample size)
# total number of genes
tg <- dat %>% filter(baseMean > 10)
length(unique(tg$gene_id))
15948
229

1 - phyper(29,82,15948-229,229)
