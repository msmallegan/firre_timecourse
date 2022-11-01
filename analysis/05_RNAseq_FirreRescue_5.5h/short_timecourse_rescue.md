Short timecourse rescue
================

For this experiment, we induced the expression of Firre with an rTTA
element by adding doxycycline to the cells. We see that Firre is indeed
expressed in the KO background after the addition of doxycycline. The
drug does instigate some gene expression changes on its own, so we will
control for the effects by using a linear model which accounts for the
effect of dox.

``` r
if(!file.exists("results/ko_rescue_short.RData")) {
  
  # Filter to ESC KO long timecourse
  ko_rescue_short_samples <- samples %>%
    filter(cell_type == "ESC",
           timecourse_length == "short",
           firre_ko == "KO")
  ko_rescue_short_counts <- salmon_gene_counts[,ko_rescue_short_samples$sample_id]
  
  # Check ordering
  stopifnot(all(rownames(ko_rescue_short_samples) == colnames(ko_rescue_short_counts)))
  stopifnot(all(rownames(ko_rescue_short_counts) == genes$gene_id))
  
  # DESeq2 -- controlling for doxycycline; likelihood ratio test
  ko_rescue_short_dds <- DESeqDataSetFromMatrix(countData = ko_rescue_short_counts, 
                                                colData = ko_rescue_short_samples, 
                                                design = ~ firre_induced + timepoint + timepoint*firre_induced)
  ko_rescue_short_dds <- DESeq(ko_rescue_short_dds, test = "LRT", reduced = ~ firre_induced + timepoint)
  
  # Compile results
  res_names <- resultsNames(ko_rescue_short_dds)
  dynamic_res <- res_names[grepl("firre_inducedfirre_induced.timepoint", res_names)]
  
  ko_rescue_short_lfc <- lapply(dynamic_res, function(x) {
    results(ko_rescue_short_dds, 
            name = x) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("firre_inducedfirre_induced.timepoint", "", result_name)))
  }) %>% bind_rows()
  
  # Shrunken LFC results
  ko_rescue_short_shrnklfc <- lapply(dynamic_res, function(x) {
    lfcShrink(ko_rescue_short_dds, 
              coef = x,
              type = "apeglm") %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("firre_inducedfirre_induced.timepoint", "", result_name)))
  }) %>% bind_rows()
  
  # Calculate the maximum fold-change in any one timepoint
  ko_rescue_short_maxfc <- ko_rescue_short_shrnklfc %>%
    group_by(gene_id) %>%
    summarize(max_fc = max(abs(log2FoldChange))) 
  
  ko_rescue_short_shrnklfc <- ko_rescue_short_shrnklfc %>%
    left_join(ko_rescue_short_maxfc)
  
  save(ko_rescue_short_lfc, ko_rescue_short_shrnklfc, file = "results/ko_rescue_short.RData")
}

load("results/ko_rescue_short.RData")
```

### Fold changes vs zero timepoint

This is without considering the control cell line.

``` r
if(!file.exists("results/ko_rescue_short_vs_zero.RData")) {
  
  ko_rescue_short_vszero_samples <- samples %>%
    filter(cell_type == "ESC",
           timecourse_length == "short",
           firre_ko == "KO",
           firre_induced == "firre_induced")
  ko_rescue_short_vszero_counts <- salmon_gene_counts[,ko_rescue_short_vszero_samples$sample_id]
  
  # Check ordering
  stopifnot(all(rownames(ko_rescue_short_vszero_samples) == colnames(ko_rescue_short_vszero_counts)))
  stopifnot(all(rownames(ko_rescue_short_vszero_counts) == genes$gene_id))
  
  # DESeq2 -- controlling for doxycycline; likelihood ratio test
  ko_rescue_short_vszero_dds <- DESeqDataSetFromMatrix(countData = ko_rescue_short_vszero_counts,
                                                       colData = ko_rescue_short_vszero_samples,
                                                       design = ~ timepoint)
  ko_rescue_short_vszero_dds <- DESeq(ko_rescue_short_vszero_dds)
  res_names <- resultsNames(ko_rescue_short_vszero_dds)
  
  vs_zero_res <- res_names[grepl("_vs_0", res_names)]
  ko_rescue_short_vszero_shrnklfc <- lapply(vs_zero_res, function(x) {
    lfcShrink(ko_rescue_short_vszero_dds, 
              coef = x,
              type = "apeglm") %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("timepoint_|_vs_0", "", result_name)))
  }) %>% bind_rows()
  
  save(ko_rescue_short_vszero_shrnklfc, file = "results/ko_rescue_short_vs_zero.RData")
}

load("results/ko_rescue_short_vs_zero.RData")
```

### Rescued genes

These are the genes that are reciprocally regulated in the KO vs WT and
the Firre induced line.

``` r
# We need the KO vs WT data
load("../04_firre_ko_vs_wt/results/wtko.RData")

wtko_res_shrnklfc <- wtko_res_shrnklfc %>%
  mutate(sig = padj <= pval_thresh & abs(log2FoldChange) > l2fc_thresh)

wtko_sig <- wtko_res_shrnklfc %>% 
  filter(sig == TRUE)

# Let's figure out which genes are rescued
rescued_genes <- wtko_sig$gene_id[wtko_sig$gene_id %in% 
                                    (ko_rescue_short_shrnklfc %>%
                                       filter(padj <= pval_thresh) %>% 
                                       pull(gene_id))]

rescue_wtko <- wtko_sig %>%
  filter(gene_id %in% rescued_genes) %>%
  dplyr::select(gene_id, gene_name, log2FoldChange) %>%
  dplyr::rename(l2fc_wtko = log2FoldChange)

rescue_ltc <- ko_rescue_short_vszero_shrnklfc %>%
  filter(gene_id %in% rescued_genes) %>%
  dplyr::select(gene_id, gene_name, log2FoldChange, timepoint) %>%
  dplyr::rename(l2fc_ltc = log2FoldChange)

rescued_df <- rescue_wtko %>%
  left_join(rescue_ltc)
```

    ## Joining, by = c("gene_id", "gene_name")

``` r
rescued_df <- rescued_df %>%
  mutate(dist = abs((l2fc_wtko * -1) - l2fc_ltc))

rescued_df_closest <- rescued_df %>%
  group_by(gene_id, gene_name) %>%
  mutate(closest = min(dist) == dist,
         max_fc = max(abs(l2fc_ltc)) == abs(l2fc_ltc)) %>%
  filter(closest == TRUE)

rescued_df_closest <-  rescued_df_closest %>%
  mutate(twenty_percent = l2fc_wtko * 0.2,
         rescued = (-1*l2fc_wtko - twenty_percent)>= l2fc_ltc & l2fc_ltc>= (-1*l2fc_wtko + twenty_percent) |
           (-1*l2fc_wtko + twenty_percent)>= l2fc_ltc & l2fc_ltc>= (-1*l2fc_wtko - twenty_percent))

rescued_genes_20 <- rescued_df_closest %>% 
  filter(rescued == TRUE) %>% 
  pull(gene_name)

ggplot(rescued_df_closest, aes(x = l2fc_wtko, y = l2fc_ltc, color = rescued)) +
  geom_abline(slope = -1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = -0.8, lty = 2) +
  geom_abline(slope = -1.2, lty = 2) +
  xlim(-2,2) +
  ylim(-2,2)
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](short_timecourse_rescue_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave(paste0("figures/short_rescue_lfc_scatter_", thresh, ".pdf"), 
       height = 3, width = 3, useDingbats = FALSE)
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

### Short timecourse call significant genes

We’ll make the p-value cutoff based on the dox controlled model and the
l2fc cutoff based on the fold change vs zero.

``` r
ko_rescue_short_dox_sig <- ko_rescue_short_shrnklfc %>% 
  filter(padj <= pval_thresh)

ko_rescue_short_vszero_sig <- ko_rescue_short_vszero_shrnklfc %>%
  filter(gene_id %in% ko_rescue_short_dox_sig$gene_id)

ko_rescue_short_vszero_maxfc <- ko_rescue_short_vszero_sig %>%
  group_by(gene_id) %>%
  summarize(max_fc = max(abs(log2FoldChange))) 

ko_rescue_short_vszero_sig <- ko_rescue_short_vszero_sig %>%
  left_join(ko_rescue_short_vszero_maxfc)
```

    ## Joining, by = "gene_id"

``` r
ko_rescue_short_vszero_sig <- ko_rescue_short_vszero_sig %>%
  filter(max_fc > l2fc_thresh)

save(ko_rescue_short_vszero_sig, file = "results/ko_rescue_short_vszero_sig.RData")
```

### Short timecourse Firre responders heatmap

``` r
# Heatmap of fold-changes for DEGs in the rescue
# Check that there are no duplicate row names.
stopifnot(all(length(unique(ko_rescue_short_vszero_sig$gene_id)) == length(unique(ko_rescue_short_vszero_sig$gene_name))))

ko_rescue_short_lfc <- ko_rescue_short_vszero_sig %>%
  dplyr::select(gene_name, timepoint, log2FoldChange) %>%
  pivot_wider(names_from = timepoint, names_sort = TRUE, values_from = log2FoldChange) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

# Add a zero column.
ko_rescue_short_lfc <- cbind(matrix(0, nrow = nrow(ko_rescue_short_lfc), ncol = 1), ko_rescue_short_lfc)
colnames(ko_rescue_short_lfc)[[1]] <- "0"

row_ha = rowAnnotation(sig_in_ko = as.numeric(rownames(ko_rescue_short_lfc) %in% wtko_sig$gene_name),
                       rescued = as.numeric(rownames(ko_rescue_short_lfc) %in% rescued_genes_20),
                       col = list(sig_in_ko = c("1" = "black", "0" = "white"),
                                  rescued = c("1" = "black", "0" = "white")))

table(as.numeric(rownames(ko_rescue_short_lfc) %in% wtko_sig$gene_name))
```

    ## 
    ##  0  1 
    ## 14  5

``` r
table(as.numeric(rownames(ko_rescue_short_lfc) %in% rescued_genes_20))
```

    ## 
    ##  0  1 
    ## 17  2

``` r
rownames(ko_rescue_short_lfc)[rownames(ko_rescue_short_lfc) %in% rescued_genes_20]
```

    ## [1] "Slc16a3" "Firre"

``` r
pdf(paste0("figures/short_firre_responders_heatmap_", thresh, ".pdf"), 
    width = 4, height = 3.5)
h1 <- Heatmap(ko_rescue_short_lfc, 
        name = "l2fc",
        cluster_columns = FALSE, show_row_names = FALSE, 
        col = colorRamp2(seq(-2,2,length.out = 100), col_pal10),
        left_annotation = row_ha)
draw(h1)
dev.off()
```

    ## png 
    ##   2

``` r
draw(h1)
```

![](short_timecourse_rescue_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
