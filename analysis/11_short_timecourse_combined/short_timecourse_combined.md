Short timecourse combined
================

For this experiment, we induced the expression of Firre with an rTTA
element by adding doxycycline to the cells. We see that Firre is indeed
expressed in the KO background after the addition of doxycycline. The
drug does instigate some gene expression changes on its own, so we will
control for the effects by using a linear model which accounts for the
effect of dox.

``` r
if(!file.exists("results/short.RData")) {
  
  # Filter to ESC KO long timecourse
  short_samples <- samples %>%
    filter(cell_type == "ESC",
           timecourse_length == "short")
  short_counts <- salmon_gene_counts[,short_samples$sample_id]
  
  # Check ordering
  stopifnot(all(rownames(short_samples) == colnames(short_counts)))
  stopifnot(all(rownames(short_counts) == genes$gene_id))
  
  # DESeq2 -- controlling for doxycycline; likelihood ratio test
  short_dds <- DESeqDataSetFromMatrix(countData = short_counts, 
                                      colData = short_samples, 
                                      design = ~ firre_ko + firre_induced + timepoint + timepoint*firre_induced)
  short_dds <- DESeq(short_dds, test = "LRT", reduced = ~ firre_ko + firre_induced + timepoint)
  
  # Compile results
  res_names <- resultsNames(short_dds)
  dynamic_res <- res_names[grepl("firre_inducedfirre_induced.timepoint", res_names)]
  
  short_lfc <- lapply(dynamic_res, function(x) {
    results(short_dds, 
            name = x) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("firre_inducedfirre_induced.timepoint", "", result_name)))
  }) %>% bind_rows()
  
  # Shrunken LFC results
  short_shrnklfc <- lapply(dynamic_res, function(x) {
    lfcShrink(short_dds, 
              coef = x,
              type = "apeglm") %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("firre_inducedfirre_induced.timepoint", "", result_name)))
  }) %>% bind_rows()
  
  # Calculate the maximum fold-change in any one timepoint
  short_maxfc <- short_shrnklfc %>%
    group_by(gene_id) %>%
    summarize(max_fc = max(abs(log2FoldChange))) 
  
  short_shrnklfc <- short_shrnklfc %>%
    left_join(short_maxfc)
  
  save(short_lfc, short_shrnklfc, file = "results/short.RData")
}

load("results/short.RData")
```

### Fold changes vs zero timepoint

This is without considering the control cell line.

``` r
if(!file.exists("results/short_vs_zero.RData")) {
  
  short_vszero_samples <- samples %>%
    filter(cell_type == "ESC",
           timecourse_length == "short",
           firre_induced == "firre_induced")
  short_vszero_counts <- salmon_gene_counts[,short_vszero_samples$sample_id]
  
  # Check ordering
  stopifnot(all(rownames(short_vszero_samples) == colnames(short_vszero_counts)))
  stopifnot(all(rownames(short_vszero_counts) == genes$gene_id))
  
  # DESeq2 -- controlling for doxycycline; likelihood ratio test
  short_vszero_dds <- DESeqDataSetFromMatrix(countData = short_vszero_counts,
                                             colData = short_vszero_samples,
                                             design = ~ firre_ko + timepoint)
  short_vszero_dds <- DESeq(short_vszero_dds)
  res_names <- resultsNames(short_vszero_dds)
  
  vs_zero_res <- res_names[grepl("_vs_0", res_names)]
  short_vszero_shrnklfc <- lapply(vs_zero_res, function(x) {
    lfcShrink(short_vszero_dds, 
              coef = x,
              type = "apeglm") %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("timepoint_|_vs_0", "", result_name)))
  }) %>% bind_rows()
  # short_vszero_shrnklfc_notadjust_for_firreko <- short_vszero_shrnklfc 
  # save(short_vszero_shrnklfc_notadjust_for_firreko, file = "results/short_vs_zero_notadjust_for_firreko.RData")
  save(short_vszero_shrnklfc, file = "results/short_vs_zero.RData")
}

load("results/short_vs_zero.RData", verbose = T)
```

    ## Loading objects:
    ##   short_vszero_shrnklfc

``` r
# load("results/short_vs_zero_notadjust_for_firreko.RData", verbose = T)
# short_vszero_shrnklfc <- short_vszero_shrnklfc_notadjust_for_firreko
```

### Short timecourse call significant genes

We’ll make the p-value cutoff based on the dox controlled model and the
l2fc cutoff based on the fold change vs zero.

``` r
short_dox_sig <- short_shrnklfc %>% 
  filter(padj <= pval_thresh)

short_vszero_sig <- short_vszero_shrnklfc %>%
  filter(gene_id %in% short_dox_sig$gene_id)

short_vszero_maxfc <- short_vszero_sig %>%
  group_by(gene_id) %>%
  summarize(max_fc = max(abs(log2FoldChange))) 

short_vszero_sig <- short_vszero_sig %>%
  left_join(short_vszero_maxfc)
```

    ## Joining, by = "gene_id"

``` r
short_vszero_sig <- short_vszero_sig %>%
  filter(max_fc > l2fc_thresh)

save(short_vszero_sig, file = "results/short_vszero_sig.RData")
```

### Short timecourse Firre responders heatmap

``` r
# Let's look at the set of genes that overlap in the two genetic backgrounds
load("../06_short_timecourse_rescue/results/ko_rescue_short_vszero_sig.RData", verbose = T)
```

    ## Loading objects:
    ##   ko_rescue_short_vszero_sig

``` r
load("../09_short_timecourse_overexpression/results/wt_overexp_short_vszero_sig.RData", verbose = T)
```

    ## Loading objects:
    ##   wt_overexp_short_vszero_sig

``` r
overlapping_genes <- unique(ko_rescue_short_vszero_sig$gene_id)[unique(ko_rescue_short_vszero_sig$gene_id) %in% unique(wt_overexp_short_vszero_sig$gene_id)]

# Heatmap of fold-changes for DEGs in the rescue
# Check that there are no duplicate row names.
stopifnot(all(length(unique(short_vszero_sig$gene_id)) == length(unique(short_vszero_sig$gene_name))))

short_lfc <- short_vszero_sig %>%
  dplyr::select(gene_name, timepoint, log2FoldChange) %>%
  pivot_wider(names_from = timepoint, names_sort = TRUE, values_from = log2FoldChange) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()


ordering <- short_vszero_sig %>%
  filter(abs(log2FoldChange) > 0.2) %>%
  group_by(gene_name) %>%
  summarize(first_tp_de = min(timepoint),
            max_fc = max(log2FoldChange)) %>%
  arrange(first_tp_de,
          -max_fc)

# Add a zero column.
short_lfc <- cbind(matrix(0, nrow = nrow(short_lfc), ncol = 1), short_lfc)
colnames(short_lfc)[[1]] <- "0"

short_lfc <- short_lfc[ordering$gene_name,]

pdf(paste0("figures/combined_responders_heatmap_", thresh, ".pdf"), 
    width = 4, height = 3.5)
ht1 <- Heatmap(short_lfc, 
               name = "l2fc",
               cluster_columns = FALSE, show_row_names = TRUE, 
               cluster_rows = FALSE,
               col = colorRamp2(seq(-4,4,length.out = 100), col_pal10))
draw(ht1)
dev.off()
```

    ## png 
    ##   2

``` r
draw(ht1)
```

![](short_timecourse_combined_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
