Long timecourse rescue
================

For this experiment, we induced the expression of Firre with an rTTA
element by adding doxycycline to the cells. We see that Firre is indeed
expressed in the KO background after the addition of doxycycline. The
drug does instigate some gene expression changes on its own, so we will
control for the effects by using a linear model which accounts for the
effect of dox.

``` r
if(!file.exists("results/ko_rescue_long.RData")) {
  
  # Filter to ESC KO long timecourse
  ko_rescue_long_samples <- samples %>%
    filter(cell_type == "ESC",
           timecourse_length == "long",
           firre_ko == "KO")
  ko_rescue_long_counts <- salmon_gene_counts[,ko_rescue_long_samples$sample_id]
  
  # Check ordering
  stopifnot(all(rownames(ko_rescue_long_samples) == colnames(ko_rescue_long_counts)))
  stopifnot(all(rownames(ko_rescue_long_counts) == genes$gene_id))
  
  # DESeq2 -- controlling for doxycycline; likelihood ratio test
  ko_rescue_long_dds <- DESeqDataSetFromMatrix(countData = ko_rescue_long_counts, 
                                          colData = ko_rescue_long_samples, 
                                          design = ~ firre_induced + timepoint + timepoint*firre_induced)
  ko_rescue_long_dds <- DESeq(ko_rescue_long_dds, test = "LRT", reduced = ~ firre_induced + timepoint)
  
  # Compile results
  res_names <- resultsNames(ko_rescue_long_dds)
  dynamic_res <- res_names[grepl("firre_inducedfirre_induced.timepoint", res_names)]
  
  ko_rescue_long_lfc <- lapply(dynamic_res, function(x) {
    results(ko_rescue_long_dds, 
            name = x) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("firre_inducedfirre_induced.timepoint", "", result_name)))
  }) %>% bind_rows()
  
  # Shrunken LFC results
  ko_rescue_long_shrnklfc <- lapply(dynamic_res, function(x) {
    lfcShrink(ko_rescue_long_dds, 
              coef = x,
              type = "apeglm") %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("firre_inducedfirre_induced.timepoint", "", result_name)))
  }) %>% bind_rows()
  
  # Calculate the maximum fold-change in any one timepoint
  ko_rescue_long_maxfc <- ko_rescue_long_shrnklfc %>%
    group_by(gene_id) %>%
    summarize(max_fc = max(abs(log2FoldChange))) 
  
  ko_rescue_long_shrnklfc <- ko_rescue_long_shrnklfc %>%
    left_join(ko_rescue_long_maxfc)
  
  save(ko_rescue_long_lfc, ko_rescue_long_shrnklfc, file = "results/ko_rescue_long_long.RData")
}

load("results/ko_rescue_long.RData")
```

### Fold changes vs zero timepoint

This is without considering the control cell line.

``` r
if(!file.exists("results/ko_rescue_long_vs_zero.RData")) {
  
  ko_rescue_long_vszero_samples <- samples %>%
    filter(cell_type == "ESC",
           timecourse_length == "long",
           firre_ko == "KO",
           firre_induced == "firre_induced")
  ko_rescue_long_vszero_counts <- salmon_gene_counts[,ko_rescue_long_vszero_samples$sample_id]
  
  # Check ordering
  stopifnot(all(rownames(ko_rescue_long_vszero_samples) == colnames(ko_rescue_long_vszero_counts)))
  stopifnot(all(rownames(ko_rescue_long_vszero_counts) == genes$gene_id))
  
  # DESeq2 -- controlling for doxycycline; likelihood ratio test
  ko_rescue_long_vszero_dds <- DESeqDataSetFromMatrix(countData = ko_rescue_long_vszero_counts,
                                                 colData = ko_rescue_long_vszero_samples,
                                                 design = ~ timepoint)
  ko_rescue_long_vszero_dds <- DESeq(ko_rescue_long_vszero_dds)
  res_names <- resultsNames(ko_rescue_long_vszero_dds)
  
  vs_zero_res <- res_names[grepl("_vs_0", res_names)]
  ko_rescue_long_vszero_shrnklfc <- lapply(vs_zero_res, function(x) {
    lfcShrink(ko_rescue_long_vszero_dds, 
              coef = x,
              type = "apeglm") %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("timepoint_|_vs_0", "", result_name)))
  }) %>% bind_rows()
  
  save(ko_rescue_long_vszero_shrnklfc, file = "results/ko_rescue_long_vs_zero.RData")
}

load("results/ko_rescue_long_vs_zero.RData", verbose = T)
```

    ## Loading objects:
    ##   ko_rescue_long_vszero_shrnklfc

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
                                    (ko_rescue_long_shrnklfc %>%
                                       filter(padj <= pval_thresh) %>% 
                                       pull(gene_id))]

rescue_wtko <- wtko_sig %>%
  filter(gene_id %in% rescued_genes) %>%
  dplyr::select(gene_id, gene_name, log2FoldChange) %>%
  dplyr::rename(l2fc_wtko = log2FoldChange)

rescue_ltc <- ko_rescue_long_vszero_shrnklfc %>%
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

![](long_timecourse_rescue_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave(paste0("figures/long_rescue_lfc_scatter_", thresh, ".pdf"), 
       height = 3, width = 3, useDingbats = FALSE)
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

### Long timecourse call significant genes

We’ll make the p-value cutoff based on the dox controlled model and the
l2fc cutoff based on the fold change vs zero.

``` r
ko_rescue_long_dox_sig <- ko_rescue_long_shrnklfc %>% 
  filter(padj <= pval_thresh)

ko_rescue_long_vszero_sig <- ko_rescue_long_vszero_shrnklfc %>%
  filter(gene_id %in% ko_rescue_long_dox_sig$gene_id)

ko_rescue_long_vszero_maxfc <- ko_rescue_long_vszero_sig %>%
  group_by(gene_id) %>%
  summarize(max_fc = max(abs(log2FoldChange))) 

ko_rescue_long_vszero_sig <- ko_rescue_long_vszero_sig %>%
  left_join(ko_rescue_long_vszero_maxfc)
```

    ## Joining, by = "gene_id"

``` r
ko_rescue_long_vszero_sig <- ko_rescue_long_vszero_sig %>%
  filter(max_fc > l2fc_thresh)

save(ko_rescue_long_vszero_sig, file = "results/ko_rescue_long_vszero_sig.RData")
```

### Long timecourse Firre responders heatmap

``` r
# Heatmap of fold-changes for DEGs in the rescue
# Check that there are no duplicate row names.
stopifnot(all(length(unique(ko_rescue_long_vszero_sig$gene_id)) == length(unique(ko_rescue_long_vszero_sig$gene_name))))

ko_rescue_long_lfc <- ko_rescue_long_vszero_sig %>%
  dplyr::select(gene_name, timepoint, log2FoldChange) %>%
  pivot_wider(names_from = timepoint, names_sort = TRUE, values_from = log2FoldChange) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

# Add a zero column.
ko_rescue_long_lfc <- cbind(matrix(0, nrow = nrow(ko_rescue_long_lfc), ncol = 1), ko_rescue_long_lfc)
colnames(ko_rescue_long_lfc)[[1]] <- "0"

row_ha = rowAnnotation(sig_in_ko = as.numeric(rownames(ko_rescue_long_lfc) %in% wtko_sig$gene_name),
                       rescued = as.numeric(rownames(ko_rescue_long_lfc) %in% rescued_genes_20),
                       col = list(sig_in_ko = c("1" = "black", "0" = "white"),
                                  rescued = c("1" = "black", "0" = "white")))

table(as.numeric(rownames(ko_rescue_long_lfc) %in% wtko_sig$gene_name))
```

    ## 
    ##   0   1 
    ## 313 102

``` r
table(as.numeric(rownames(ko_rescue_long_lfc) %in% rescued_genes_20))
```

    ## 
    ##   0   1 
    ## 348  67

``` r
pdf(paste0("figures/long_responders_heatmap_", thresh, ".pdf"), width = 3, height = 5)
h1 <- Heatmap(ko_rescue_long_lfc, 
        name = "l2fc",
        cluster_columns = FALSE, show_row_names = FALSE, 
        col = colorRamp2(seq(-3,3,length.out = 100), col_pal10),
        left_annotation = row_ha)
draw(h1)
dev.off()
```

    ## png 
    ##   2

``` r
draw(h1)
```

![](long_timecourse_rescue_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Firre rescued genes heatmap

We can visualize the expression levels of the genes that were
significantly differentially expressed in the Firre knockout and
dynamically changing in response to Firre’s induction in a knockout
background.

``` r
# Retrieve the set of genes that is DE in the KO and is significantly changing
# in the KO rescue experiment. 
# Retrieve the TPM -- take the mean TPM in each condition.
rescued_genes_20 <- rescued_df_closest %>% 
  filter(rescued == TRUE) %>% 
  pull(gene_id)

wtko_samples <- samples %>% 
    filter(cell_type == "ESC", 
           timepoint == 0, 
           timecourse_length == "long")

wtko_matrix <- tpm %>% 
  filter(gene_id %in% rescued_genes_20) %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(wtko_samples) %>%
  drop_na() %>%
  group_by(gene_name, firre_ko) %>%
  summarize(tpm = mean(tpm)) %>%
  pivot_wider(names_from = "firre_ko", 
              values_from = "tpm")
```

    ## Joining, by = "sample_id"

``` r
ko_rescue_long_samples <- samples %>%
    filter(cell_type == "ESC",
           timecourse_length == "long",
           firre_ko == "KO")

ko_rescue_matrix <- tpm %>% 
  filter(gene_id %in% rescued_genes_20) %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(ko_rescue_long_samples) %>%
  drop_na() %>%
  filter(timepoint != 0) %>%
  group_by(gene_name, timepoint) %>%
  summarize(tpm = mean(tpm)) %>%
  pivot_wider(names_from = "timepoint", 
              values_from = "tpm")
```

    ## Joining, by = "sample_id"

``` r
comb_matrix <- wtko_matrix %>% 
  left_join(ko_rescue_matrix) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
```

    ## Joining, by = "gene_name"

``` r
comb_matrix <- comb_matrix[,c("KO", "720", "1440", "2880", "5760", "WT")]

comb_matrix_scaled <- t(scale(t(comb_matrix)))

p1 <- pheatmap::pheatmap(comb_matrix_scaled, cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE,
                   color = col_pal10,
                   clustering_callback = olo_seriate,
                   clustering_distance_rows = "euclidean",
                   border_color = NA,
                   fontsize_row = 8,
                   breaks = seq(-1.5, 1.5, length.out = length(col_pal10)),
                   treeheight_row = 25)
```

![](long_timecourse_rescue_files/figure-gfm/rescue_heatmap-1.png)<!-- -->

``` r
ggsave(paste0("figures/long_rescue_heatmap_", thresh, ".pdf"), 
       plot = p1, width = 3.5, height = 5)
```

### Rescue PCA

``` r
pca_dat <- prcomp(t(comb_matrix))
proportion_of_variance <- summary(pca_dat)$importance[2,1:2]

pca_df <- pca_dat$x[,1:2] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(PC1 = PC1 * -1)

ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
  geom_point() +
  geom_text_repel() +
  xlab(paste0("PC1: ", round(proportion_of_variance["PC1"]*100), "%")) +
  ylab(paste0("PC2: ", round(proportion_of_variance["PC2"]*100), "%"))
```

![](long_timecourse_rescue_files/figure-gfm/rescue_pca-1.png)<!-- -->

``` r
ggsave(paste0("figures/long_rescue_pca_", thresh, ".pdf"), 
       width = 3.5, height = 0.9, 
       useDingbats = FALSE)
```
