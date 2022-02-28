library(tidyverse)
library(DESeq2)

deseq_vs_zero <- function(ct, tc_len, fko, fi) {
  
  load("../01_setup/results/rnaseq_data.RData")
  
  samples <- samples %>%
    filter(cell_type == ct,
           timecourse_length == tc_len,
           firre_ko == fko,
           firre_induced == fi)
  counts <- salmon_gene_counts[,samples$sample_id]
  
  # Check ordering
  stopifnot(all(rownames(samples) == colnames(counts)))
  stopifnot(all(rownames(counts) == genes$gene_id))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = samples, 
                                design = ~ timepoint)
  dds <- DESeq(dds)
  
  # Compile results
  res_names <- resultsNames(dds)
  dynamic_res <- res_names[grepl("_vs_0", res_names)]
  
  lfc <- lapply(dynamic_res, function(x) {
    results(dds, 
            name = x) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("timepoint_|_vs_0", "", result_name)))
  }) %>% bind_rows()
  
  # Shrunken LFC results
  shrnklfc <- lapply(dynamic_res, function(x) {
    lfcShrink(dds, 
              coef = x,
              type = "apeglm") %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>% 
      merge(g2s) %>%
      mutate(result_name = x,
             timepoint = as.numeric(gsub("timepoint_|_vs_0", "", result_name)))
  }) %>% bind_rows()
  
  res_df <- lfc %>%
    left_join(shrnklfc %>%
                dplyr::rename(l2fc_shrunken = log2FoldChange) %>%
                dplyr::select(gene_id, result_name, timepoint, l2fc_shrunken))
  
  # Calculate the maximum fold-change in any one timepoint
  maxfc <- res_df %>%
    group_by(gene_id) %>%
    summarize(max_abs_l2fc_shrunken = max(abs(l2fc_shrunken)))
  
  res_df <- res_df %>%
    left_join(maxfc)
  
  return(res_df)
}