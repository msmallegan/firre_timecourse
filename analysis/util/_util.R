run_timecourse_deseq <- function(experiment,
                                 counts, samples, genes,
                                 design_formula,  
                                 reduced_formula, 
                                 ncores = 12,
                                 save_dds = TRUE) {
  
  ct <- unlist(strsplit(experiment, "_"))[[1]]
  fko <- unlist(strsplit(experiment, "_"))[[2]]
  
  samples <- samples %>% filter(cell_type == ct,
                                firre_ko == fko)
  
  counts <- counts[,samples$sample_id]
  mode(counts) <- "integer"
  
  # Reorder gencode genes
  genes <- genes[rownames(counts)]
  
  # Check ordering
  stopifnot(all(rownames(samples) == colnames(counts)))
  stopifnot(all(names(gencode_genes) == rownames(counts)))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = samples, 
                                design = design_formula,
                                rowData = gencode_genes)
  dds <- DESeq(dds, 
               test="LRT", 
               reduced = reduced_formula,
               parallel=TRUE,
               BPPARAM=MulticoreParam(ncores))
  
  if(save_dds == TRUE) {
    saveRDS(dds, paste0("results/",
                        ct, "_",
                        fko, "_dds.rds"))
  }
  
  # TODO: Fix base mean filtering since I previously set
  # independent filtering to false.
  # https://support.bioconductor.org/p/76144/
  # res$pvalue[res$baseMean < x] <- NA
  # res$padj <- p.adjust(res$pvalue, method="BH")
  
  # Compile results
  res_names <- resultsNames(dds)
  res <- results(dds, 
                 name = res_names[1],
                 independentFiltering=FALSE,
                 cooksCutoff=FALSE) %>% 
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    merge(g2s) %>%
    mutate(result_name = res_names[1]) 
  
  for(i in 2:length(res_names)) {
    tmp_res <- results(dds, 
                       name = res_names[i],
                       independentFiltering=FALSE) %>% 
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      merge(g2s) %>%
      mutate(result_name = res_names[i]) 
    res <- bind_rows(res, tmp_res)
  }
  
  # Label results
  res$cell_type <- ct
  res$firre_ko <- fko
  return(res)
}