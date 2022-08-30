run_timecourse_deseq <- function(experiment,
                                 counts, samples, genes,
                                 design_formula,  
                                 reduced_formula, 
                                 independent_filtering = TRUE,
                                 ncores = 12,
                                 save_dds = TRUE,
                                 save_res = TRUE) {
  
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
  stopifnot(all(names(genes) == rownames(counts)))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = samples, 
                                design = design_formula,
                                rowData = genes)
  
  # Filter low counts
  dds <- dds[rowSums(counts(dds)) >= 100,]
  
  
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
  
  # Compile results
  res_names <- resultsNames(dds)
  res <- results(dds, 
                 name = res_names[2],
                 independentFiltering=independent_filtering) 
  
  res_shrunken <- lfcShrink(dds = dds,
                            coef = res_names[2],
                            res = res,
                            type = "apeglm",
                            parallel = TRUE,
                            BPPARAM=MulticoreParam(ncores)) %>% 
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    merge(g2s) %>%
    mutate(result_name = res_names[2]) 
  
  if(save_res == TRUE) {
    saveRDS(res, paste0("results/",
                        ct, "_",
                        fko, "_res.rds"))
  }
  
  # Convert to data frame
  res <- res %>% 
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    merge(g2s) %>%
    mutate(result_name = res_names[2]) 
  
  for(i in 2:length(res_names)) {
    tmp_res <- results(dds, 
                       name = res_names[i],
                       independentFiltering=independent_filtering)
    
    tmp_res_shrunken <- suppressMessages(lfcShrink(dds = dds,
                                                   coef = res_names[i],
                                                   res = tmp_res,
                                                   type = "apeglm",
                                                   parallel = TRUE,
                                                   BPPARAM=MulticoreParam(ncores))) %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      merge(g2s) %>%
      mutate(result_name = res_names[i]) 
    
    tmp_res <- tmp_res %>% 
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      merge(g2s) %>%
      mutate(result_name = res_names[i]) 

    res_shrunken <- bind_rows(res_shrunken, tmp_res_shrunken)
    res <- bind_rows(res, tmp_res)
  }
  
  # Label results
  res_shrunken$cell_type <- ct
  res_shrunken$firre_ko <- fko
  res$cell_type <- ct
  res$firre_ko <- fko
  
  # Label the comparisons
  res$comparison <- "intercept"
  res[grep("timepoint_minutes_",res$result_name),"comparison"] <- "vs_zero"
  res[res$result_name == "firre_induced_firre_induced_vs_control",
        "comparison"] <- "static_firre_induction_vs_control"
  res[res$result_name == "firre_ko_KO_vs_WT",
        "comparison"] <- "KO_vs_WT"
  res[grep(".firre_induced",res$result_name,
             fixed = T),"comparison"] <- "dynamic_firre_induction_vs_control"
  
  combined_results <- list("res" = res, "res_shrunken" = res_shrunken)
  
  return(combined_results)
}

run_control_deseq <- function(experiment,
                                 counts, samples, genes,
                                 design_formula,  
                                 reduced_formula, 
                                 independent_filtering = TRUE,
                                 ncores = 12,
                                 save_dds = TRUE,
                                 save_res = TRUE) {
  
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
  stopifnot(all(names(genes) == rownames(counts)))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = samples, 
                                design = design_formula,
                                rowData = genes)
  
  # Filter low counts
  dds <- dds[rowSums(counts(dds)) >= 100,]
  
  
  dds <- DESeq(dds, 
               test="LRT", 
               reduced = reduced_formula,
               parallel=TRUE,
               BPPARAM=MulticoreParam(ncores))
  
  if(save_dds == TRUE) {
    saveRDS(dds, paste0("results/dox_control_",
                        ct, "_",
                        fko, "_dds.rds"))
  }
  
  # Compile results
  res_names <- resultsNames(dds)
  res <- results(dds, 
                 name = res_names[2],
                 independentFiltering=independent_filtering)
  
  res_shrunken <- lfcShrink(dds = dds,
                            coef = res_names[2],
                            res = res,
                            type = "apeglm",
                            parallel = TRUE,
                            BPPARAM=MulticoreParam(ncores)) %>% 
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    merge(g2s) %>%
    mutate(result_name = res_names[2]) 

  if(save_res == TRUE) {
    saveRDS(res, paste0("results/dox_control_",
                        ct, "_",
                        fko, "_res.rds"))
  }
  
  # Convert to data frame
  res <- res %>% 
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    merge(g2s) %>%
    mutate(result_name = res_names[2]) 
  
  for(i in 3:length(res_names)) {
    tmp_res <- results(dds, 
                       name = res_names[i],
                       independentFiltering=independent_filtering)
    
    tmp_res_shrunken <- suppressMessages(lfcShrink(dds = dds,
                              coef = res_names[i],
                              res = tmp_res,
                              type = "apeglm",
                              parallel = TRUE,
                              BPPARAM=MulticoreParam(ncores))) %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      merge(g2s) %>%
      mutate(result_name = res_names[i]) 
    
    tmp_res <- tmp_res %>% 
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      merge(g2s) %>%
      mutate(result_name = res_names[i]) 
    
    res_shrunken <- bind_rows(res_shrunken, tmp_res_shrunken)
    res <- bind_rows(res, tmp_res)
  }
  
  # Label results
  res_shrunken$cell_type <- ct
  res_shrunken$firre_ko <- fko
  res$cell_type <- ct
  res$firre_ko <- fko
  
  combined_results <- list("res" = res, "res_shrunken" = res_shrunken)
  
  return(combined_results)
}

get_dynamic_comparison_df <- function(res) {
  
  dyn_res <- res %>% 
    filter(comparison == "dynamic_firre_induction_vs_control")
  dyn_res$timepoint <- sapply(dyn_res$result_name,
                              function(x) {
                                x <- gsub("timepoint_minutes", "", x)
                                gsub(".firre_inducedfirre_induced", "", x,
                                     fixed = T)
                              }) %>%
    as.numeric()
  # Add zeros for visualization
  dyn_res_zeros <- dyn_res %>% dplyr::select(gene_id,
                                             baseMean,
                                             stat,
                                             pvalue,
                                             padj,
                                             gene_name,
                                             cell_type,
                                             firre_ko,
                                             comparison) %>%
    distinct()
  dyn_res_zeros$timepoint <- 0
  dyn_res_zeros$log2FoldChange <- 0
  dyn_res_zeros$lfcSE <- 0
  dyn_res_zeros$result_name <- "zeros"
  #Put in order for merging
  dyn_res_zeros <- dyn_res_zeros %>% dplyr::select(colnames(dyn_res))
  dyn_res <- bind_rows(dyn_res, dyn_res_zeros)
  dyn_res$timepoint <- factor(dyn_res$timepoint, 
                              levels = seq(0,360,30))
  return(dyn_res)
}


olo_seriate <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "euclidean")
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}

calculate_overlap_pval = function(list1, list2, total.size, lower.tail = FALSE, adjust = FALSE) {
  
  # calculate actual overlap
  actual.overlap <- length(intersect(list1, list2));
  
  # calculate expected overlap
  # need to cast to avoid integer overflow when length(list1) * length(list2) is extremely large
  expected.overlap <- as.numeric(length(list1)) * length(list2) / total.size;
  
  adjust.value <- 0;
  
  # adjust actual.overlap to reflect P[X >= x]
  if (adjust & !lower.tail) {
    adjust.value <- 1;
    warning('Calculating P[X >= x]');
  }
  
  # calculate significance of the overlap
  overlap.pvalue <- phyper(
    q = actual.overlap - adjust.value,
    m = length(list1),
    n = total.size - length(list1),
    k = length(list2),
    lower.tail = lower.tail
  );
  
  # return values
  return(overlap.pvalue);
  
}

# Credit: https://www.biostars.org/p/171766/
#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}



