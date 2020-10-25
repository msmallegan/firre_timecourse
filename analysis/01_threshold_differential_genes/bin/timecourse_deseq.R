#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

options(stringsAsFactors = FALSE)
library(tidyverse)
library(DESeq2)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  f <- args[1]
} else {
  stop("Too many args")
}

f <- paste0("../", f)
exp <- gsub(".csv", "",gsub("../data//","",f))

samples <- read_csv("../../../rnaseq/rnaseq_samplesheet.csv") %>%
  unite(experiment, cell_type, firre_ko, firre_induced, timecourse_length, 
        experiment_replicate, sep = "-",
        remove = FALSE) %>%
  filter(experiment == exp) %>%
  mutate(rn = sample_id) %>%
  column_to_rownames("rn")

t <- samples$timepoint %>%
  unique() %>%
  sort()


# Factorfy the timepoint for deseq
samples$timepoint <- factor(samples$timepoint, 
                                    levels = t)


counts <- read_csv(f) %>%
  dplyr::select(gene_id, sample_id, count) %>%
  pivot_wider(names_from = "sample_id", values_from = "count") %>%
  column_to_rownames("gene_id") %>%
  as.matrix()


# Check ordering
stopifnot(all(rownames(samples) %in% colnames(counts)))
counts <- counts[,rownames(samples)]
mode(counts) <- "integer"
stopifnot(all(rownames(samples) == colnames(counts)))

# DEseq
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = samples, 
                              design = ~ timepoint)
# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10,]
# Run deseq -- individual timepoints
dds <- DESeq(dds)

res_wrapper <- function(x) {
  contrast <- unlist(str_split(gsub("_vs", "", x), "_"))
  res <- results(dds, contrast = contrast) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")
  return(res)
}

shrink_wrapper <- function(x) {
  contrast <- unlist(str_split(gsub("_vs", "", x), "_"))
  res <- lfcShrink(dds, contrast = contrast, type = "ashr") %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")
  return(res)
}

# Retrieve single timepoint comparisons
# st -- single_timepoint
# tc -- timecourse
st <- tibble(test = c(t[-1], t[3:length(t)]),
                      control = c(rep(t[1], length(t[-1])), t[2:(length(t)-1)])) %>%
  mutate(result_name = map2(test, control, ~ paste0("timepoint_", .x, "_vs_", .y))) %>%
  unnest(result_name) %>%
  mutate(unshrunken = map(result_name, res_wrapper),
         shrunken = map(result_name, shrink_wrapper)) %>%
  unnest(c(unshrunken, shrunken), names_sep = "_") %>%
  dplyr::rename(gene_id = unshrunken_gene_id,
                unshrunken = unshrunken_log2FoldChange,
                shrunken = shrunken_log2FoldChange,
                pvalue = unshrunken_pvalue,
                padj = unshrunken_padj) %>%
  dplyr::select(result_name, gene_id, unshrunken, 
                shrunken, pvalue, padj) %>%
  pivot_longer(cols = c(unshrunken, shrunken), 
               names_to = "shrunken", values_to = "l2fc")


# Run deseq -- timecourse model
dds <- DESeq(dds, test="LRT", reduced = ~1)

tc <- data.frame("result_name" = resultsNames(dds)[-1]) %>%
  mutate(unshrunken = map(result_name, res_wrapper),
         shrunken = map(result_name, shrink_wrapper)) %>%
  unnest(c(unshrunken, shrunken), names_sep = "_") %>%
  dplyr::rename(gene_id = unshrunken_gene_id,
                pvalue = unshrunken_pvalue,
                padj = unshrunken_padj) %>%
  dplyr::select(gene_id, pvalue, padj) %>%
  distinct() %>%
  mutate("result_name" = "lrt_tc") %>%
  dplyr::select(result_name, everything()) %>%
  mutate(shrunken = NA, l2fc = NA)


# Combine results
res <- bind_rows(st, tc)
write_csv(res, paste0("../results/", exp, ".csv"))


# Get significantly DEG genes
pval_thresh <- 0.05

is_consecutive <- function(tps) {
  if(length(tps) == 1) {
    return(FALSE)
  }
  tps <- tps[order(tps)]
  idx <- which(t %in% tps)
  # Will take only two consecutive timepoints
  if(min(diff(idx)) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


vs_zero <- res %>% filter(grepl("vs_0", result_name),
                          padj < pval_thresh,
                          shrunken == "shrunken") %>%
  mutate(timepoint = as.numeric(str_split(result_name, "_", simplify = TRUE)[, 2])) %>%
  group_by(gene_id) %>%
  mutate(sig_count = n()) %>%
  filter(sig_count > 1 | timepoint == max(t)) %>%
  group_by(gene_id) %>%
  mutate(consecutive = is_consecutive(timepoint)) %>%
  filter(consecutive == TRUE | timepoint == max(t)) %>%
  mutate(sig = TRUE,
         result_name = gsub("timepoint_", "t_", result_name)) %>%
  dplyr::select(gene_id, result_name, sig)

vs_prev <- res %>% filter(result_name != "lrt_tc", 
                          !grepl("vs_0", result_name),
                          padj < pval_thresh,
                          shrunken == "shrunken") %>%
  mutate(timepoint = as.numeric(str_split(result_name, "_", simplify = TRUE)[, 2])) %>%
  group_by(gene_id) %>%
  mutate(sig_count = n()) %>%
  filter(sig_count > 1) %>%
  mutate(consecutive = is_consecutive(timepoint)) %>%
  filter(consecutive == TRUE) %>%
  mutate(sig = TRUE,
         result_name = gsub("timepoint_", "t_", result_name)) %>%
  dplyr::select(gene_id, result_name, sig)

tc <- res %>% filter(result_name == "lrt_tc",
                     padj < pval_thresh) %>%
  mutate(sig = TRUE) %>%
  dplyr::select(gene_id, result_name, sig)

sig_matrix <- bind_rows(vs_zero, vs_prev, tc) 
write_csv(sig_matrix, paste0("../results/", exp, "_sig_matrix.csv"))

res <- res %>% filter(gene_id %in% unique(sig_matrix$gene_id))
write_csv(sig_matrix, paste0("../results/", exp, "_sig.csv"))

