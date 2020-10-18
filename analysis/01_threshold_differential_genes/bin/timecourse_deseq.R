options(stringsAsFactors = FALSE)
library(tidyverse)
f <- "../data/ESC-KO-control-long-1.csv"
exp <- gsub(".csv", "",gsub("../data/","",f))

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

# TODO: write method of determining whether it's a significant gene.
