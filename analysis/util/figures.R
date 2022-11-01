library(tidyverse)

make_timecourse_lfc_plot <- function(deseq_res, genes_to_include, draw_plot = TRUE,
                                     x_adjust_scale = 20,
                                     y_lim = c(-5,5)) {
  # x_adjust_scale of 20 works for the long timecourse, 1 works for the short
  
  deg <- deseq_res %>%
    filter(gene_id %in% genes_to_include) %>%
    dplyr::select(gene_id, gene_name, baseMean, l2fc_shrunken, timepoint)
  
  # Add a zero timepoint
  deg_zero <- deg %>%
    dplyr::select(gene_id, gene_name, baseMean) %>%
    distinct() %>%
    mutate(l2fc_shrunken = 0,
           timepoint = 0) 
  
  deg <- deg %>%
    bind_rows(deg_zero)
  
  if(length(unique(deg$gene_id)) > 1) {
    deg <- deg %>%
      mutate(x_adjust = scale(log10(baseMean)),
             tp_adjust = timepoint + x_adjust*x_adjust_scale)
  } else {
    deg <- deg %>%
      mutate(tp_adjust = timepoint)
  }
  
  num_deg <- deg  %>%
    group_by(gene_id) %>%
    summarize(max_fc = l2fc_shrunken[which.max(abs(l2fc_shrunken))]) %>%
    mutate(direction = ifelse(max_fc > 0, "up", "down")) %>%
    group_by(direction) %>%
    summarize(num_deg = length(unique(gene_id)))
  
  label_x <- max(deg$timepoint / 2)
  
  g1 <- ggplot(deg, aes(x = tp_adjust, y = l2fc_shrunken, group = gene_id)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.1, color = "gray70") +
    geom_point(alpha = 0.4) +
    annotate("text", x = label_x, y = y_lim[[2]]/2, 
             label = paste0("n = ", num_deg %>% 
                              filter(direction == "up") %>%
                              pull(num_deg))) +
    annotate("text", x = label_x, y = y_lim[[1]]/2, 
             label = paste0("n = ", num_deg %>% 
                              filter(direction == "down") %>%
                              pull(num_deg))) +
    ylim(y_lim[[1]], y_lim[[2]])
  if(draw_plot == TRUE) {
    show(g1)
  }
  return(g1)
}
