library(bcbioRNASeq)
library(DESeq2)
library(DEGreport)
library(knitr)
library(tidyverse)

# Set seed for reproducibility
set.seed(1454944673L)

opts_chunk[["set"]](
  autodep = TRUE,
  bootstrap.show.code = FALSE,
  cache = TRUE,
  cache.lazy = TRUE,
  dev = c("png", "pdf"),
  fig.height = 10L,
  fig.retina = 2L,
  fig.width = 10L,
  highlight = TRUE,
  prompt = TRUE,
  tidy = FALSE
)

theme_paperwhite <- function(
  base_size = 14L,
  base_family = "",
  face = c("bold", "plain"),
  aspect_ratio = NULL,
  legend_position = c("right", "bottom", "top", "none"),
  grid = FALSE,
  minimal = FALSE
) {
  
  face <- match.arg(face)
  legend_position <- match.arg(legend_position)
  
  
  gray <- "gray95"
  
  text <- element_text(
    family = base_family,
    face = face,
    colour = "black"
  )
  
  # Include the grid lines.
  if (isTRUE(grid)) {
    panel_grid_major <- element_line(colour = gray, size = 0.5)
  } else {
    panel_grid_major <- element_blank()
  }
  
  # Remove panel border and axis ticks.
  if (isTRUE(minimal)) {
    axis_ticks <- element_blank()
    panel_border <- element_blank()
  } else {
    axis_ticks <- element_line(colour = "black")
    panel_border <- element_rect(colour = "black", fill = NA)
  }
  
  theme_linedraw(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      text = text,
      aspect.ratio = aspect_ratio,
      axis.line = element_blank(),
      axis.text = text,
      axis.text.x = element_text(angle = 90L, hjust = 1L, vjust = 0.5),
      axis.ticks = axis_ticks,
      panel.background = element_blank(),
      panel.border = panel_border,
      panel.grid.major = panel_grid_major,
      panel.grid.minor = element_blank(),
      legend.background = element_blank(),
      legend.position = legend_position,
      strip.background = element_rect(colour = NA, fill = "white"),
      strip.text = text,
      complete = TRUE,
      validate = TRUE
    )
}

theme_set(
  theme_paperwhite(
    base_size = 14L,
    legend_position = "right"
  )
)

# genes <- rtracklayer::import("../../camp_firre/analysis/Mus_musculus.GRCm38.81.gtf.gz")
# genes <- genes[which(genes$type == "gene")]
# rtracklayer::export(genes, "Mus_musculus.GRCm38.81_genes.gtf")

