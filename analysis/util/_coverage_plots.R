# Credit: signag (https://satijalab.org/signac/articles/visualization.html)
library(fastmatch)
library(GenomicRanges)
library(RcppRoll)

AnnotationPlot <- function(annotation, region, mode = "gene") {
  if(mode == "gene") {
    collapse_transcript <- TRUE
    label <- "gene_name"
  } else if (mode == "transcript") {
    collapse_transcript <- FALSE
    label <- "tx_id"
  } else {
    stop("Unknown mode requested, choose either 'gene' or 'transcript'")
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  if (mode == "gene") {
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
      fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
    ]
  } else {
    tx.keep <- unique(x = annotation.subset$tx_id)
    annotation.subset <- annotation[
      fmatch(x = annotation$tx_id, table = tx.keep, nomatch = 0L) > 0L
    ]
  }
  
  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      collapse_transcript = collapse_transcript
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 3
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.2
    p <- p + geom_text(
      data = annotation_df_list$labels,
      mapping = aes_string(x = "position", y = "dodge", label = label),
      size = 2.5
    )
    y_limit <- c(0.9, n_stack + 0.4)
  }
  p <- p +
    theme_classic() +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(start.pos, end.pos) +
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}


reformat_annotations <- function(
  annotation,
  start.pos,
  end.pos,
  collapse_transcript = TRUE
) {
  total.width <- end.pos - start.pos
  tick.freq <- total.width / 50
  annotation <- annotation[annotation$type == "exon"]
  exons <- base::as.data.frame(x = annotation)
  if (collapse_transcript) {
    annotation <- split(
      x = annotation,
      f = annotation$gene_name
    )
  } else {
    annotation <- split(
      x = annotation,
      f = annotation$tx_id
    )
  }
  annotation <- lapply(X = annotation, FUN = base::as.data.frame)
  
  # add gene total start / end
  gene_bodies <- list()
  for (i in seq_along(annotation)) {
    df <- data.frame(
      seqnames = annotation[[i]]$seqnames[[1]],
      start = min(annotation[[i]]$start),
      end = max(annotation[[i]]$end),
      strand = annotation[[i]]$strand[[1]],
      tx_id = annotation[[i]]$tx_id[[1]],
      gene_name = annotation[[i]]$gene_name[[1]],
      gene_biotype = annotation[[i]]$gene_biotype[[1]],
      type = "body"
    )
    # trim any that extend beyond region
    df$start <- ifelse(
      test = df$start < start.pos,
      yes = start.pos,
      no = df$start
    )
    df$end <- ifelse(
      test = df$end > end.pos,
      yes = end.pos,
      no = df$end
    )
    breaks <- split_body(df = df, width = tick.freq)
    df <- rbind(df, breaks)
    gene_bodies[[i]] <- df
  }
  gene_bodies <- do.call(what = rbind, args = gene_bodies)
  
  # record if genes overlap
  overlap_idx <- record_overlapping(
    annotation = gene_bodies,
    min.gapwidth = 1000,
    collapse_transcript = collapse_transcript
  )
  # overlap_idx <- overlap_idx
  if (collapse_transcript) {
    gene_bodies$dodge <- overlap_idx[gene_bodies$gene_name]
    exons$dodge <- overlap_idx[exons$gene_name]
  } else {
    gene_bodies$dodge <- overlap_idx[gene_bodies$tx_id]
    exons$dodge <- overlap_idx[exons$tx_id]
  }
  
  label_df <- gene_bodies[gene_bodies$type == "body", ]
  label_df$width <- label_df$end - label_df$start
  label_df$position <- label_df$start + (label_df$width / 2)
  
  onplus <- gene_bodies[gene_bodies$strand %in% c("*", "+"), ]
  onminus <- gene_bodies[gene_bodies$strand == "-", ]
  
  return(
    list(
      "labels" = label_df,
      "exons" = exons,
      "plus" = onplus,
      "minus" = onminus
    )
  )
}

#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
record_overlapping <- function(
  annotation,
  min.gapwidth = 1000,
  collapse_transcript = TRUE
) {
  # convert back to granges
  annotation$strand <- "*"
  gr <- makeGRangesFromDataFrame(
    df = annotation[annotation$type == "body", ], keep.extra.columns = TRUE
  )
  # work out which ranges overlap
  collapsed <- GenomicRanges::reduce(
    x = gr, with.revmap = TRUE, min.gapwidth = min.gapwidth
  )$revmap
  idx <- seq_along(gr)
  for (i in seq_along(collapsed)) {
    mrg <- collapsed[[i]]
    for (j in seq_along(mrg)) {
      idx[[mrg[[j]]]] <- j
    }
  }
  if (collapse_transcript) {
    names(x = idx) <- gr$gene_name
  } else {
    names(x = idx) <- gr$tx_id
  }
  return(idx)
}

GetGRangesFromEnsDb <- function(
  ensdb,
  standard.chromosomes = TRUE,
  biotypes = c("protein_coding", "lincRNA", "rRNA", "processed_transcript"),
  verbose = TRUE
) {
  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("Please install biovizBase\n",
         "https://www.bioconductor.org/packages/biovizBase/")
  }
  # convert seqinfo to granges
  whole.genome <-  as(object = seqinfo(x = ensdb), Class = "GRanges")
  if (standard.chromosomes) {
    whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
  }
  
  # extract genes from each chromosome
  my_lapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  tx <- my_lapply(X = seq_along(whole.genome), FUN = function(x){
    suppressMessages(expr = biovizBase::crunch(
      obj = ensdb,
      which = whole.genome[x],
      columns = c("tx_id", "gene_name", "gene_id", "gene_biotype")))
  })
  
  # combine
  tx <- do.call(what = c, args = tx)
  tx <- tx[tx$gene_biotype %in% biotypes]
  return(tx)
}


split_body <- function(df, width = 1000) {
  wd <- df$end - df$start
  nbreak <- wd / width
  if (nbreak > 1) {
    steps <- 0:(nbreak)
    starts <- (width * steps) + df$start
    starts[starts > df$end] <- NULL
  } else {
    starts <- df$end
  }
  breaks <- data.frame(
    seqnames = df$seqnames[[1]],
    start = starts,
    end = starts + 1,
    strand = df$strand[[1]],
    tx_id = df$tx_id[[1]],
    gene_name = df$gene_name[[1]],
    gene_biotype = df$gene_biotype[[1]],
    type = "arrow"
  )
  return(breaks)
}

CombineTracks <- function(
  plotlist,
  expression.plot = NULL,
  heights = NULL,
  widths = NULL
) {
  # remove any that are NULL
  nullplots <- sapply(X = plotlist, FUN = is.null)
  plotlist <- plotlist[!nullplots]
  heights <- heights[!nullplots]
  
  if (length(x = plotlist) == 1) {
    return(plotlist[[1]])
  }
  
  # remove x-axis from all but last plot
  for (i in 1:(length(x = plotlist) - 1)) {
    plotlist[[i]] <- plotlist[[i]] + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    )
  }
  
  # combine plots
  if (is.null(x = heights)) {
    # set height of first element to 10x more than other elements
    n.plots <- length(x = plotlist)
    heights <- c(8, rep(1, n.plots - 1))
  } else {
    if (length(x = heights) != length(x = plotlist)) {
      stop("Relative height must be supplied for each plot")
    }
  }
  if (!is.null(x = expression.plot)) {
    # align expression plot with the first element in plot list
    p <- (plotlist[[1]] + expression.plot) +
      plot_layout(widths = widths)
    
    n <- length(x = plotlist)
    heights.2 <- heights[2:n]
    p2 <- wrap_plots(plotlist[2:n], ncol = 1, heights = heights.2)
    
    p <- p + p2 + guide_area() + plot_layout(
      ncol = 2, heights = c(heights[[1]], sum(heights.2)),
      guides = "collect")
  } else {
    p <- wrap_plots(plotlist, ncol = 1, heights = heights)
  }
  return(p)
}

BigwigTrack <- function(
  region,
  bigwig,
  smooth = 200,
  extend.upstream = 0,
  extend.downstream = 0,
  type = "coverage",
  y_label = "bigWig",
  bigwig.scale = "common",
  ymax = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1
) {
  if (!inherits(x = bigwig, what = "list")) {
    bigwig <- list("bigWig" = bigwig)
  }
  possible_types <- c("line", "heatmap", "coverage")
  if (!(type %in% possible_types)) {
    stop(
      "Invalid type requested. Choose ",
      paste(possible_types, collapse = ", ")
    )
  }
  if (.Platform$OS.type == "windows") {
    message("BigwigTrack not supported on Windows")
    return(NULL)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
    return(NULL)
  }
  region <- FindRegion(
    object = NULL,
    region = region,
    sep = c("-", "-"),
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  if (!inherits(x = region, what = "GRanges")) {
    stop("region should be a GRanges object")
  }
  all.data <- data.frame()
  for (i in seq_along(bigwig)) {
    region_data <- rtracklayer::import(
      con = bigwig[[i]],
      which = region,
      as = "NumericList"
    )[[1]]
    if (!is.null(x = smooth)) {
      region_data <- roll_mean(x = region_data, n = smooth, fill = 0L)
    }
    region_data <- data.frame(
      position = start(x = region):end(x = region),
      score = region_data,
      stringsAsFactors = FALSE,
      bw = names(x = bigwig)[[i]]
    )
    if (bigwig.scale == "separate") {
      # scale to fraction of max for each separately
      file.max <- max(region_data$score, na.rm = TRUE)
      region_data$score <- region_data$score / file.max
    }
    all.data <- rbind(all.data, region_data)
  }
  all.data$bw <- factor(x = all.data$bw, levels = names(x = bigwig))
  window.size = width(x = region)
  sampling <- max(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = all.data, n = sampling)
  
  covmax <- signif(x = max(coverages$score, na.rm = TRUE), digits = 2)
  if (is.null(x = ymax)) {
    ymax <- covmax
  } else if (is.character(x = ymax)) {
    if (!startsWith(x = ymax, prefix = "q")) {
      stop("Unknown ymax requested. Must be NULL, a numeric value, or 
           a quantile denoted by 'qXX' with XX the desired quantile value,
           e.g. q95 for 95th percentile")
    }
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(x = ymax))
    ) / 100
    ymax <- covmax * percentile.use
  }
  ymin <- 0
  
  # perform clipping
  coverages$score[coverages$score > ymax] <- ymax 
  
  if (type == "line") {
    p <- ggplot(
      data = coverages,
      mapping = aes_string(x = "position", y = "score", color = "bw")
    ) + geom_line() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1) +
      scale_color_grey()
  } else if (type == "heatmap") {
    # different downsampling needed for heatmap
    # cut into n bins and average within each bin
    all.data$bin <- floor(x = all.data$position / smooth)
    all.data <- group_by(all.data, bin, bw)
    all.data <- mutate(all.data, score = mean(x = score))
    all.data <- ungroup(all.data)
    all.data <- unique(x = all.data[, c("bin", "score", "bw")])
    p <- ggplot(
      data = all.data,
      mapping = aes_string(x = "bin", y = 1, fill = "score")
    ) + geom_tile() + scale_fill_viridis_c() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1)
  } else if (type == "coverage") {
    p <- ggplot(
      data = coverages,
      mapping = aes_string(x = "position", y = "score", fill = "bw")
    ) + geom_area() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1) +
      scale_fill_grey()
  }
  chromosome <- as.character(x = seqnames(x = region))
  p <- p + theme_browser(axis.text.y = TRUE) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = y_label)
  return(p)
}

BigwigTrack_stranded_allreps <- function(
  region,
  bigwig,
  smooth = 200,
  extend.upstream = 0,
  extend.downstream = 0,
  type = "coverage",
  y_label = "bigWig",
  bigwig.scale = "common",
  ymax = NULL,
  ymin = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1
) {
  if (!inherits(x = bigwig, what = "list")) {
    bigwig <- list("bigWig" = bigwig)
  }
  possible_types <- c("line", "heatmap", "coverage")
  if (!(type %in% possible_types)) {
    stop(
      "Invalid type requested. Choose ",
      paste(possible_types, collapse = ", ")
    )
  }
  
  region <- FindRegion(
    object = NULL,
    region = region,
    sep = c("-", "-"),
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  if (!inherits(x = region, what = "GRanges")) {
    stop("region should be a GRanges object")
  }
  all.data <- data.frame()
  for (i in seq_along(bigwig)) {
    for(strand in c("pos", "neg")) {
      region_data <- rtracklayer::import(
        con = bigwig[[i]][[strand]],
        which = region,
        as = "NumericList"
      )[[1]]
      if (!is.null(x = smooth)) {
        region_data <- roll_mean(x = region_data, n = smooth, fill = 0L)
      }
      region_data <- data.frame(
        position = start(x = region):end(x = region),
        score = region_data,
        stringsAsFactors = FALSE,
        bw = names(x = bigwig)[[i]],
        strand = strand
      )
      if (bigwig.scale == "separate") {
        # scale to fraction of max for each separately
        # TODO: this may need to be fixed for strandedness.
        file.max <- max(abs(region_data$score), na.rm = TRUE)
        region_data$score <- region_data$score / file.max
      }
      all.data <- rbind(all.data, region_data)
    }
  }
  all.data$bw <- factor(x = all.data$bw, levels = names(x = bigwig))
  window.size = width(x = region)
  sampling <- max(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = all.data, n = sampling)
  
  covmax <- signif(x = max(coverages$score, na.rm = TRUE), digits = 2)
  covmin <- signif(x = min(coverages$score, na.rm = TRUE), digits = 2)
  if (is.null(x = ymax)) {
    ymax <- covmax
  } else if (is.character(x = ymax)) {
    if (!startsWith(x = ymax, prefix = "q")) {
      stop("Unknown ymax requested. Must be NULL, a numeric value, or 
           a quantile denoted by 'qXX' with XX the desired quantile value,
           e.g. q95 for 95th percentile")
    }
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(x = ymax))
    ) / 100
    ymax <- covmax * percentile.use
  }
  
  if (is.null(x = ymin)) {
    ymin <- covmin
  } else if (is.character(x = ymin)) {
    if (!startsWith(x = ymin, prefix = "q")) {
      stop("Unknown ymin requested. Must be NULL, a numeric value, or 
           a quantile denoted by 'qXX' with XX the desired quantile value,
           e.g. q95 for 95th percentile")
    }
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(x = ymin))
    ) / 100
    ymax <- covmin * percentile.use
  }
  
  # perform clipping
  coverages$score[coverages$score > ymax] <- ymax 
  coverages$score[coverages$score < ymin] <- ymin
  
  coverages <- coverages %>%
    separate(bw, into = c("timepoint", "replicate"), remove = FALSE)
  
  if (type == "line") {
    p <- ggplot(
      data = coverages,
      mapping = aes_string(x = "position", y = "score", color = "bw")
    ) + geom_line() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1) +
      scale_color_grey()
  } else if (type == "heatmap") {
    # different downsampling needed for heatmap
    # cut into n bins and average within each bin
    all.data$bin <- floor(x = all.data$position / smooth)
    all.data <- group_by(all.data, bin, bw)
    all.data <- mutate(all.data, score = mean(x = score))
    all.data <- ungroup(all.data)
    all.data <- unique(x = all.data[, c("bin", "score", "bw")])
    p <- ggplot(
      data = all.data,
      mapping = aes_string(x = "bin", y = 1, fill = "score")
    ) + geom_tile() + scale_fill_viridis_c() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1)
  } else if (type == "coverage") {
    p <- ggplot(
      data = coverages,
      mapping = aes_string(x = "position", y = "score", fill = "bw")
    ) + geom_area(data = coverages %>% filter(strand == "pos"), fill = discrete_pal1_sns[[1]]) +
      geom_area(data = coverages %>% filter(strand == "neg"), fill = discrete_pal1_sns[[4]]) +
      facet_nested_wrap(timepoint + replicate ~ ., strip = strip_nested(bleed = TRUE), ncol = 1,
                        strip.position = "left", nest_line = element_line(), shrink = FALSE) +
      scale_fill_grey()
  }
  chromosome <- as.character(x = seqnames(x = region))
  p <- p + theme_browser(axis.text.y = FALSE, legend = TRUE) +
    ylab(label = y_label) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    theme(panel.spacing.y = unit(x = 0.1, units = "line"))
  return(p)
}

FindRegion <- function(
  object,
  region,
  sep = c("-", "-"),
  assay = NULL,
  extend.upstream = 0,
  extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # first try to convert to coordinates, if not lookup gene
    region <- tryCatch(
      expr = suppressWarnings(
        expr = StringToGRanges(regions = region, sep = sep)
      ),
      error = function(x) {
        region <- LookupGeneCoords(
          object = object,
          assay = assay,
          gene = region
        )
        return(region)
      }
    )
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}


#' Extend
#'
#' Resize GenomicRanges upstream and or downstream.
#' From \url{https://support.bioconductor.org/p/78652/}
#'
#' @param x A range
#' @param upstream Length to extend upstream
#' @param downstream Length to extend downstream
#' @param from.midpoint Count bases from region midpoint,
#' rather than the 5' or 3' end for upstream and downstream
#' respectively.
#'
#' @importFrom GenomicRanges trim
#' @importFrom BiocGenerics start strand end width
#' @importMethodsFrom GenomicRanges strand start end width
#' @importFrom IRanges ranges IRanges "ranges<-"
#' @export
#' @concept utilities
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#' @examples
#' Extend(x = blacklist_hg19, upstream = 100, downstream = 100)
Extend <- function(
  x,
  upstream = 0,
  downstream = 0,
  from.midpoint = FALSE
) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- midpoints + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  } else {
    new_start <- start(x = x) - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- end(x = x) + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  }
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}