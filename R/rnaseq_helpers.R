
#' Single Gene Plot
#'
#' Make a plot of a single gene.  If more genes are provided will facet by gene.  Generally useful
#' for less than about 30 genes
#'
#' @param genes character vector of gene identifers that you wish to plot.  ID type must match gene
#'   ids in the dds object.
#' @param dds valid DESeq2 object used to extract the gene expression values.  Uses the
#'   \code{\link{plotCounts}} function from DESeq2 to extract the expression values.
#' @param xaxis character, name of the sample column to put on the x-axis, must be a column in
#'   \code{colData(dds)}.
#' @param fill character, name of the sample column to use for colouring the plot, must be a column
#'   in \code{colData(dds)}.
#' @param normalized whether the counts should be normalized by size factor (default is TRUE).
#'   Passed to \code{\link{plotCounts}}
#' @param transform whether to present log2 counts (TRUE) or to present the counts on the log scale
#'   (FALSE, default).  Passed to \code{\link{plotCounts}}.
#' @param gene_key an optional named list of gene symbols with the gene ids as names.  Used as a
#'   lookup key to replace gene ids with symbols in the plot
#' @param palette character, name of RColorBrewer palette to use. Default is "Set1"
#' @param type character, type of plot to make.  Defaults to boxplot.
#'
#' @export
#'
single_gene_plot = function(genes, dds, xaxis, fill, normalized = TRUE, transform = FALSE, gene_key = NULL,
                            palette = "Set1", type = c("box", "dot")) {

  gene_data = lapply(genes, function(x) {
    plotCounts(dds, x, c(xaxis, fill), returnData = TRUE, normalized = normalized, transform = transform)
  })


  if (!is.null(gene_key)) {
    names(gene_data) = gene_key[genes]
  } else {
    names(gene_data) = genes
    }

  gene_data = bind_rows(gene_data, .id = 'Gene')

  if (transform) {
    gene_data = gene_data %>%
      mutate(count = log2(count))
  }

  if (type == "box") {
	  p = ggplot(gene_data, aes_string(x = xaxis, y = "count", fill = fill)) +
  	 geom_boxplot() +
	   scale_fill_brewer(palette = palette) +
  	 facet_wrap(~Gene, scales = "free_y")


  } else if (type == "dot") {
    p = ggplot(gene_data, aes_string(x = xaxis, y = "count", colour = fill)) +
      geom_point() +
      geom_jitter(width = 0.25, height = 0) +
      scale_colour_brewer(palette = palette) +
      facet_wrap(~Gene, scales = "free_y")
  } else {
  	stop("Unrecognized plot type")
  }

  if (transform) {
    p = p + labs(y = expression(Log[2]~Normalized~Counts))
  } else {
    p = p + labs(y = "Normalized Counts") + scale_y_continuous(trans = "log2")
  }
 	return(p)

}


#' Make Heatmap
#'
#' Create a heatmap with the provide list of genes
#'
#' @param genes character vector of gene identifers that you wish to plot.  ID type must match gene
#'   ids in the dds object.
#' @param ddr valid DESeqTransform object used to extract the gene expression values.  Recommended transform for visualization is the rlog transform.  See the DESeq2 vignette for more information.
#' @param gene_key a named list of gene symbols with the gene ids as names.  Used as a lookup key to replace gene ids with symbols in the plot
#' @param ... additional parameters passed to \code{\link{pheatmap}}
#'
#' @export
#'
make_heatmap = function(genes, ddr, gene_key = NULL, ...) {

  genes = genes[genes %in% rownames(assay(ddr))]
  mat = assay(ddr)[genes, ]
  if (!is.null(gene_key)) {
  	rownames(mat) = unname(gene_key[rownames(mat)])
  }

  pheatmap(mat, ...)

}


#' Diagnostic plots for a DESeq differential expression analysis
#'
#' Given a valid DESeq results object this function will make a p-value
#' histogram and a volcano plot.  The p-value histogram is constructed with the
#' p-values prior to multiple-testing correction and is used to diagnose if
#' something went wrong with your testing procedure.  See
#' \url{http://varianceexplained.org/statistics/interpreting-pvalue-histogram/}
#' for more details on interpreting the p-value histogram.  The volcano plot
#' gives you an idea of how many genes are differentially expressed in both
#' directions compared to their expression levels.
#' @param res A DESeq results object
#' @param cutoff The p value to use as a cutoff for significantly differentially
#'   expressed genes.  Defaults to 0.05
#' @param show_title If specified the description of the results object will be
#'   printed as the plot's title
#'
#' @export
#'
plot_diagnositics = function(res, cutoff = 0.05, show_title = TRUE){

  if (missing(res)) {
		stop("Please provide DESeq results object")
	}
	if (!inherits(res, "DESeqResults")) {
		stop("Not a valid DESeq results object")
	}

  # create a new column to define genes that pass the cutoff for colouring
  res$sigcolour = ifelse(res$padj < cutoff, paste0("< ", cutoff), paste0("> ", cutoff))

  # p-value histogram plot
  p_hist = ggplot(as.data.frame(res), aes(x = pvalue)) + geom_histogram()

  # need to define x-axis limits manually so they are symetrical around zero
  max_lfc = as.data.frame(res) %>% summarise(max(abs(log2FoldChange), na.rm = TRUE))
  xlims = c(-max_lfc - (max_lfc*0.1), max_lfc + (max_lfc*0.1))
  xlims = round(as.numeric(xlims))

  # volcano plot
  volcano = ggplot(as.data.frame(res), aes(y = -log10(padj), x = log2FoldChange, colour = sigcolour)) +
    geom_point() + scale_colour_discrete("Cutoff") + xlim(xlims)

  # put them together
  if (show_title) {
    title = res@elementMetadata@listData$description[2]
    gridExtra::grid.arrange(p_hist, volcano, ncol = 2, top = title)
  } else {
    gridExtra::grid.arrange(p_hist, volcano, ncol = 2)
  }
}

#' Plot stats from featureCounts output
#'
#' Takes the count summary data from a \code{\link[Rsubread]{featureCounts}}
#' object. Alternitavley you can load the counts.txt.summary file as a dataframe
#' and provide this if using the command line version of featureCounts.
#'
#' @param stats The count summary data
#' @param palette A valid RColorBrewer palette
#'
#' @export
#'
plot_featureCount_stats = function(stats, palette = "Set3"){

  stats %>% gather(Sample, count, -starts_with("Status")) %>%
    group_by(Sample) %>%
    mutate(total_reads = sum(count)) %>%
    group_by(Status, Sample) %>%
    mutate(percentage = count/total_reads) %>%
    ggplot(aes(x = Sample, y = percentage, fill = Status)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_brewer(palette = palette) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(y = "Percentage of total reads") +
      scale_y_continuous(labels = percent)

}

#' Plot PCA (modified)
#'
#' A slight modification from the DESeq function \code{\link[DESeq2]{plotPCA}} to return the 3 PC also
#'
#' @export
#'
plotPCA2 = function (x, intgroup = "condition", ntop = 500, returnData = FALSE)
{

    rv <- rowVars(assay(x))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
        length(rv)))]
    pca <- prcomp(t(assay(x)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(x)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
    group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group,
        intgroup.df, names = colnames(x))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:3]
        return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *
        100), "% variance"))
}
