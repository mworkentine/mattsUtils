#'   Diagnostic plots for a DESeq differential expression analysis
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
#' @export

plot_diagnositics = function(res, cutoff = 0.05, show_title = TRUE){

  if(missing(res)){
		stop("Please provide DESeq results object")
	}
	if(!inherits(res, "DESeqResults")){
		stop("Not a valid DESeq results object")
	}
  # create a new column to define genes that pass the cutoff for colouring
  res$sigcolour = ifelse(res$padj < cutoff, paste0("< ", cutoff), paste0("> ", cutoff))
  p_hist = ggplot(as.data.frame(res), aes(x = pvalue)) + geom_histogram()
  # need to define x-axis limits manually so they are symetrical around zero
  max_lfc = as.data.frame(res) %>% summarise(max(abs(log2FoldChange), na.rm = TRUE))
  xlims = c(-max_lfc - (max_lfc*0.1), max_lfc + (max_lfc*0.1))
  xlims = round(as.numeric(xlims))
  volcano = ggplot(as.data.frame(res), aes(y = -log10(padj), x = log2FoldChange, colour = sigcolour)) +
    geom_point() + scale_colour_discrete("Cutoff") + xlim(xlims)
  if(show_title){
    title = res@elementMetadata@listData$description[2]
    grid.arrange(p_hist, volcano, ncol = 2, top = title)
  } else {
    grid.arrange(p_hist, volcano, ncol = 2)
  }
}
