#' Make a heatmap of a list of genes from a DESeq analysis
#'
#' @param rld A summarized experiment object from DESeq that has been transformed with rlog
#' @param genes A character vector of genes to plot.  Must match gene names in the rld object.
#' @param anno_var The metadat variable to annotate the columns with.  Must be a column in rld
#'   column data
#' @param subset_samples A character vector of samples to subset the data.  Names must match the rld
#'   object
#' @param title A plot title
#' @param ... Other arguments passed to Heatmap
#'
#' @import ComplexHeatmap
#' @export

make_heatmap = function(rld, genes, anno_var, subset_samples = NULL, title = "", ...){

  if(missing(anno_var)){
    stop("Please provide annotation variable")
  }

  if(!is.null(subset_samples)){
    sample_data = as.data.frame(rld@colData[subset_samples,])[anno_var]
    mat = assay(rld)[genes, subset_samples]
  } else {
    sample_data = as.data.frame(rld@colData)[anno_var]
    mat = assay(rld)[genes, ]
  }

  mat = mat - rowMeans(mat)
  #rownames(mat) = convertIDs(row.names(mat), from = "ENTREZID", "SYMBOL", db)
  pal = brewer.pal(9, "Set1")

  colours = unique(sample_data[[1]])
  if(is.factor(colours)){
    colours = as.character(colours)
  }

  pal = c(pal[1:length(colours)])
  names(pal) = colours
  sample_colours = list(pal)
  names(sample_colours) = anno_var
  ann = HeatmapAnnotation(df = sample_data, col = sample_colours)

  breaks = seq(min(mat),max(mat), 0.01)
  make_colours = colorRampPalette(wes_palette("Zissou", n = 10, type = "continuous"))

  Heatmap(mat,
          name = "Normalized Counts",
          row_names_gp = gpar(fontsize = 8),
          top_annotation = ann,
          col = circlize::colorRamp2(breaks, make_colours(length(breaks))),
          column_title = title,
          ...)
}
