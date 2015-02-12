# Utility functions written specifically for rna-seq analysis

#' Create a gene list from a DESeq results object
#'
#' This is a function that takes a DESeqResults object and a data frame gene list
#' with only the genes that have and adjusted p-valu above cutoff.  Can include a link
#' to NCBI (inspired by the Reporting tools package)
#' @param results DESeq results object
#' @param db An optional AnnotationDB object for gene identifiers
#' @param cutoff Cutoff for the adjusted p-value (or p-value if there is no adjusted). Defaults to 0.05
#' @param ncbi_link Logical, include NCBI link in the results. Defaults to FALSE.
#' @export

create_gene_list = function(results, db=NULL, cutoff = 0.05, ncbi_link=FALSE){

  require(hwriter)

  stopifnot(inherits( results, "DESeqResults" ))
  if(!is.null(db)){
    stopifnot(inherits( db, "AnnotationDb" ) )
  }

  results = results[!(is.na(results$padj)),]

  if ("padj" %in% colnames(results)){
    results = results[!(is.na(results$padj)),]
    results = results[order(results$padj),]
    results = results[results$padj < cutoff, ]
  } else if ("pvalue" %in% colnames(results)){
      results = results[order(results$pvalue),]
      results = results[results$pvalue < cutoff, ]
      results = results[!(is.na(results$pvalue)),]
  } else {
    stop("Error: Couldn't find a padj or pvalue column.")
  }

  # Add in the gene info
  results$GeneID = rownames(results)
  results = as.data.frame(results)

  if(ncbi_link){
    results$Symbol = convertIDs(results$GeneID, "ENTREZID", "SYMBOL", db)
    results$GeneName = convertIDs(results$GeneID, "ENTREZID", "GENENAME", db)
    results$Link = hwrite(as.character(results$GeneID),
                          link=paste("http://www.ncbi.nlm.nih.gov/gene/",
                                     as.character(results$GeneID), sep=''),
                          table=FALSE, target = "_blank")
  }
  return(results)

}

#' Convert between gene IDs
#'
#' from \url{http://www.bioconductor.org/help/course-materials/2014/SeattleOct2014/B02.1.1_RNASeqLab.html}
#'@export
#'
convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {


  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


#' Plot genes
#'
#' A function to plot normalized counts of genes from a DESeq object
#'@export
plotGenes = function(geneList, dds, intgroup, db, title = NULL, plotType = c("box", "point", "justData")) {


  stopifnot( inherits( db, "AnnotationDb" ) )
  stopifnot(inherits( dds, "DESeqDataSet" ))

  plotType = match.arg(plotType)

  require(ggplot2)

  # generate the count data for plotting
  data = data.frame()
  for (gene in geneList) {
    temp = plotCounts(dds, gene, intgroup = intgroup, returnData = TRUE)
    temp$Sample = rownames(temp)
    temp = cbind(temp, Gene = gene)
    data = rbind(data, temp)
  }

  colour = intgroup[1]
  if (length(intgroup) > 1) {
    shape = intgroup[2]
  } else {
    shape = NULL
  }

  data$symbol = convertIDs(as.character(data$Gene), "ENTREZID", "SYMBOL", db, "useFirst")
  #data$Sample = factor(rownames(data))

  colours = brewer.pal(n = 3, name = "YlOrRd")
  if (plotType == "box") {
    plot = ggplot(data, aes_string(x = intgroup[1], y = "count")) +
      geom_point(aes_string(colour = colour, shape = shape)) +
      geom_boxplot(aes_string(fill = colour), alpha = 0.5) +
      facet_wrap(~symbol, scales = "free") +
      theme_bw() +
      scale_colour_brewer(palette = "Set2") +
      scale_fill_brewer(palette = "Set2") +
      labs(title = title, y = "Normalized Counts", x = "")
  } else if (plotType == "point"){
    plot = ggplot(data, aes_string(x = intgroup[1], y = "count")) +
      geom_point(aes_string(colour = colour, shape = shape)) +
      facet_wrap(~symbol, scales = "free") +
      theme_bw() +
      scale_colour_brewer(palette = "Set2") +
      scale_fill_brewer(palette = "Set2") +
      labs(title = title, y = "Normalized Counts", x = "")
  } else {
    return(data)
  }

  return(plot)
}


#' Plot most variable genes
#'
#' This function takes a DESeq results object that has been transformed (rlog)
#' and plots the top variable genes.  Also mean centers the data
#' @export
plot_variable_genes = function(rld_obj, annotation_vars = NULL, top = 100, onlyGenes = FALSE, ...) {


  require(genefilter)
  require(RColorBrewer)

  stopifnot(inherits(rld_obj, "SummarizedExperiment"))

  topVarGenes = head(order(-rowVars(assay(rld_obj))),top)
  colors = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)


  mat = assay(rld_obj)[ topVarGenes, ]
  mat = mat - rowMeans(mat)
  #rownames(mat) = convertIDs(row.names(mat), from = "ENTREZID", "SYMBOL", org.Hs.eg.db)

  if(onlyGenes){
    return(row.names(mat))
  }

  if (!is.null(annotation_vars)){

    colour_vars = as.character(as.data.frame(unique((rld_obj@colData)[annotation_vars[1]]))[,1])
    colours = brewer.pal(9, "Set1")[1:length(colour_vars)]
    names(colours) = colour_vars
    ann_colours = list(colours)
    names(ann_colours) = annotation_vars[1]

    pheatmap(mat, scale = "row", border_color = NA,
             annotation = as.data.frame(rld_obj@colData)[annotation_vars],
             annotation_colors = ann_colours, ...)
  } else {
     pheatmap(mat, scale = "row", border_color = NA, ...)
  }
}


#' Plot differentially expressed genes
#' This function takes a DESeqResults object and contrasts argument and plots
#' the differentially expressed genes at padj cutoff
#'@export
plot_de_genes = function(dds, contrast, cutoff = 0.05, ...){

  require(pheatmap)
  require(RColorBrewer)

  if (missing(contrast)){
    stop("The contrast argument is required so the plot can be annotated correctly")
  }

  stopifnot(inherits(dds, "DESeqDataSet"))
  results = results(dds, contrast)
  results = results[!(is.na(results$padj)),]
  results = results[results$padj < cutoff, ]

  colour_vars = contrast[2:length(contrast)]
  colours = brewer.pal(9, "Set1")[1:length(colour_vars)]
  names(colours) = colour_vars
  ann_colours = list(colours)
  names(ann_colours) = contrast[1]

  rld = rlog(dds)
  mat = assay(rld)[rownames(results),]
  mat = mat - rowMeans(mat)
  #rownames(mat) = convertIDs(row.names(mat), from = "ENTREZID", "SYMBOL", org.Hs.eg.db)

  pheatmap(mat, scale = "row", border_color = NA,
       annotation = as.data.frame(rld@colData)[contrast[1]],
       annotation_colors = ann_colours, ...)


}

#' Plot PCA (modified)
#'
#' A slight modification from the DESeq function \code{\link[DESeq2]{plotPCA}} to return the 3 PC also
#' @export
plotPCA2 = function (x, intgroup = "condition", ntop = 500, returnData = FALSE)
{
    require(genefilter)

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

