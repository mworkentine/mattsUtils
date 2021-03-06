% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_helpers.R
\name{make_heatmap}
\alias{make_heatmap}
\title{Make Heatmap}
\usage{
make_heatmap(genes, ddr, gene_key = NULL, ...)
}
\arguments{
\item{genes}{character vector of gene identifers that you wish to plot.  ID
type must match gene ids in the dds object.}

\item{ddr}{valid DESeqTransform object used to extract the gene expression
values.  Recommended transform for visualization is the rlog transform.
See the DESeq2 vignette for more information.}

\item{gene_key}{a named list of gene symbols with the gene ids as names.
Used as a lookup key to replace gene ids with symbols in the plot}

\item{...}{additional parameters passed to \code{\link{pheatmap}}}
}
\description{
Create a heatmap with the provide list of genes
}
