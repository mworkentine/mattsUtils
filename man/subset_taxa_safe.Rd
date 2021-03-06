% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microbiome_helpers.R
\name{subset_taxa_safe}
\alias{subset_taxa_safe}
\title{Subset taxa 2}
\usage{
subset_taxa_safe(physeq, rank, taxa)
}
\arguments{
\item{physeq}{a valid phyloseq object}

\item{rank}{the taxonomic rank of the taxa to subset}

\item{taxa}{the taxa name to subset to}
}
\description{
A programming safe version of \code{\link{phyloseq::subset_taxa}}
}
