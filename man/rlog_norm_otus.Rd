% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microbiome_helpers.R
\name{rlog_norm_otus}
\alias{rlog_norm_otus}
\title{rlog normalize an OTU table}
\usage{
rlog_norm_otus(physeq, remove_negs = FALSE)
}
\arguments{
\item{physeq}{valid phyloseq object}

\item{remove_negs}{remove negative values from the normalized OTU table (TRUE). Defaults to FALSE}
}
\value{
The original phyloseq object with the normalized OTU table
}
\description{
Use the \code{\link{rlog}} function from DESeq2 to normalize an OTU table
}
