% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microbiome_helpers.R
\name{add_taxonomy_column}
\alias{add_taxonomy_column}
\title{Add taxonomy label}
\usage{
add_taxonomy_column(physeq, num_species = 2)
}
\arguments{
\item{physeq}{a valid phyloseq object that contains a taxonomy table}

\item{num_species}{the number of species to retain if more than one are identified}
}
\value{
a phyloseq object with an additional column on the taxonomy
        table called "Taxonomy"
}
\description{
add a column to the taxonomy table of a phyloseq object that lists the
lowest rank taxonomy assigned to that OTU along with a prefix indicating
the taxonomic rank.
}
\details{
Example: g:Pseudomonas
}
