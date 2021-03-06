% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_monophyletic_subsets.R
\name{get_monophyletic_subsets}
\alias{get_monophyletic_subsets}
\title{Get monophyletic subsets}
\usage{
get_monophyletic_subsets(tree, tips)
}
\arguments{
\item{tree}{phylo, a phylogentic tree}

\item{tips}{character vector of tips in the tree}
}
\value{
a list
}
\description{
Given a phylogenetic tree and a subset of tips in that tree, find all
monophyletic groups within that subset
}
