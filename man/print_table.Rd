% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_table.R
\name{print_table}
\alias{print_table}
\title{Print table}
\usage{
print_table(table, fname = "table", pageLength = 25L, digits = 3)
}
\arguments{
\item{table}{dataframe}

\item{fname}{filename to use for excel download}

\item{pageLength}{integer, how many rows to display by default}

\item{digits}{number of signficant digits to show}
}
\value{
datatable
}
\description{
Prints a nicely formatted DT::datatable with an option to download the excel file
}
