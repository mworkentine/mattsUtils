# ----------------------------------
#  tidy phyloseq
# ----------------------------------

#' Tidy phyloseq
#'
#' A tidy and fast version of phyloseq::psmelt, albeit with less checking
#'
#' @param physeq, a valid phyloseq object, sample names must not include
#'   'Sample' or 'OTU' or any of the taxonomic ranks (kingdom,  phylum, etc)
#'
#' @return a tidy tibble
#'
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq taxa_are_rows
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq tax_table
#' @importFrom tidyr gather
#' @importFrom dplyr as_tibble
#' @importFrom dplyr %>%
#' @importFrom dplyr full_join
#'
#' @export
tidy_phyloseq = function(physeq) {

  # --- otu table
  otutab = otu_table(physeq)@.Data
  if (!taxa_are_rows(physeq)) {
    otutab = t(otutab)
  }
  otutab = as_tibble(otutab, rownames = "OTU") %>%
    gather(Sample, Abundance, -OTU)

  # --- sample data
  samp = as_tibble(sample_data(physeq), rownames = "Sample")

  # --- taxonomy
  tax = as_tibble(tax_table(physeq)@.Data, rownames = "OTU")

  otutab %>%
    full_join(samp) %>%
    full_join(tax)
}
