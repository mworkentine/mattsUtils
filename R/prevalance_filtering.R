############################################################
#                                                          #
#                   Prevalance functions                   #
#                                                          #
############################################################


#' Make a prevalance dataframe
#'
#' Function that takes a phyloseq object and creates a dataframe indicating
#' the prevalance of each OTU
#'
#' @param physeq a valid phyloseq object
#' @param prev_def prevalance definition, the minimum count for an OTU to defined
#'   as present
#' @export
#'
make_prev_df = function(physeq, prev_def) {
  physeq %>% psmelt() %>%
    group_by(Phylum, Genus, OTU) %>%
    summarise(prev = sum(Abundance >= prev_def), tot = sum(Abundance))
}

#' Prevalance filter
#'
#' Filter OTUs base on prevalance, defined as the number of samples the OTU is
#' considered 'present'
#'
#' @inheritParams make_prev_df
#' @param cutoff the proportion of samples that an OTU must be present in to be
#'   included
#'
#' @return a phyloseq object with OTUs below the cutoff filtered out
#' @export
#'
prev_filt = function(physeq, prev_def, cutoff) {

  prevdf = make_prev_df(physeq, prev_def)

  # define prevelance threshold and filter data
  cutoff_num = cutoff * nsamples(physeq)
  keeps = dplyr::filter(prevdf, prev >= cutoff_num)$OTU

  new_data = prune_taxa(keeps, physeq)
  message("Filtered ", ntaxa(physeq) - ntaxa(new_data), " taxa from original ",
          ntaxa(physeq), ", leaving ", ntaxa(new_data), " taxa")
  return(new_data)
}


#' Prevelance Plot
#'
#' Plot the prevalance of each OTU, faceted by phyla
#'
#' @inheritParams prev_filt
#' @param phylum phylum/phyla to subset
#' @return ggplot object
#' @export
#'
prev_plot = function(physeq, cutoff = 0.10, prev_def = 1, phylum = NULL) {
  d = make_prev_df(physeq, prev_def)
  if (!is.null(phylum)) {
    d = dplyr::filter(d, Phylum %in% phylum)
  }
  d %>%
    ggplot(aes(x = tot, y = prev / nsamples(physeq))) +
      geom_point(aes(colour = Phylum), alpha = 0.6) +
      geom_hline(yintercept = cutoff, alpha = 0.5, linetype = 2) +
      facet_wrap(~Phylum) +
      scale_x_log10() +
      scale_color_viridis(discrete = TRUE, option = "plasma") +
      theme(strip.background = element_blank(), legend.position = "none") +
      labs(x = "Total Abundance", y = "Fraction of total samples")
}
