
#' rlog normalize an OTU table
#'
#' Use the \code{\link{rlog}} function from DESeq2 to normalize an OTU table
#'
#' @param physeq valid phyloseq object
#' @param remove_negs remove negative values from the normalized OTU table (TRUE). Defaults to FALSE
#' @return The original phyloseq object with the normalized OTU table
#' @export
#'
rlog_norm_otus = function(physeq, remove_negs = FALSE) {


  physeq_norm = physeq
  otu_table(physeq_norm) = otu_table(physeq) + 1

  dds = phyloseq_to_deseq2(physeq_norm, ~1)
  ddr = rlog(dds, blind = TRUE, fitType = "local")
  counts_norm = as.matrix(assay(ddr))

  if (remove_negs) {
    counts_norm[counts_norm < 0] = 0
  }

  otu_table(physeq_norm) = otu_table(counts_norm,taxa_are_rows = TRUE)

  return(physeq_norm)

}


#' Log fold-change plot
#'
#' Create plot that shows the log-fold change for various taxa
#'
#' @param x results from DESeq2 differential analysis, processed with biobroom
#' @param tax_level the taxonomy level to group the y-axis by
#' @param colour character, the taxonomy level to colour the points
#' @param errors logical, should the error bars be included, defaults to FALSE
#'
#' @export
#'
lfc_plot = function(x, tax_level = "Genus", colour = "Phylum", errors = FALSE) {

  new_levels = unique(x[[tax_level]][order(x$estimate)])
  x[[tax_level]] = factor(x[[tax_level]], levels = new_levels)
  #x %>% mutate(Genus = factor(Genus, levels = unique(Genus[order(estimate)])))  %>%
   p = x %>%
     ggplot(aes_string(x = "estimate", y = tax_level, colour = colour)) + geom_point(size = 3) +
      geom_vline(xintercept = 0) +
      scale_color_brewer(palette = "Set1")

   if (errors) {
     p = p + geom_errorbarh(aes(xmin = estimate - stderror, xmax = estimate + stderror), height = 0.5)
   }

   return(p)
}


#' Get a list of significant taxa from a DESeq analysis
#'
#' @param dds valid DESeq2 object
#' @param physeq valid phyloseq object
#' @param cutoff pvalue cutoff
#' @param ... other parameters passed to \code{\link{DESeq2::results}}
#'
#' @export
#'
get_sig_taxa = function(dds, physeq, cutoff = 0.05, ...) {

  res = results(dds, ...) %>% biobroom::tidy() %>% filter(p.adjusted < cutoff) %>%
    dplyr::rename(OTU = gene)
  tax = tax_table(physeq) %>% data.frame() %>% tibble::rownames_to_column("OTU")

  res = res %>% left_join(tax, by = "OTU")

  return(res)

}


#' Summarize taxonomy
#'
#' Summarize a specified taxonomy within a higher level, usually phylum
#'
#' @param physeq a valid phyloseq object
#' @param grouping_tax the higher level taxonomy to group within
#' @param summary_tax the lower level taxonomy to summarize
#' @param filter numeric, filter out members of the summary taxonomy below this fraction
#'
#' @export
#'
summarize_taxonomy = function(physeq, grouping_tax = "Phylum", summary_tax = "Genus", filter = .01) {
  physeq %>% psmelt() %>% group_by_("OTU", grouping_tax, summary_tax) %>%
    summarise(Abundance = median(Abundance)) %>%
    group_by_(grouping_tax, summary_tax) %>%
    summarise(tot = sum(Abundance)) %>%
    mutate(rel = tot / sum(tot)) %>%
    filter(rel > filter) %>%
    arrange_(grouping_tax, ~desc(rel))
}

#' Subset taxa 2
#'
#' A programming safe version of \code{\link{phyloseq::subset_taxa}}
#'
#' @param physeq a valid phyloseq object
#' @param rank the taxonomic rank of the taxa to subset
#' @param taxa the taxa name to subset to
#'
#' @export
#'
subset_taxa2 = function(physeq, rank, taxa)
{
    if (is.null(tax_table(physeq))) {
        cat("Nothing subset. No taxonomyTable in physeq.\n")
        return(physeq)
    }
    else {
        oldMA <- as(tax_table(physeq), "matrix")
        oldDF <- data.frame(oldMA)
        filtDF = oldDF %>% tbl_df() %>% add_rownames("OTU") %>%
          filter_(interp(~r == taxa , r = as.name(rank)))
        newDF = filtDF[-1]
        rownames(newDF) = filtDF$OTU
        newMA <- as(newDF, "matrix")
        if (inherits(physeq, "taxonomyTable")) {
            return(tax_table(newMA))
        }
        else {
            tax_table(physeq) <- tax_table(newMA)
            return(physeq)
        }
    }
}


#' Plot OTUs
#'
#' Makes a box plot from a list of OTUs
#'
#' @param physeq a valid phyloseq object
#' @param otus a character vector of OTU ids, present in physeq.  Can be "all" to plot all OTUs
#' @param xaxis character, the sample column to plot on the x-axis
#' @param fill character, the sample column to use for colouring the boxes
#' @param labeller character, taxonomy rank to label the OTUs with
#' @param glom character, the taxonomic rank to glom at
#' @param dds valid DESeq2 object.  If present the normalized count data will be used instead of the raw count value.
#' @param justDf logical, should just the data be returned instead of plotting
#'
#' @export
#'
plot_OTUs = function(physeq, otus, xaxis, fill, labeller = "Genus", glom = NULL, dds = NULL, justDf = FALSE) {

  if ("all" %in% otus) {
    subset_data = physeq
  } else {
    subset_data = prune_taxa(otus, physeq)
    names(otus) = otus
  }

  if (!is.null(glom)) {
    subset_data = tax_glom(subset_data, glom)
  }

  subset_data = subset_data %>% psmelt()

  if (!is.null(dds)) {
    count_data = lapply(otus, function(x) plotCounts(dds, x, c(xaxis, fill), returnData = TRUE) %>%
    											tibble::rownames_to_column("Sample")) %>%
      bind_rows(.id = 'OTU')

    subset_data = subset_data %>% dplyr::select(-Abundance) %>%
     left_join(count_data) %>% dplyr::rename(Abundance = count)

  }

   subset_data  = subset_data %>%
    mutate(Log_Abundance = log(Abundance + 1))

   if (justDf) return(subset_data)

   p = subset_data %>%
    ggplot(aes_string(x = xaxis, y = "Log_Abundance", fill = fill)) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Set1") +
      labs(y = "Log Abundance") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

   facet_names = as.character(subset_data[[labeller]])

  if (length(facet_names) == 0) {
    warning("No taxonomic rank found for the requested labeller")
    p = p + facet_wrap(~OTU, scales = "free_y")
  } else {
    facet_names = str_replace_na(facet_names, "Unassigned")
    names(facet_names) = subset_data$OTU
    p = p + facet_wrap(~OTU, scales = "free_y", labeller = labeller(OTU = facet_names))
  }

  return(p)

}


#' Plot genus by phylum
#'
#' @export
plot_genus_by_phylum = function(physeq, phylum, filter = NULL, facet_by_fam = FALSE, xlabel = NULL) {

	tax = as.data.frame(tax_table(physeq))
  otus = rownames(tax[tax$Phylum == phylum & !is.na(tax$Phylum), ])
	phylum_physeq = prune_taxa(otus, physeq)

	if (is.numeric(filter)) {
		topff = filterfun_sample(topf(filter))
		phylum_physeq = prune_taxa(genefilter_sample(phylum_physeq, topff), phylum_physeq)
	}

	gen = get_taxa_unique(phylum_physeq, "Genus")
	gen_cols  = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(gen))
	p = phylum_physeq %>% tax_glom("Genus") %>% psmelt() %>%
		ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
			geom_bar_interactive(aes(tooltip = Genus, data_id = Genus), stat = "identity") +
			scale_fill_manual(values = gen_cols)


	if (facet_by_fam) {
		p = p + facet_wrap(~Family)
	}

	if (!is.null(xlabel)) {
		xlab = setNames(sample_data(phylum_physeq)[[xlabel]], rownames(sample_data(phylum_physeq)))
		p = p + scale_x_discrete(label = xlab)
	}

	return(p)
}


#' Plot number of sequences
#'
#' @export
plot_num_seqs = function(physeq) {

	num_seqs_gg = data_frame(num_seqs = sample_sums(physeq), Sample = names(sample_sums(physeq))) %>%
		mutate(pretty_num = format(num_seqs, big.mark = ",")) %>%
		ggplot(aes(x = reorder(Sample, num_seqs), y = num_seqs)) +
		geom_point_interactive(aes(tooltip = pretty_num, data_id = pretty_num), size = 4) +
		labs(x = "Sample", y = "Number of sequences") +
		scale_y_continuous(breaks = seq(0,
					max(sample_sums(data)) +  0.1*max(sample_sums(data)), by = 10000), labels = scales::comma) +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

	return(num_seqs_gg)

}


#' Plot bar 2
#'
#' An alternate version of \code{\link{phyloseq::plot_bar}}
#'
#'
#' @export
plot_bar2 = function(physeq, rank = "Phylum", x = "Sample", xlabs = NULL) {

	levels = phyloseq::get_taxa_unique(physeq, rank)
	cols  = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(levels))

	plot_data = physeq %>% tax_glom(rank) %>% psmelt()

	p = plot_data %>%
		ggplot(aes_string(x = x, y = "Abundance", fill = rank)) +
		geom_bar_interactive(aes_string(tooltip = rank, data_id = rank), stat = "identity") +
		theme_bw() +
	 	scale_fill_manual(values = cols) +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

	if (!is.null(xlabs)) {
		xlabs = setNames(plot_data[[xlabs]], plot_data$Sample)
	 	p = p + scale_x_discrete(labels = xlabs)
	}

	return(p)
}

#' Convert phyloseq otu table to vegan otu table
#'
#' @export
veganotu = function(physeq) {

  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}


#' Make ordination plots
#'
#' Function to quickly plot a variety of oridation methods Uses
#' \code{\link{plot_ordination}} from phyloseq
#'
#' @param physeq (Required). A phyloseq object
#' @param ord_methods A character vector containing the names of the ordinations
#'   to plot. Default is \code{c("NMDS", "PCoA")}.
#' @param distance A character string. Default is "bray".  Name of the distance
#'   method to use.  Value is passed to \code{\link{ordinate}}.
#' @param colour_var A character string containing the name of the variable to
#'   use for colouring the plot.  Must be a valid column of
#'   \code{sample_data(physeq)}
#'
#' @return Returns a \code{ggplot} object with the plots of interest
#'
#' @export
#'
make_ordination_plots = function(physeq, ord_methods = c("NMDS", "PCoA"), distance = "bray",
                                 colour_var = NULL){

  ords = lapply(ord_methods, function(x) {
    ord = ordinate(physeq = physeq, method = x, distance = distance)
    ord_data = plot_ordination(physeq, ord, justDF = TRUE)
    ord_data$Method = x
    colnames(ord_data)[1:2] = c("Axis1", "Axis2")
    return(ord_data)
  })
  all_data = do.call(rbind, ords)
  ggplot(all_data, aes(x = Axis1, y = Axis2)) +
    geom_point(aes_string(colour = colour_var), size = 3) +
    facet_wrap(~Method, scales = "free") +
    labs(x = "Axis 1", y = "Axis 2")
}



#' Get a dataframe with richness estimates
#'
#' Wrapper around \code{\link{estimate_richness}} to get an easier to use dataframe
#' @param physeq A phyloseq object
#' @param ... Arguments passed to \code{estimate_richness}
#' @return A dataframe with richness estimates and sample data
#' @export
get_richness = function(physeq, ...) {
  rich = estimate_richness(physeq, ...)
  rich %<>% add_rownames("SampleID") %>%
    left_join(sample_data(physeq))
  return(rich)
}


#' A dplyr version of psmelt
#'
#' borrowed from http://chuckpr.github.io/blog/melt.html
#'
#'
#' @export
psmelt_dplyr = function(physeq) {
    sd = data.frame(sample_data(physeq))
    tt = data.frame(tax_table(physeq)) %>% add_rownames("OTU")
    otu_tab = data.frame(otu_table(physeq), check.names = FALSE) %>% add_rownames("OTU")
    otu_tab %>%
        left_join(tt) %>%
        gather_("SampleID", "Abundance", setdiff(colnames(otu_tab), "OTU")) %>%
        left_join(sd)
}

