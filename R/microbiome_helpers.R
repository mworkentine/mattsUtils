
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
#' @import biobroom
#' @import broom
#'
get_sig_taxa = function(dds, physeq, cutoff = 0.05, ...) {

  res = results(dds, ...) %>% tidy() %>% filter(p.adjusted < cutoff) %>%
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
#' @return a data frame
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
subset_taxa_safe = function(physeq, rank, taxa) {
	tx = as.data.frame(tax_table(physeq))
	otus = rownames(tx[tx[[rank]] == taxa, ])
	phyloseq::prune_taxa(otus, physeq)
}

#' Subset samples 2
#'
#' A programming safe version of \code{\link{phyloseq::subset_samples}}
#'
#' @param physeq a valid phyloseq object
#' @param var variable (column name) of sample data
#' @param value filter to keep only this value of the variable
#'
#' @export
#'
subset_samples_safe = function(physeq, var, value) {
  keeps = rownames(phylseq::sample_data(physeq)[phyloseq::sample_data(physeq)[[var]] == value,])
  phyloseq::prune_samples(keeps, physeq)
}


#' Replace counts
#'
#' Replace OTU counts in a phyloseq object with the normalized counts from a DESeq object containing
#' the same OTUS
#'
#' @export
#'
replace_counts = function(physeq, dds) {

  dds_counts = counts(dds, normalized = TRUE)
  if (!identical(taxa_names(physeq), rownames(dds_counts))) {
    stop("OTU ids don't match")
  }
  otu_table(physeq) = otu_table(dds_counts, taxa_are_rows = TRUE)
  return(physeq)

}


#' Plot OTUs
#'
#' Makes a box plot from a list of OTUs
#'
#' @param physeq a valid phyloseq object with raw OTU counts that have not been transformed.
#' @param otus a character vector of OTU ids, present in physeq.  Can be "all" to plot all OTUs
#' @param xaxis character, the sample column to plot on the x-axis
#' @param fill character, the sample column to use for colouring the boxes
#' @param labeller character, taxonomy rank to label the OTUs with.  Can also be "label" which will
#'   label with genus if available and family if genus is unassigned.  In this case a prefix for the
#'   taxonomic rank is addes as well
#' @param scales fixed or free scales, passed to facet_wrap
#' @param palette RColorBrewer palette to use
#' @param glom character, the taxonomic rank to glom at
#' @param dds valid DESeq2 object.  If present the normalized count data will be used instead of the
#'   raw count value.
#' @param justDf logical, should just the data be returned instead of plotting
#' @param y_scale character, the abundance scale on the y-axis.  Counts are just the raw OTU counts
#'   unless a DESeq object is provided with the dds argument, in which case the normalized counts
#'   will be used.  log_counts will transform the counts by log2 using a pseudcount of +0.5.
#'   relative will use the relative abundance values instead of counts.
#'   @param nrow number of facet rows, passed to facet_wrap
#'   @param ncol same as nrow but for columns
#'
#' @importFrom stringr str_c
#' @export
#'
plot_OTUs = function(physeq, otus, xaxis, fill, labeller = "Genus", scales = "free_y",
                     palette = "Set1", glom = NULL, dds = NULL, justDf = FALSE,
                     y_scale = c("log_counts", "counts", "relative"), nrow = NULL, ncol = NULL)
  {

  y_scale = match.arg(y_scale)

  # replace with normalized counts if dds is present
  if (!is.null(dds)) {
    physeq = replace_counts(physeq, dds)
  }

  # transform to appropriate scale
  if (y_scale == "relative") {
    physeq = transform_sample_counts(physeq, function(x) x/sum(x))
    y_axis_title = "Relative Abundance"
  } else if (y_scale == "log_counts") {
    physeq = transform_sample_counts(physeq, function(x) log2(x + 0.5))
    y_axis_title = "Log2 Abundance"
  } else {
    y_axis_title = "Abundance"
  }

  if ("all" %in% otus) {
    subset_data = physeq
  } else {
    subset_data = prune_taxa(otus, physeq)
    names(otus) = otus
  }


  if (!is.null(glom)) {
    subset_data = tax_glom(subset_data, glom)
  }

  subset_data = subset_data %>% psmelt() %>%
    mutate(label = case_when(
      is.na(Phylum) ~ "Unassigned",
      is.na(Class) ~ str_c("p:", Phylum),
      is.na(Order) ~ str_c("c", Class),
      is.na(Family) ~ str_c("o:", Order),
      is.na(Genus) ~ str_c("f:", Family),
      TRUE ~ str_c("g:", Genus)
    ))

  if (justDf) return(subset_data)

   p = subset_data %>%
    ggplot(aes_string(x = xaxis, y = "Abundance", fill = fill)) +
      geom_boxplot() +
      scale_fill_brewer(palette = palette) +
      labs(y = y_axis_title)

   facet_names = as.character(subset_data[[labeller]])

  if (length(facet_names) == 0) {
    warning("No taxonomic rank found for the requested labeller")
    p = p + facet_wrap(~OTU, scales = "free_y")
  } else {
    facet_names = stringr::str_replace_na(facet_names, "Unassigned")
    names(facet_names) = subset_data$OTU
    p = p + facet_wrap(~OTU, scales = scales, labeller = labeller(OTU = facet_names),
                       nrow = nrow, ncol = ncol)
  }

  p = p + theme(strip.background = element_blank())

  return(p)

}


#' Plot genus by phylum
#'
#' Bar plot of all the genera within a specified phylum.
#'
#' @param physeq A valid phyloseq object
#' @param phylum character, the phylum to plot
#' @param x character, the variable to plot on the x-axis, should be one of
#'   sample_variables(physeq).  Note that if this is anything other than
#'   "Samples", the default, the mean abundance for each value will be plotted.
#' @param filter numeric, top fraction of abundance to return
#' @param  facet_by_fam logical, facet by family, default is FALSE
#' @param xlabel character, when plotting individual samples the variable to use for labels
#'
#' @export
plot_genus_by_phylum = function(physeq, phylum, x = "Sample", filter = NULL,
                                facet_by_fam = FALSE, xlabel = NULL) {

	tax = as.data.frame(tax_table(physeq))
  otus = rownames(tax[tax$Phylum == phylum & !is.na(tax$Phylum), ])
	phylum_physeq = prune_taxa(otus, physeq)

	if (is.numeric(filter)) {
		topff = filterfun_sample(topf(filter))
		phylum_physeq = prune_taxa(genefilter_sample(phylum_physeq, topff), phylum_physeq)
	}

	p_data = phylum_physeq %>% tax_glom("Genus") %>% psmelt()

	if (x != "Sample") {
	   p_data = p_data %>%
	     group_by_("Phylum", "Genus", x) %>%
	     summarise(Abundance = mean(Abundance))
	}

	gen = get_taxa_unique(phylum_physeq, "Genus")
	gen_cols  = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(gen))
	p =  p_data %>%
		ggplot(aes_string(x = x, y = "Abundance", fill = "Genus")) +
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
		ggiraph::geom_point_interactive(aes(tooltip = pretty_num, data_id = pretty_num), size = 4) +
		labs(x = "Sample", y = "Number of sequences") +
		scale_y_continuous(breaks = seq(0,
					max(sample_sums(physeq)) +  0.1*max(sample_sums(physeq)), by = 10000), labels = scales::comma) +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

	return(num_seqs_gg)

}


#' Make colour palettes for taxa
#'
#' @param rank Taxonomic rank to use
#' @param physeq valid phyloseq object
#' @param palette Either a valid RColorBrewer palette or "iwanthue" to generate
#'   a random palette
#' @export
make_taxa_colours = function(rank, physeq, palette = "iwanthue") {
  taxa = na.omit(get_taxa_unique(physeq, rank))
  ncols = length(taxa)
  if (!palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    pal = hues::iwanthue(ncols, 0, 360, 40, 70, 50, 95, random = TRUE)
  } else {
    palsize = RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    pal = colorRampPalette(RColorBrewer::brewer.pal(palsize, palette))(ncols)
  }
  return(setNames(pal, taxa))
}

#' Plot bar 2
#'
#' An alternate version of \code{\link{phyloseq::plot_bar}}
#'
#' @param physeq valid phyloseq object
#' @param rank character, the taxa rank to plot, must be one of \code{rank_names(physeq)}
#' @param glom logical, should the taxa be aggregated at the level of \code{rank}
#' @param x character, the sample variable to plot on the x-axis
#' @param xlabs character vector, labels to use on the x-axis
#' @param position character, for geom_bar
#' @param palette Either a valid RColorBrewer palette or "iwanthue" to generate
#'   a random palette
#' @param pal_seed Seed for generating palette.  Set to NULL for random
#'
#' @return a ggplot object
#'
#' @export
#'
plot_bar2 = function(physeq, rank = "Phylum", glom = TRUE, x = "Sample",
                     xlabs = NULL, position = "stack", palette = "iwanthue",
                     pal_seed = 5858) {

  if (!is.null(pal_seed)) set.seed(pal_seed)
	cols = make_taxa_colours(rank, physeq, palette)

	if (glom) {
	  plot_data = physeq %>% tax_glom(rank) %>% psmelt()
	} else {
	  plot_data = physeq %>% psmelt()
	}

	p = plot_data %>%
		ggplot(aes_string(x = x, y = "Abundance", fill = rank)) +
		ggiraph::geom_bar_interactive(aes_string(tooltip = rank, data_id = rank),
		                              stat = "identity", position = position) +
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

#' VST Tranform a phyloseq object
#'
#' Uses the \code{\link{getVarianceStabilizedData}} function from DESeq2 to vst
#' transform an OTU table
#'
#' @param physeq a phyloseq object where the OTU table is raw counts
#' @return a phyloseq object with vst normalized OTU table
#' @export
vst_transform = function(physeq, size_factors = NULL) {
  vst = phyloseq_to_deseq2(physeq, ~1)
  if (!is.null(size_factors)) {
    sizeFactors(vst) = size_factors
  } else {
    vst = estimateSizeFactors(vst, type = "poscounts")
  }

  vst = estimateDispersions(vst, fitType = "local")
  vst = getVarianceStabilizedData(vst)
  physeq_vst = physeq
  otu_table(physeq_vst) = otu_table(vst, taxa_are_rows = TRUE)
  return(physeq_vst)
}


split_species = function(string, n = 2) {
  splits = str_split(string, "/", n + 1)
  res = map_if(splits, ~length(.x) > 2, ~.x[1:n]) %>%
    map_chr(str_c, collapse = "/")
  return(res)
}


#' Add taxonomy label
#'
#' add a column to the taxonomy table of a phyloseq object that lists the
#' lowest rank taxonomy assigned to that OTU along with a prefix indicating
#' the taxonomic rank.
#'
#' Example: g:Pseudomonas
#'
#' @param physeq a valid phyloseq object that contains a taxonomy table
#' @param num_species the number of species to retain if more than one are identified
#' @return a phyloseq object with an additional column on the taxonomy
#'         table called "Taxonomy"
#' @export
add_taxonomy_column = function(physeq, num_species = 2) {
  tax_df = as.data.frame(tax_table(physeq)) %>%
    rownames_to_column("OTU") %>%
    mutate(Species = split_species(Species, n = num_species)) %>%
    mutate(Taxonomy =
      case_when(
        is.na(Class)  ~ str_c("o:", Phylum),
        is.na(Order)  ~ str_c("o:", Class),
        is.na(Family)  ~ str_c("o:", Order),
        is.na(Genus)   ~ str_c("f:", Family),
        is.na(Species) ~ str_c("g:", Genus),
        TRUE ~ str_c(Genus, " ", Species)
      )
    )

  tax = as.matrix(tax_df[, -1])
  rownames(tax) = tax_df$OTU
  tax_table(physeq) = tax_table(tax)

  return(physeq)
}


#' Sequence stats
#' @export
get_seq_stats = function(physeq) {
  return(list(
    ntaxa = ntaxa(physeq),
    sum = sum(sample_sums(physeq)),
    nsamples = nsamples(physeq),
    median = median(sample_sums(physeq)),
    sd  = sd(sample_sums(physeq)),
    min = min(sample_sums(physeq)),
    max = max(sample_sums(physeq))
    )
  )
}


