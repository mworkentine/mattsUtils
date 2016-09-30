# --- Some useful functions for microbiome analysis

#' ggvis taxa plot
#'
#' @export
ggvis_taxa_plot = function(physeq, x_axis_var, taxa_level, interactive=FALSE) {
  # this is a function that takes a phyloseq object and creates an interactive plot
  # that can be viewed in RStudio or in your browser

  if (!(x_axis_var %in% sample_variables(physeq))) {
    stop("Couldn't find sample variable in physeq object!")
  }

  if (!(taxa_level %in% rank_names(physeq))) {
    stop("Couldn't find taxa level in physeq object!")
  }


  # ---to format the tooltip correctly
  all_values <- function(x) {
    if(is.null(x)) return(NULL)
    show = x
    show$xmin_ = NULL
    show$xmax_ = NULL
    show$stack_upr_ = NULL
    show$stack_lwr_ = NULL
    show$Percent = round(100 * (x$stack_upr_ - x$stack_lwr_), 2)
    paste0(names(show), ": ", format(show), collapse = "<br />")
  }

  # --- need these as formulas for ggvis
  taxa_level = as.formula(paste0("~", taxa_level))
  x_axis_var = as.formula(paste0("~", x_axis_var))

  data_melt = psmelt(physeq)

  plot = data_melt %>% group_by_(taxa_level) %>%
    ggvis::ggvis(x = x_axis_var, y = ~Abundance, fill = taxa_level) %>%
    ggvis::layer_bars()

  if (interactive) {
    plot = plot %>% ggvis::add_tooltip(all_values, "hover")
  }

  plot
}

#' ggvis ordination plot
#'
#' @export
ggvis_ord_plot = function(physeq, ord_obj, colour, interactive=FALSE){

  all_values <- function(x) {
    if(is.null(x)) return(NULL)
    paste0(names(x), ": ", format(x), collapse = "<br />")
    }

  # extract the scores from the ordination
  ord_data = plot_ordination(physeq, ord_obj, justDF = TRUE)

  # format for ggvis
  x = as.formula(paste0("~", colnames(ord_data)[1]))
  y = as.formula(paste0("~", colnames(ord_data)[2]))
  colour = as.formula(paste0("~", as.character(colour)))

  plot = ord_data %>% ggvis::ggvis(x = x, y = y, strokeWidth := 0, fill = colour) %>%
    ggvis::layer_points()

  if (interactive) {
    plot = plot %>% ggvis::add_tooltip(all_values, "hover")
  }

  plot
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
make_ordination_plots = function(physeq, ord_methods = c("NMDS", "PCoA"), distance = "bray", colour_var = NULL){

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
    theme_bw(base_size = 16) + labs(x = "Axis 1", y = "Axis 2")
}



#' Expand ggthemr swatch
#'
#' Utility function to interpolate a ggthemr swatch to create more colours. Uses the swatch from the
#' currently activated theme.
#'
#' @param num_cols integer The number of colours you would like to create
#'
#' @usage scale_colour_manual(values = expand_swatch(num_cols))
#'
#' @return A vector of colours
#' @export
expand_swatch = function(num_cols){
  c(swatch()[1], colorRampPalette(swatch()[-1])(num_cols))
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


#' Get significant taxa
#'
#' Get significantly differentially abundant taxa from a DESeq results object
#'
#' @param res DESeq results object obtained from calling \code{results} on a DESeq object
#' @param alpha P-value cutoff
#' @param physeq The original phyloseq object
#' @return A dataframe with the significant taxa
#' @export
get_sig_taxa = function(res, alpha, physeq){
	if(missing(alpha)){
		stop("Please provide cutoff value")
	}
	if(missing(res)){
		stop("Please provide DESeq results object")
	}
	if(!inherits(res, "DESeqResults")){
		stop("Not a valid DESeq results object")
	}
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
  return(sigtab)
}


#' Plot log-fold changes
#'
#' Plot log-fold changes for taxa in a data frame that contains DESeq results
#' Can be used with output of \code{get_sig_taxa}
#' @param sigs_obj Dataframe with DESeq results
#' @param order string How to order the plot.  One of "lfc" or "Phylum".  Defaults to "lfc"
#' @param xtitle string Title of the x-axis. Defaults to "".
#' @param colour string Variable to colour the points by.  Must be a column in sigs_obj. Defaults to "Phylum"
#'
#' @return A ggplot object
#'
#' @export
lfc_plot = function(sigs_obj, order = "lfc", xtitle = "", colour = "Phylum"){
	if(!("OTU" %in% colnames(sigs_obj))){
		sigs_obj = sigs_obj %>% add_rownames("OTU")
		}
	plot_data = sigs_obj %>%
		filter(!is.na(Family)) %>%
		mutate_each(funs(as.character), Phylum, Class, Order, Family, Genus, Species) %>%
		arrange(desc(Phylum)) %>%
		mutate(label = paste(Family, Genus, Species, sep = " :: "))

	meanByLabel = plot_data %>% group_by(label) %>% summarise(meanlfc = mean(log2FoldChange))

	if(order == "Phylum"){
		plot_data = plot_data %>% mutate(label = factor(label, levels = unique(label)))
	} else {
		plot_data = plot_data %>%
			mutate(label = factor(label, levels = meanByLabel$label[order(meanByLabel$meanlfc)]))
	}

	plot_data %>%
		ggplot(aes_string(y = "label", x = "log2FoldChange", colour = colour)) +
			geom_point(size = 3) + geom_vline(x = 0) +
			theme(axis.text.y = element_text(hjust = 0)) +
		labs(x = xtitle, y = "")
}





