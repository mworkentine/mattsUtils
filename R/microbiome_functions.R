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
    ggvis(x = x_axis_var, y = ~Abundance, fill = taxa_level) %>%
    layer_bars()

  if (interactive) {
    plot = plot %>% add_tooltip(all_values, "hover")
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

  plot = ord_data %>% ggvis(x = x, y = y, strokeWidth := 0, fill = colour) %>%
    layer_points()

  if (interactive) {
    plot = plot %>% add_tooltip(all_values, "hover")
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








