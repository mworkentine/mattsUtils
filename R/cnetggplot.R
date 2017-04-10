#' cnetggplot
#'
#' This function makes a similar plot to the 'cnetplot' in the Bioconductor package 'clusterProfiler but does it with ggplot using the ggraph package
#'
#' @param res A clusterProfiler results object
#' @param dds_res DESeq results cleaned up with biobroom::tidy
#' @param algorithm the layout algorithm for igraph
#' @param exclude any pathways to exclude from the path
#'
#' @import ggraph
#' @import igraph
#' @importFrom viridis scale_color_viridis
#'
#' @export
cnetggplot = function(res, dds_res, algorithm = "kk", exclude = NULL) {
  # get pathway-gene data
  gc = geneInCategory(res)
  if (!is.null(exclude)) {
    gc = gc[!names(gc) %in% exclude]
  }
  resdf = as.data.frame(res)
  y = resdf[resdf$ID %in% names(gc), c("ID", "Description", "pvalue")]
  y$Description = str_wrap(y$Description, 15)
  gc = set_names(gc[as.character(y$ID)], y$Description)

  # vertex attributes
  category_data = y %>% dplyr::select(name = Description, pvalue = pvalue) %>%
    mutate(Category = "Category")
  gene_data = dds_res %>% dplyr::select(name = gene, lfc = estimate) %>%
    mutate(Category = "gene")
  v_data = full_join(category_data, gene_data)

  # convert to graph and create layout
  gdata = gc %>% map_df(~data_frame(to = .), .id = "from")
  v_data = v_data %>% filter(name %in% c(gdata$from, gdata$to))
  g = gdata %>% graph_from_data_frame(directed = FALSE, vertices = v_data)
  layout = create_layout(g, layout = "igraph", algorithm = algorithm)

  ggraph(layout) +
    geom_edge_link0(colour = "grey70", alpha = 0.25) +
    geom_node_point(aes(filter = Category == "Category"), size = 2, colour = "Goldenrod") +
    geom_node_point(aes(filter = Category == "gene", colour = lfc), size = 0.75) +
    geom_node_label(aes(filter = Category == "Category", label = name), size = 2, alpha = 0.75, repel = TRUE) +
    scale_color_viridis("Log2 fold-change") +
    scale_size(trans = "reverse") +
    ggforce::theme_no_axes() +
    theme(panel.border = element_blank())
}
