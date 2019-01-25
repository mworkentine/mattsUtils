# ----------------------------------
#  ggheatmap
# ----------------------------------

#' Cluster Order
#'
#' Order columns in a dataframe by a value using hierchical clustering.  Used primarily when
#' making a heatmap style plot so the function takes two columns, one intended to the be rows
#' of the heatmap and the other intended to be the column.
#'
#' @param x dataframe
#' @param row bare name of column to be the rows
#' @param column bare name of column to be the columns
#' @param value bare name of column to order by
#' @param dist_method distance method, passed to stats::dist
#' @param hclust_method hierchical clustering method, passed to stats::hclust
#'
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#'
#' @export
cluster_order = function(x, row, column, value,
												 dist_method = "euclidean", hclust_method = "complete") {

	row = rlang::enquo(row)
	column = rlang::enquo(column)
	value = rlang::enquo(value)

	form = stats::as.formula(paste(quo_name(row), quo_name(column), sep = " ~ "))
	rmat = reshape2::acast(x, form, value.var = quo_name(value), fill = 0)
	cmat = t(rmat)

	rd = stats::dist(rmat, method = dist_method)
	cd = stats::dist(cmat, method = dist_method)
	rh = stats::hclust(rd, method = hclust_method)
	ch = stats::hclust(cd, method = hclust_method)

	rcld = tibble::tibble(!!row := rh$labels[rh$order], row_order = seq(1, length(rh$labels)))
	ccld = tibble(!!column := ch$labels[ch$order], column_order = seq(1, length(ch$labels)))

	x %>%
		dplyr::inner_join(rcld) %>%
		dplyr::inner_join(ccld) %>%
		  dplyr::mutate(
		 	!!row := forcats::fct_reorder(!!row, row_order),
		 	!!column := forcats::fct_reorder(!!column, column_order)
		 ) %>%
		dplyr::select(-row_order, -column_order)
}

#' Make a ggplot heatmap
#'
#' Function to create a ggplot heatmap where the rows and columns are ordered with
#' heirarchical clustering, as in a traditional heatmap
#'
#' @param x dataframe
#' @param row bare name of column to be the rows
#' @param column bare name of column to be the columns
#' @param value bare name of column with values to fill the cells and order the rows and columns
#' @param dist_method distance method, passed to stats::dist
#' @param hclust_method hierchical clustering method, passed to stats::hclust
#' @param scico_pal colour palette to use from the package scico
#'
#' @export
#'
ggheatmap = function(x, row, column, value, scale = c("row", "column"), filter_empty = TRUE,
                     dist_method = "euclidean", hclust_method = "complete", scico_pal = "vik") {

	scale = match.arg(scale)

	row = rlang::enquo(row)
	column = rlang::enquo(column)
	value = rlang::enquo(value)


	if (filter_empty) {
		x = x %>%
			dplyr::group_by(!!row) %>% dplyr::filter(sum(!!value) > 0) %>% dplyr::ungroup() %>%
			dplyr::group_by(!!column) %>% dplyr::filter(sum(!!value) > 0) %>% dplyr::ungroup()
	}

	if (scale == "row") {
		pdata = x %>%
			dplyr::group_by(!!row) %>%
			dplyr::mutate(!!value := scale(!!value)) %>%
			dplyr::ungroup()
	} else if (scale == "column") {
		pdata = x %>%
			dplyr::group_by(!!column) %>%
			dplyr::mutate(!!value := scale(!!value)) %>%
			dplyr::ungroup()
	}


	top = max(abs(pdata[[rlang::quo_name(value)]]), na.rm = TRUE)

	pdata %>%
		cluster_order(!!row, !!column, !!value, dist_method, hclust_method) %>%
		ggplot2::ggplot(ggplot2::aes(x = !!column, y = !!row, fill = !!value)) +
		ggplot2::geom_tile(colour = "grey50") +
		scico::scale_fill_scico(palette = scico_pal, limits = c(-top, top))
}
