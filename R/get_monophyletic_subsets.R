#' Get monophyletic subsets
#'
#' Given a phylogenetic tree and a subset of tips in that tree, find all
#' monophyletic groups within that subset
#'
#' @param tree phylo, a phylogentic tree
#' @param tips character vector of tips in the tree
#'
#' @return a list
#' @export
get_monophyletic_subsets = function(tree, tips) {

  if (!inherits(tree, "phylo")) stop("tree must be of class phylo")
  if (!all(tips %in% tree$tip.label)) stop("Provided tips not found in tree")

	# get node ids of tips
	nodes = match(tips, tree$tip.label)

	# MRCA matrix of whole tree, subset to nodes of interest
	mrca_mat = phangorn::mrca.phylo(tree)[nodes, nodes]

	# get the unique set of MRCAs from matrix and subset to internal nodes
	mrca_vals = unique(as.vector(mrca_mat))
	mrca_vals = setdiff(mrca_vals, 1:length(tree$tip.label))

	# find decendents of these
	descen = phangorn::Descendants(tree, mrca_vals)

	# intersect with tips of interest to see if in our list
	# these are our monophyletic groups
	hits = descen[purrr::map_lgl(descen, ~all(.x %in% nodes))]

	# collapse to largest monophyletic groups
	is_subset = purrr::map(hits, function(x) purrr::map_lgl(hits, ~all(x %in% .x) ))
	hits = hits[purrr::map_lgl(is_subset, ~sum(.x) == 1)]

	return(purrr::map(hits, function(x) tree$tip.label[x]))

}
