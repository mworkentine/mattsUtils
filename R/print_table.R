# ----------------------------------
#  print a nice datatable
# ----------------------------------

#' Print table
#'
#' Prints a nicely formatted DT::datatable with an option to download the excel file
#'
#' @param table dataframe
#' @param fname filename to use for excel download
#' @param pageLength integer, how many rows to display by default
#' @param digits number of signficant digits to show
#'
#' @return datatable
#'
#' @export
print_table = function(table, fname = "table", pageLength = 25L, digits = 3) {
  DT::datatable(
    dplyr::mutate_if(table, is.numeric, prettyNum, digits = digits),
    extensions = 'Buttons', class = "compact stripe",
    options = list(
      dom = 'Blfrtip',
      buttons = list(
        list(extend = 'excel',
             filename = fname,
             text = "Download Excel")
      ),
      pageLength = pageLength,
      lengthMenu = c(10, 25, 50, 100, 150)
    )
  )
}
