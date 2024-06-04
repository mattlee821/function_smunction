# function to calculate percent change between corresponding cells in two dataframes and insert these as new columns

#' @title precent change and insertion
#' @description You need two data frames of equal size. The function assumes column 1 is an unused ID column.
#' @param data1 your original dataframe
#' @param data2 your second dataframe
#' @param col_name the name of the column to calculate percent change for
#' @returns a dataframe
#' @note
#' for (col_name in colnames(data2)) {
#'   data2 <- calculate_and_insert_percent_change(data1, data2, col_name)}
#' @export
calculate_and_insert_percent_change <- function(data1, data2, col_name) {
  new_col_name <- paste0(col_name, "_change")
  percent_change <- -((data2[[col_name]] - data1[[col_name]]) / data1[[col_name]]) * 100
  data2[[new_col_name]] <- round(percent_change, 2) # Round to 2 decimal places
  return(data2)
}
