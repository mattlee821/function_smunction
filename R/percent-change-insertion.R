#' @title Percent Change and Insertion
#' @description This function calculates the percent change between corresponding columns in two dataframes and inserts these changes as new columns in the second dataframe. It assumes that the first column in each dataframe is an unused ID column.
#' @param data1 The original dataframe.
#' @param data2 The second dataframe to which the percent change columns will be added.
#' @param col_name The name of the column (as a character string) to calculate the percent change for.
#' @returns A dataframe with the new percent change columns added to `data2`.
#' @note Both `data1` and `data2` must have the same structure and number of rows. Column 1 should be excluded from calculations.
#' @export
calculate_and_insert_percent_change <- function(data1, data2, col_name) {
  new_col_name <- paste0(col_name, "_change")
  percent_change <- -((data2[[col_name]] - data1[[col_name]]) / data1[[col_name]]) * 100
  data2[[new_col_name]] <- round(percent_change, 2) # Round to 2 decimal places
  return(data2)
}
