#' Read File Based on Extension
#'
#' This function reads a file based on its extension and returns the data in a suitable format.
#' It supports reading text files, CSV files, Stata files, Excel files, and SAS files.
#'
#' @param file_path A string specifying the path to the file.
#' @return A data frame or data table containing the data read from the file.
#' @details
#' The function uses different packages to read different file types:
#' - \code{data.table::fread} for reading text and CSV files.
#' - \code{readstata13::read.dta13} for reading Stata files.
#' - \code{readxl::read_xlsx} for reading Excel files.
#' - \code{haven::read_sas} for reading SAS files.
#'
#' @export
read_file <- function(file_path) {
  # Get the file extension
  file_extension <- tools::file_ext(file_path)

  # Read the file based on its extension
  if (file_extension == "txt" || file_extension == "csv") {
    # Use data.table::fread for reading text and CSV files
    data <- data.table::fread(file_path)
  } else if (file_extension == "dta") {
    # Use readstata13::read.dta13 for reading Stata files
    data <- readstata13::read.dta13(file_path)
  } else if (file_extension == "xlsx") {
    # Use readxl::read_xlsx for reading Excel files
    data <- readxl::read_xlsx(file_path)
  } else if (file_extension == "sas7bdat" || file_extension == "xpt") {
    # Use haven::read_sas for reading SAS files
    data <- haven::read_sas(file_path)
  } else {
    # Stop and return an error message if the file extension is not supported
    stop("Unsupported file extension: ", file_extension)
  }

  return(data)
}
