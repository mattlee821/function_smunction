#' Perform Imputation Test
#'
#' This function performs an imputation test using the specified qmd file.
#' Your data should have samples as rows and features as columns. There should
#' be no extra columns.
#'
#' @param data your data variable e.g., mtcars
#' @param output_dir where you want to save the output
#' @param subtitle what you want to call this imputation test
#'
#' @return The result of the imputation test.
#'
#' @examples
#' \dontrun{
#' imputation_test(data = mtcars,
#'                 output_dir = "project1/data/",
#'                 subtitle = "imputation testing for mtcars")
#' }
#'
#' @export
#' @importFrom rmarkdown render
#' @importFrom here here
imputation_test <- function(data, output_dir, subtitle) {
  rmarkdown::render(input = paste0(system.file(package = "functions"), "/rmd/imputation_test.qmd"),
                    output_format = "html_document",
                    output_dir = output_dir,
                    params = list(data = data, subtitle = subtitle))
}

#' Perform Imputation Test
#'
#' This function performs an imputation test using the specified qmd file.
#' Your data should have samples as rows and features as columns. There should
#' be no extra columns. This function is for testing.
#'
#' @param rmd_location Location of the Rmd file within the package.
#' @param data         Data to be passed to the Rmd file.
#' @param output_dir   Output directory for the rendered HTML file.
#' @param subtitle     Subtitle to be used in the rendered HTML file.
#'
#' @return The result of the imputation test.
#'
#' @examples
#' \dontrun{
#' imputation_test_package(rmd_location = "rmd/imputation_test.qmd",
#'                         data = mtcars,
#'                         output_dir = "project1/data/",
#'                         subtitle = "imputation testing for mtcars")
#' }
#'
#' @export
#' @importFrom rmarkdown render
#' @importFrom here here
imputation_test_package <- function(rmd_location, data, output_dir, subtitle) {
  rmd_file <- here::here(rmd_location)

  if (!file.exists(rmd_file)) {
    stop("Rmd file not found:", rmd_file)
  }

  rmarkdown::render(input = rmd_file,
                    output_format = "html_document",
                    output_dir = output_dir,
                    params = list(data = data, subtitle = subtitle))
}
