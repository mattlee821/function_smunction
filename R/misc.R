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

#' Retrieve Batched Data from BioMart with Filters and Attribute Chunks
#'
#' This function retrieves data from BioMart in batches by splitting attributes into smaller chunks.
#' It fetches data using specified filters and joins the results sequentially, with optional filtering for protein-coding genes
#' and standard chromosomes.
#'
#' @param mart A `Mart` object created with `useMart`, representing the BioMart dataset connection.
#' @param attributes A character vector of attributes to retrieve from BioMart.
#' @param filters A single attribute (string) used as a filter in the BioMart query, often the unique identifier.
#' @param values A vector of values to match with the specified `filters` attribute.
#' @param chunk_size An integer defining the number of attributes per batch request.
#' @return A data frame containing the retrieved and filtered data, with each row representing a unique entity
#' and columns corresponding to the requested attributes.
#'
#' @details This function processes the retrieval in chunks, making sequential calls to BioMart and avoiding large single queries.
#' The resulting data is filtered to include only protein-coding genes and standard chromosomes (1–22, X, Y).
#' If an attribute specified as `filters` is also in `attributes`, it will be removed from `attributes` to prevent duplication.
#'
#' @examples
#' \dontrun{
#' library(biomaRt)
#' mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#' attributes <- c("chromosome_name", "start_position", "end_position", "hgnc_symbol")
#' filters <- "hgnc_symbol"
#' values <- c("BRCA1", "TP53", "EGFR")
#' chunk_size <- 2
#' result <- biomaRt_getBM_batch(mart, attributes, filters, values, chunk_size)
#' }
#'
#' @export
biomaRt_getBM_batch <- function(mart, attributes, filters, values, chunk_size) {
  # Remove `filters` from `attributes` if it's already present
  if (filters %in% attributes) {
    attributes <- setdiff(attributes, filters)
  }

  # Split attributes into chunks, each including the filters column
  attribute_chunks <- lapply(
    split(attributes, ceiling(seq_along(attributes) / chunk_size)),
    function(chunk) c(filters, chunk)
  )

  # Loop over each chunk and fetch data, storing results in a list
  result_list <- list()
  for (i in seq_along(attribute_chunks)) {
    result <- tryCatch(
      biomaRt::getBM(
        attributes = attribute_chunks[[i]],
        filters = filters,
        values = values,
        mart = mart
      ),
      error = function(e) NULL # If a query fails, return NULL for this chunk
    )

    # Only add result to the list if it's not NULL and contains the identifier column
    if (!is.null(result) && filters %in% colnames(result)) {
      result_list[[i]] <- result
    }
  }

  # Sequentially join each data frame in result_list
  sequential_join <- function(result_list, filters, values) {
    # Start with a data frame containing only the 'values' column
    final_result <- data.frame(setNames(list(values), filters))

    # Sequentially join each data frame in result_list
    for (i in seq_along(result_list)) {
      final_result <- dplyr::full_join(final_result, result_list[[i]], by = filters)

      # Free up memory after each join
      gc()
    }

    return(final_result)
  }

  # Perform sequential join on result_list
  map <- sequential_join(result_list, filters = filters, values = values)

  # Apply filtering on the resulting data frame
  map <- map %>%
    dplyr::filter(gene_biotype == "protein_coding") %>%
    dplyr::filter(chromosome_name %in% c(as.character(1:22), "X", "Y", "")) %>%
    dplyr::filter(uniprot_gn_symbol == entrezgene_accession |
                    uniprot_gn_symbol == hgnc_symbol |
                    uniprot_gn_symbol == external_gene_name) %>%
    dplyr::rename(CHR = chromosome_name) %>%
    distinct()

  return(map)
}

