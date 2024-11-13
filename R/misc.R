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

#' check the alignment of GWAS to reference panel and flip allele beta where required
#'
#' This function performs an alignment check for GWAS data by comparing the input
#' data (`df`) with a reference dataset (`reference`). It checks if the alleles
#' in the data are aligned and flips the alleles when necessary to ensure
#' consistency with the reference. The function also computes an LD matrix,
#' performs a kriging-based procedure to adjust the z-scores, and generates a
#' series of plots to visualize the alignment.
#'
#' @param df A data frame containing the GWAS summary statistics.
#' The following columns should be present:
#' \itemize{
#'   \item \strong{CHR}: Chromosome number (integer).
#'   \item \strong{POS}: Position of the SNP (integer).
#'   \item \strong{SNP}: SNP identifier (character).
#'   \item \strong{EA}: Effect allele (character).
#'   \item \strong{OA}: Other allele (character).
#'   \item \strong{EAF}: Effect allele frequency (numeric).
#'   \item \strong{BETA}: Effect size estimate (numeric).
#'   \item \strong{SE}: Standard error of the effect size (numeric).
#'   \item \strong{P}: P-value for the association (numeric).
#'   \item \strong{N}: Sample size (integer).
#'   \item \strong{phenotype}: Phenotype identifier (character).
#' }
#'
#' @param reference A data frame containing the reference data for comparison.
#' The following columns should be present:
#' \itemize{
#'   \item \strong{Predictor}: SNP identifier (character).
#'   \item \strong{A1}: Allele 1 (character).
#'   \item \strong{A2}: Allele 2 (character).
#'   \item \strong{A1_Mean}: Mean value for allele 1 (numeric).
#'   \item \strong{MAF}: Minor allele frequency (numeric).
#'   \item \strong{Call_Rate}: Call rate (integer).
#'   \item \strong{Info}: Information score (integer).
#' }
#'
#' @param bfile file path for reference population (built for using 1kG,
#' e.g., /path/EUR/EUR).
#'
#' @return A list containing:
#' \itemize{
#'   \item \strong{plots}: A list of plots, including the alignment plot and
#'   observed vs expected z-score plots.
#'   \item \strong{list_df}: The final data frame after allele flipping and
#'   adjustments.
#'   \item \strong{lambda}: A list of lambda estimates for the adjusted z-scores.
#' }
#'
#' @details
#' The function first merges the GWAS summary statistics with the reference data
#' based on the SNP identifier and flips the effect sizes if the effect allele
#' (EA) does not match allele 1 (A1) in the reference. It then computes an LD
#' matrix using the `ieugwasr::ld_matrix_local` function, which is used in
#' subsequent analyses. The function also runs a kriging procedure using
#' `susieR::kriging_rss` to adjust the z-scores based on the LD matrix, and
#' generates plots comparing the observed and expected z-scores before and
#' after allele flipping. The function will iteratively flip alleles and update
#' the data until no moreallele flips are needed (based on a log likelihood
#' ratio test).
#' @export
alignment_check <- function(df, reference, bfile) {
  cat("# START \n")
  list_plot <- list()
  list_lambda <- list()

  # Data ====
  df <- df %>%
    dplyr::inner_join(reference, by = c("SNP" = "Predictor")) %>%
    dplyr::mutate(BETA = ifelse(EA != A1, -BETA, BETA)) %>%  # Flip beta if EA != A1 to align GWAS with reference
    dplyr::distinct(SNP, .keep_all = TRUE)
  # Check if the df has been created correctly
  if (nrow(df) == 0) {
    stop("Data frame 'df' is empty")
  }
  cat("## DF made \n")

  # LD matrix ====
  ld <- ieugwasr::ld_matrix_local(
    df$SNP,
    with_alleles = FALSE,
    bfile = bfile,
    plink_bin = plinkbinr::get_plink_exe()
  )
  # Check if the LD matrix has been created correctly
  if (is.null(ld) || nrow(ld) == 0) {
    stop("LD matrix is NULL or empty")
  }
  cat("## LD matrix made \n")

  # list_df ====
  list_df <- list(
    CHR = df$CHR,
    POS = df$POS,
    SNP = rownames(ld),
    EA = df$EA,
    OA = df$OA,
    MAF = df$EAF,
    BETA = df$BETA,
    SE = df$SE,
    P = df$P,
    varbeta = df$SE^2,
    N = min(df$N),
    phenotype = unique(df$phenotype),
    Z = df$BETA / df$SE,
    LD = ld,
    type = "quant"
  )
  cat("## list_df made \n")

  # alignment plot ====
  plot_alignment <- list_df %>%
    purrr::set_names(~ case_when(
      . == "SNP" ~ "snp",  # Rename SNP to snp
      . == "P" ~ "p",      # Rename P to p
      . == "BETA" ~ "beta",      # Rename P to p
      TRUE ~ .              # Leave other names unchanged
    )) %>%
    functions::coloc_check_alignment(do_plot = TRUE)
  # Check if the alignment plot is not NULL
  if (is.null(plot_alignment)) {
    stop("Alignment plot is NULL")
  }
  cat("## plot_alignment made \n")
  list_plot[[length(list_plot) + 1]] <- plot_alignment

  # VAR_lambda ====
  VAR_lambda <- susieR::estimate_s_rss(z = list_df$Z,
                                       R = list_df$LD,
                                       n = list_df$N)
  list_lambda[[length(list_lambda) + 1]] <- VAR_lambda
  cat("## S:", VAR_lambda, "\n")
  cat("### larger S indicates inconsistency between Z and LD matrix \n")

  # kriging_rss ====
  res_kriging <- susieR::kriging_rss(
    z = list_df$Z,
    R = list_df$LD,
    n = list_df$N,
    s = VAR_lambda)
  # Check for NULL or empty values in the conditional distribution
  if (is.null(res_kriging$conditional_dist)) {
    stop("kriging_rss conditional distribution is NULL")
  }
  cat("## kriging_rss made \n")

  # Create plot for observed vs expected ====
  plot_observed_expected <- res_kriging$plot +
    ggplot2::labs(title = "Distribution of z-scores") +
    ggplot2::labs(subtitle = paste0("Before allele flipping; S = ", round(VAR_lambda,4)))
  # Append the plot to the list
  list_plot[[length(list_plot) + 1]] <- plot_observed_expected
  cat("## plot_observed_expected made \n")

  # id for flipping ====
  id <- res_kriging$conditional_dist
  id <- which(id$logLR > 2 & abs(id$z) > 2)
  # If no rows to flip, exit the loop and return the plots and list_df
  if (length(id) == 0) {
    cat("No alleles to flip, exiting loop.\n")
    return(list(plots = list_plot, list_df = list_df, S = VAR_lambda))
  }
  cat("## Alleles to flip: ", length(id), "\n")

  # allele flip loop ====
  iteration_counter <- 1  # Initialize the counter for iterations
  repeat {
    # Print iteration number
    cat("## Starting iteration: ", iteration_counter, "\n")

    # id for flipping ====
    id <- susieR::kriging_rss(
      z = list_df$Z,
      R = list_df$LD,
      n = list_df$N,
      s = VAR_lambda)$conditional_dist
    id <- which(id$logLR > 2 & abs(id$z) > 2)
    # If no rows to flip, exit the loop and return the plots and list_df
    if (length(id) == 0) {
      cat("No alleles to flip, exiting loop.\n")
      return(list(plots = list_plot, list_df = list_df, S = list_lambda))
    }
    cat("## Alleles to flip: ", length(id), "\n")

    # flip alleles ====
    df$BETA[id] <- -df$BETA[id]
    cat("## alleles flipped \n")

    # LD matrix ====
    ld <- ieugwasr::ld_matrix_local(
      df$SNP,
      with_alleles = FALSE,
      bfile = bfile,
      plink_bin = plinkbinr::get_plink_exe()
    )
    # Check if the LD matrix has been created correctly
    if (is.null(ld) || nrow(ld) == 0) {
      stop("LD matrix is NULL or empty")
    }
    cat("## LD matrix made \n")

    # list_df ====
    list_df <- list(
      CHR = df$CHR,
      POS = df$POS,
      SNP = rownames(ld),
      EA = df$EA,
      OA = df$OA,
      MAF = df$EAF,
      BETA = df$BETA,
      SE = df$SE,
      P = df$P,
      varbeta = df$SE^2,
      N = min(df$N),
      phenotype = unique(df$phenotype),
      Z = df$BETA / df$SE,
      LD = ld,
      type = "quant"
    )
    cat("## list_df made \n")

    # alignment plot ====
    plot_alignment <- list_df %>%
      purrr::set_names(~ case_when(
        . == "SNP" ~ "snp",  # Rename SNP to snp
        . == "P" ~ "p",      # Rename P to p
        . == "BETA" ~ "beta",      # Rename P to p
        TRUE ~ .              # Leave other names unchanged
      )) %>%
      functions::coloc_check_alignment(do_plot = TRUE)
    # Check if the alignment plot is not NULL
    if (is.null(plot_alignment)) {
      stop("Alignment plot is NULL")
    }
    cat("## plot_alignment made \n")
    list_plot[[length(list_plot) + 1]] <- plot_alignment

    # VAR_lambda ====
    VAR_lambda <- susieR::estimate_s_rss(z = list_df$Z,
                                         R = list_df$LD,
                                         n = list_df$N)
    list_lambda[[length(list_lambda) + 1]] <- VAR_lambda
    cat("## S:", VAR_lambda, "\n")
    cat("### larger S indicates inconsistency between Z and LD matrix \n")

    # kriging_rss ====
    res_kriging <- susieR::kriging_rss(
      z = list_df$Z,
      R = list_df$LD,
      n = list_df$N,
      s = VAR_lambda)
    # Check for NULL or empty values in the conditional distribution
    if (is.null(res_kriging$conditional_dist)) {
      stop("kriging_rss conditional distribution is NULL")
    }
    cat("## kriging_rss made \n")

    # observed vs expected plot ====
    plot_observed_expected <- res_kriging$plot +
      ggplot2::labs(title = "Distribution of z-scores") +
      ggplot2::labs(subtitle = paste0("Allele flipping iteration = ", iteration_counter, "; lambda = ", round(VAR_lambda,4)))
    # Append the plot to the list
    list_plot[[length(list_plot) + 1]] <- plot_observed_expected
    cat("## plot_observed_expected made \n")

    # iteration counter ====
    iteration_counter <- iteration_counter + 1
  }
}
