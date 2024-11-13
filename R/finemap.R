#' Perform LD clumping on a dataframe based on p-value threshold
#'
#' This function performs linkage disequilibrium (LD) clumping on a given dataframe of genetic association results. The function filters the dataframe based on the provided p-value threshold and then applies the LD clumping for different values of `r2`. If clumping fails, the function returns a dataframe with `NA` values in the `SNP` and `test` columns.
#'
#' @param df A data frame containing genetic association results with at least the following columns: `pval` (p-value) and `rsid` (SNP identifier).
#' @param bfile A character string specifying the prefix for the PLINK binary files (e.g., `.bed`, `.bim`, `.fam`).
#' @param pval_threshold A numeric value specifying the p-value threshold for filtering SNPs. Default is `5E-8`.
#' @param clump_r2 A numeric vector specifying the r2 threshold(s) for LD clumping. Default is `0.001`. The function will perform clumping for each r2 value.
#'
#' @return A tibble with columns `SNP`, `test`, and `LD`, where:
#' - `SNP`: the SNPs selected after clumping.
#' - `test`: a character string indicating the test, always "p_ld".
#' - `LD`: the r2 threshold used for clumping.
#'
#' @examples
#' # Example usage
#' \dontrun{
#' result <- finemap_pval_LD(df = my_data, bfile = "path_to_bfile", pval_threshold = 5E-8, clump_r2 = c(0.01, 0.1))
#' }
#' @export
finemap_pval_LD <- function(df, bfile, pval_threshold = 5E-8, clump_r2 = 0.001) {
  # Allow clump_r2 to be a vector and map over each value
  purrr::map_dfr(clump_r2, function(r2) {
    finemap_p_ld <- df %>%
      as.data.frame() %>%
      dplyr::filter(pval <= pval_threshold)

    # Clumping step with error handling
    finemap_p_ld <- tryCatch({
      # Modify column names and perform LD clumping
      ieugwasr::ld_clump(
        dat = finemap_p_ld,
        clump_kb = 10000,
        clump_r2 = r2,
        clump_p = pval_threshold,
        plink_bin = genetics.binaRies::get_plink_binary(),
        bfile = bfile
      )
    }, error = function(e) {
      cat("Clumping failed. Error:", conditionMessage(e), "\n")
      return(NULL)  # Return NULL if clumping fails
    })

    # If clumping fails or returns empty, return a dataframe with NA in 'snp' and 'test'
    if (is.null(finemap_p_ld) || nrow(finemap_p_ld) == 0) {
      tibble::tibble(SNP = NA, test = "p_ld", LD = r2)
    } else {
      # Create the dataframe with snp, test, and LD columns
      tibble::tibble(SNP = finemap_p_ld$rsid, test = "p_ld", LD = r2)
    }
  })
}

#' Create a Table of Credible Sets from SusieR Model
#'
#' This function processes a `susieR` model object to extract the SNPs and their corresponding Posterior Inclusion Probabilities (PIP) for each credible set. It also includes a list of other SNPs in the same credible set.
#'
#' @param model A `susieR` model object, which includes the following components:
#' \itemize{
#'   \item \code{sets} A list containing credible sets, where each set is represented by indices of SNPs.
#'   \item \code{X_column_scale_factors} A vector of SNP identifiers (e.g., SNP names).
#'   \item \code{pip} A named vector of Posterior Inclusion Probability (PIP) values, indexed by SNP names.
#' }
#' @return A tibble (data frame) containing the following columns:
#' \itemize{
#'   \item \code{SNP} The SNP identifiers for each SNP in the credible set.
#'   \item \code{PIP} The Posterior Inclusion Probability for each SNP.
#'   \item \code{cs_snps} A string listing the other SNPs in the same credible set (NA if only one SNP).
#'   \item \code{test} The label "susie", indicating the test used for generating the table.
#' }
#' @examples
#' \dontrun{
#'   result <- susieR_cs_table(susie_model)
#'   print(result)
#' }
#' @export
susieR_cs_table <- function(model) {
  # Extract the credible set list from the model object
  cs_list <- model$sets$cs

  # Extract SNP IDs and PIP values for each credible set
  credible_sets <- purrr::map(cs_list, ~ {
    # Extract the SNP indices for each credible set
    snp_indices <- .x

    # Get the corresponding SNP IDs (assumes that the SNP IDs are the names of the vectors, based on your structure)
    snp_ids <- names(model$X_column_scale_factors)[snp_indices]

    # Get the corresponding PIP values for these SNPs
    pip_values <- model$pip[snp_ids]

    # If there is more than one SNP in the set, list the other SNPs
    other_snps <- if (length(snp_ids) > 1) {
      # For each SNP, list other SNPs in the same set
      purrr::map(snp_ids, ~ paste(setdiff(snp_ids, .x), collapse = ", "))
    } else {
      # If only one SNP, set to NA
      rep(NA, length(snp_ids))
    }

    tibble::tibble(
      SNP = snp_ids,
      PIP = pip_values,
      cs_snps = unlist(other_snps),
      test = "susie"  # Adding the 'test' column with value "susie"
    )
  })

  # Combine the list of tibbles into one dataframe
  result_df <- dplyr::bind_rows(credible_sets)

  return(result_df)
}
