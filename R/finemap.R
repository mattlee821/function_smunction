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
#' This function processes a `susieR` model object to extract the SNPs and their corresponding Posterior Inclusion Probabilities (PIP) for each credible set. It also includes a list of other SNPs in the same credible set, ensuring that only the SNP with the highest PIP in each credible set receives a label in the `label` column.
#'
#' @param susieR_model A `susieR` model object, which includes the following components:
#' \itemize{
#'   \item \code{sets} A list containing credible sets, where each set is represented by indices of SNPs.
#'   \item \code{X_column_scale_factors} A vector of SNP identifiers (e.g., SNP names).
#'   \item \code{pip} A named vector of Posterior Inclusion Probability (PIP) values, indexed by SNP names.
#' }
#' @param df A data frame containing SNP information, with the following columns:
#' \itemize{
#'   \item \code{SNP} The SNP identifiers.
#'   \item \code{POS} The position of each SNP.
#'   \item \code{P} The p-value associated with each SNP.
#' }
#' @return A tibble (data frame) containing the following columns:
#' \itemize{
#'   \item \code{SNP} The SNP identifiers.
#'   \item \code{POS} The position of each SNP.
#'   \item \code{P} The p-value for each SNP.
#'   \item \code{PIP} The Posterior Inclusion Probability for each SNP.
#'   \item \code{cs_snps} A string listing the other SNPs in the same credible set (NA if only one SNP).
#'   \item \code{cs} The credible set identifier, represented as a factor.
#'   \item \code{test} A string indicating the test used for generating the table (always "susie").
#'   \item \code{label} The SNP identifier for the SNP with the lowest P value within each credible set (NA for other SNPs).
#' }
#' @examples
#' \dontrun{
#'   result <- susieR_cs_table(susie_model, df)
#'   print(result)
#' }
#' @export
susieR_cs_table <- function(susieR_model, df) {
  # Extract the credible sets
  cs_list <- susieR_model$sets$cs
  PIP <- data.frame(SNP = names(susieR_model$pip), PIP = susieR_model$pip)

  table <- data.frame(
    SNP = df$SNP,
    POS = df$POS,
    P = df$P
  ) %>%
    left_join(PIP, by = "SNP")

  # Initialize the `cs` column to store credible set labels
  table$cs <- NA

  # Assign credible set labels (e.g., "L1", "L2") to the `cs` column
  for (i in seq_along(cs_list)) {
    row_indices <- cs_list[[i]]  # Extract row numbers for this list item
    table$cs[row_indices] <- paste0("L", i)  # Assign the label to these rows
  }

  # Initialize the `cs_snps` column to store SNPs in the same credible set
  table$cs_snps <- NA

  # Populate the `cs_snps` column with comma-separated SNPs in the same credible set
  for (i in seq_along(cs_list)) {
    row_indices <- cs_list[[i]]  # Row indices for the current credible set
    snps <- table$SNP[row_indices]  # Extract SNP values for these rows

    # Assign concatenated SNPs excluding the current SNP for each row
    for (row in row_indices) {
      other_snps <- setdiff(snps, table$SNP[row])
      if (length(other_snps) > 0) {
        table$cs_snps[row] <- paste(other_snps, collapse = ", ")
      }
    }
  }

  # Add a label for the SNP with the smallest p-value in each credible set
  table <- table %>%
    dplyr::filter(!is.na(cs)) %>% # Exclude rows where cs is NA
    dplyr::group_by(cs) %>%
    dplyr::mutate(label = if_else(PIP == max(PIP), SNP, NA_character_)) %>%  # Identify SNP with largest PIP
    dplyr::ungroup() %>%
    dplyr::bind_rows(table %>% dplyr::filter(is.na(cs))) %>% # Add back rows with NA cs
    dplyr::mutate(test = "SuSiE")

  # Return the final table
  return(table)
}

#' Create a Table of Credible Sets from Finimom Model
#'
#' This function processes a Finimom model to generate a table summarizing credible sets,
#' their associated SNPs, positions, and posterior inclusion probabilities (PIPs).
#' Additionally, it identifies the SNP with the smallest p-value in each credible set.
#'
#' @param finimom_model A Finimom model object containing credible sets and PIPs.
#' @param df A data frame containing SNP information with columns:
#'   \describe{
#'     \item{SNP}{Character vector of SNP identifiers.}
#'     \item{POS}{Numeric vector of SNP positions.}
#'     \item{P}{Numeric vector of p-values.}
#'   }
#'
#' @return A data frame summarizing the credible sets with the following columns:
#' \itemize{
#'   \item SNP: SNP identifier.
#'   \item POS: SNP position.
#'   \item P: SNP p-value.
#'   \item PIP: Posterior inclusion probability for the SNP.
#'   \item cs: Credible set label (e.g., "L1", "L2").
#'   \item cs_snps: Comma-separated list of other SNPs in the same credible set.
#'   \item label: SNP identifier of the SNP with the smallest p-value in each credible set. We choose P because PIP is equal for SNPs in a CS.
#' }
#'
#' @export
finimom_cs_table <- function(finimom_model, df) {

  # Extract the credible sets from the Finimom model
  cs_list <- finimom_model$sets

  # Create the initial table with SNP details and PIP values
  table <- data.frame(
    SNP = df$SNP,  # SNP identifiers
    POS = df$POS,  # SNP positions
    P = df$P,      # SNP p-values
    PIP = finimom_model$pip  # Posterior inclusion probabilities
  )

  # Initialize the `cs` column to store credible set labels
  table$cs <- NA

  # Assign credible set labels (e.g., "L1", "L2") to the `cs` column
  for (i in seq_along(cs_list)) {
    row_indices <- cs_list[[i]]  # Extract row numbers for this list item
    table$cs[row_indices] <- paste0("L", i)  # Assign the label to these rows
  }

  # Initialize the `cs_snps` column to store SNPs in the same credible set
  table$cs_snps <- NA

  # Populate the `cs_snps` column with comma-separated SNPs in the same credible set
  for (i in seq_along(cs_list)) {
    row_indices <- cs_list[[i]]  # Row indices for the current credible set
    snps <- table$SNP[row_indices]  # Extract SNP values for these rows

    # Assign concatenated SNPs excluding the current SNP for each row
    for (row in row_indices) {
      table$cs_snps[row] <- paste(setdiff(snps, table$SNP[row]), collapse = ", ")
    }
  }

  # Add a label for the SNP with the smallest p-value in each credible set
  table <- table %>%
    dplyr::filter(!is.na(cs)) %>% # Exclude rows where cs is NA
    dplyr::group_by(cs) %>%
    dplyr::mutate(label = if_else(P == min(P), SNP, NA_character_)) %>%  # Identify SNP with smallest p-value
    dplyr::ungroup() %>%
    dplyr::bind_rows(table %>% dplyr::filter(is.na(cs))) %>% # Add back rows with NA cs
    dplyr::mutate(test = "Finimom")

  # Return the final table
  return(table)
}

