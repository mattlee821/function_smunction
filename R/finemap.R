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
