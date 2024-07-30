#' Generate random strings
#'
#' This function generates random strings of specified length.
#'
#' @param n Number of random strings to generate.
#' @param len Length of each random string.
#' @return A vector of randomly generated strings.
#' @export
random_string <- function(n = 1, len = 6) {
  randomString <- c(1:n)
  for (i in 1:n) {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS), len, replace = TRUE), collapse = "")
  }
  return(randomString)
}

#' Create IDs from input data
#'
#' This function creates IDs from the input data by converting the input to factors and assigning random strings as levels.
#'
#' @param x Input data.
#' @return A vector of IDs.
#' @export
create_ids <- function(x) {
  a <- as.factor(x)
  levels(a) <- random_string(length(levels(a)))
  a <- as.character(a)
  return(a)
}

#' Format data for Mendelian Randomization (MR) analysis
#'
#' This function formats input data for Mendelian Randomization (MR) analysis, ensuring required columns are present and have appropriate data types.
#'
#' @param dat Input data frame.
#' @param type Type of data (e.g., "exposure", "outcome").
#' @param snps Vector of SNP IDs to subset data.
#' @param header Logical indicating if the input data has headers.
#' @param phenotype_col Name of the phenotype column.
#' @param snp_col Name of the SNP column.
#' @param beta_col Name of the beta column.
#' @param se_col Name of the standard error column.
#' @param eaf_col Name of the effect allele frequency column.
#' @param effect_allele_col Name of the effect allele column.
#' @param other_allele_col Name of the other allele column.
#' @param pval_col Name of the p-value column.
#' @param units_col Name of the units column.
#' @param ncase_col Name of the case count column.
#' @param ncontrol_col Name of the control count column.
#' @param samplesize_col Name of the sample size column.
#' @param gene_col Name of the gene column.
#' @param id_col Name of the ID column.
#' @param min_pval Minimum p-value threshold.
#' @param z_col Name of the z-score column.
#' @param info_col Name of the INFO column.
#' @param chr_col Name of the chromosome column.
#' @param pos_col Name of the position column.
#' @param log_pval Logical indicating if the p-value should be log-transformed.
#' @return Formatted data frame for MR analysis.
#' @export
format_data <- function(dat, type = "exposure", snps = NULL, header = TRUE,
                        phenotype_col = "Phenotype", snp_col = "SNP",
                        beta_col = "beta", se_col = "se", eaf_col = "eaf",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele", pval_col = "pval",
                        units_col = "units", ncase_col = "ncase",
                        ncontrol_col = "ncontrol", samplesize_col = "samplesize",
                        gene_col = "gene", id_col = "id", min_pval = 1e-200,
                        z_col = "z", info_col = "info", chr_col = "chr",
                        pos_col = "pos", log_pval = FALSE) {
  # Extract columns based on specified column names
  all_cols <- c(phenotype_col, snp_col, beta_col, se_col, eaf_col, effect_allele_col, other_allele_col, pval_col, units_col, ncase_col, ncontrol_col, samplesize_col, gene_col, id_col, z_col, info_col, chr_col, pos_col)
  i <- names(dat) %in% all_cols
  if (sum(i) == 0) {
    stop("None of the specified columns present")
  }
  dat <- dat[, ..i]

  # Check if SNP column exists
  if (!snp_col %in% names(dat)) {
    stop("SNP column not found")
  }

  # Clean and subset data
  names(dat)[names(dat) == snp_col] <- "SNP"
  snp_col <- "SNP"
  dat$SNP <- tolower(gsub("[[:space:]]", "", dat$SNP))
  dat <- subset(dat, !is.na(SNP))
  if (!is.null(snps)) {
    dat <- subset(dat, SNP %in% snps)
  }

  # Handle phenotype column
  if (!phenotype_col %in% names(dat)) {
    message("No phenotype name specified, defaulting to '", type, "'.")
    dat[[type]] <- type
  } else {
    dat[[type]] <- dat[[phenotype_col]]
    if (phenotype_col != type) {
      dat <- dat[, -which(names(dat) == phenotype_col)]
    }
  }

  # Handle p-value transformation
  if (log_pval) {
    dat$pval <- 10^-dat[[pval_col]]
  }

  # Remove duplicated SNPs
  dat <- plyr::ddply(dat, type, function(x) {
    x <- plyr::mutate(x)
    dup <- duplicated(x$SNP)
    if (any(dup)) {
      warning("Duplicated SNPs present in exposure data for phenotype '", x[[type]][1], ". Just keeping the first instance:\n", paste(x$SNP[dup], collapse="\n"))
      x <- x[!dup,]
    }
    return(x)
  })

  # Check if columns required for MR are present
  # Additional checks and transformations...

  return(dat)
}
