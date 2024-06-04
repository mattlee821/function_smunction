#' Format Data for Mendelian Randomization Analysis
#'
#' This function formats input data for Mendelian Randomization (MR) analysis.
#' It extracts columns based on specified column names, cleans and subsets the data,
#' handles phenotype column, transforms p-values if required, removes duplicated SNPs,
#' and performs additional checks and transformations required for MR analysis.
#'
#' @param dat A data frame containing the input data.
#' @param type A character string specifying the name of the phenotype column. Default is "exposure".
#' @param snps A vector of SNP IDs to subset the data. Default is NULL.
#' @param header Logical indicating whether the data has a header row. Default is TRUE.
#' @param phenotype_col A character string specifying the name of the phenotype column in the input data. Default is "Phenotype".
#' @param snp_col A character string specifying the name of the SNP column in the input data. Default is "SNP".
#' @param beta_col A character string specifying the name of the beta coefficient column in the input data. Default is "beta".
#' @param se_col A character string specifying the name of the standard error column in the input data. Default is "se".
#' @param eaf_col A character string specifying the name of the effect allele frequency column in the input data. Default is "eaf".
#' @param effect_allele_col A character string specifying the name of the effect allele column in the input data. Default is "effect_allele".
#' @param other_allele_col A character string specifying the name of the other allele column in the input data. Default is "other_allele".
#' @param pval_col A character string specifying the name of the p-value column in the input data. Default is "pval".
#' @param units_col A character string specifying the name of the units column in the input data. Default is "units".
#' @param ncase_col A character string specifying the name of the number of cases column in the input data. Default is "ncase".
#' @param ncontrol_col A character string specifying the name of the number of controls column in the input data. Default is "ncontrol".
#' @param samplesize_col A character string specifying the name of the sample size column in the input data. Default is "samplesize".
#' @param gene_col A character string specifying the name of the gene column in the input data. Default is "gene".
#' @param id_col A character string specifying the name of the ID column in the input data. Default is "id".
#' @param min_pval A numeric value specifying the minimum p-value threshold. Default is 1e-200.
#' @param z_col A character string specifying the name of the Z-score column in the input data. Default is "z".
#' @param info_col A character string specifying the name of the INFO column in the input data. Default is "info".
#' @param chr_col A character string specifying the name of the chromosome column in the input data. Default is "chr".
#' @param pos_col A character string specifying the name of the position column in the input data. Default is "pos".
#' @param log_pval Logical indicating whether to transform p-values to -log10 scale. Default is FALSE.
#' @return A formatted data frame ready for Mendelian Randomization analysis.
#' @export
format_data <- function (dat, type = "exposure", snps = NULL, header = TRUE,
                         phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                         se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                         other_allele_col = "other_allele", pval_col = "pval", units_col = "units",
                         ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize",
                         gene_col = "gene", id_col = "id", min_pval = 1e-200, z_col = "z",
                         info_col = "info", chr_col = "chr", pos_col = "pos", log_pval = FALSE)
{
  all_cols <- c(phenotype_col, snp_col, beta_col, se_col,
                eaf_col, effect_allele_col, other_allele_col, pval_col,
                units_col, ncase_col, ncontrol_col, samplesize_col,
                gene_col, id_col, z_col, info_col, chr_col, pos_col)
  i <- names(dat) %in% all_cols
  if (sum(i) == 0) {
    stop("None of the specified columns present")
  }
  dat <- dat[, ..i]
  if (!snp_col %in% names(dat)) {
    stop("SNP column not found")
  }
  names(dat)[names(dat) == snp_col] <- "SNP"
  snp_col <- "SNP"
  dat$SNP <- tolower(dat$SNP)
  dat$SNP <- gsub("[[:space:]]", "", dat$SNP)
  dat <- subset(dat, !is.na(SNP))
  if (!is.null(snps)) {
    dat <- subset(dat, SNP %in% snps)
  }
  if (!phenotype_col %in% names(dat)) {
    message("No phenotype name specified, defaulting to '", type, "'.")
    dat[[type]] <- type
  } else {
    dat[[type]] <- dat[[phenotype_col]]
  }
  if (log_pval) {
    dat$pval <- 10^-dat[[pval_col]]
  }

  # dat <- plyr::ddply(dat, type, function(x) {
  #   x <- plyr::mutate(x)
  #   dup <- duplicated(x$SNP)
  #   if (any(dup)) {
  #     warning("Duplicated SNPs present in exposure data for phenotype '",
  #             x[[type]][1], ". Just keeping the first instance:\n",
  #             paste(x$SNP[dup], collapse = "\n"))
  #     x <- x[!dup, ]
  #   }
  #   return(x)
  # })

  dat <- dat %>%
    group_by(!!sym(type)) %>%
    distinct(!!sym(snp_col), .keep_all = TRUE) %>%
    ungroup()

  mr_cols_required <- c(snp_col, beta_col, se_col, effect_allele_col)
  mr_cols_desired <- c(other_allele_col, eaf_col)
  if (!all(mr_cols_required %in% names(dat))) {
    warning("The following columns are not present and are required for MR analysis\n",
            paste(mr_cols_required[!mr_cols_required %in% names(dat)]),
            collapse = "\n")
    dat$mr_keep.outcome <- FALSE
  }
  else {
    dat$mr_keep.outcome <- TRUE
  }
  if (!all(mr_cols_desired %in% names(dat))) {
    warning("The following columns are not present but are helpful for harmonisation\n",
            paste(mr_cols_desired[!mr_cols_desired %in% names(dat)]),
            collapse = "\n")
  }
  i <- which(names(dat) == beta_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "beta.outcome"
    if (!is.numeric(dat$beta.outcome)) {
      warning("beta column is not numeric. Coercing...")
      dat$beta.outcome <- as.numeric(dat$beta.outcome)
    }
    index <- !is.finite(dat$beta.outcome)
    index[is.na(index)] <- TRUE
    dat$beta.outcome[index] <- NA
  }
  i <- which(names(dat) == se_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "se.outcome"
    if (!is.numeric(dat$se.outcome)) {
      warning("se column is not numeric. Coercing...")
      dat$se.outcome <- as.numeric(dat$se.outcome)
    }
    index <- !is.finite(dat$se.outcome) | dat$se.outcome <=
      0
    index[is.na(index)] <- TRUE
    dat$se.outcome[index] <- NA
  }
  i <- which(names(dat) == eaf_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "eaf.outcome"
    if (!is.numeric(dat$eaf.outcome)) {
      warning("eaf column is not numeric. Coercing...")
      dat$eaf.outcome <- as.numeric(dat$eaf.outcome)
    }
    index <- !is.finite(dat$eaf.outcome) | dat$eaf.outcome <=
      0 | dat$eaf.outcome >= 1
    index[is.na(index)] <- TRUE
    dat$eaf.outcome[index] <- NA
  }
  i <- which(names(dat) == effect_allele_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "effect_allele.outcome"
    if (is.logical(dat$effect_allele.outcome)) {
      dat$effect_allele.outcome <- substr(as.character(dat$effect_allele.outcome),
                                          1, 1)
    }
    if (!is.character(dat$effect_allele.outcome)) {
      warning("effect_allele column is not character data. Coercing...")
      dat$effect_allele.outcome <- as.character(dat$effect_allele.outcome)
    }
    dat$effect_allele.outcome <- toupper(dat$effect_allele.outcome)
    index <- !(grepl("^[ACTG]+$", dat$effect_allele.outcome) |
                 dat$effect_allele.outcome %in% c("D", "I"))
    index[is.na(index)] <- TRUE
    if (any(index)) {
      warning("effect_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded.")
      dat$effect_allele.outcome[index] <- NA
      dat$mr_keep.outcome[index] <- FALSE
    }
  }
  i <- which(names(dat) == other_allele_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "other_allele.outcome"
    if (is.logical(dat$other_allele.outcome)) {
      dat$other_allele.outcome <- substr(as.character(dat$other_allele.outcome),
                                         1, 1)
    }
    if (!is.character(dat$other_allele.outcome)) {
      warning("other_allele column is not character data. Coercing...")
      dat$other_allele.outcome <- as.character(dat$other_allele.outcome)
    }
    dat$other_allele.outcome <- toupper(dat$other_allele.outcome)
    index <- !(grepl("^[ACTG]+$", dat$other_allele.outcome) |
                 dat$other_allele.outcome %in% c("D", "I"))
    index[is.na(index)] <- TRUE
    if (any(index)) {
      warning("other_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded")
      dat$other_allele.outcome[index] <- NA
      dat$mr_keep.outcome[index] <- FALSE
    }
  }
  i <- which(names(dat) == pval_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "pval.outcome"
    if (!is.numeric(dat$pval.outcome)) {
      warning("pval column is not numeric. Coercing...")
      dat$pval.outcome <- as.numeric(dat$pval.outcome)
    }
    index <- !is.finite(dat$pval.outcome) | dat$pval.outcome <
      0 | dat$pval.outcome > 1
    index[is.na(index)] <- TRUE
    dat$pval.outcome[index] <- NA
    index <- dat$pval.outcome < min_pval
    index[is.na(index)] <- FALSE
    dat$pval.outcome[index] <- min_pval
    dat$pval_origin.outcome <- "reported"
    if (any(is.na(dat$pval.outcome))) {
      if ("beta.outcome" %in% names(dat) & "se.outcome" %in%
          names(dat)) {
        index <- is.na(dat$pval.outcome)
        dat$pval.outcome[index] <- stats::pnorm(abs(dat$beta.outcome[index])/dat$se.outcome[index],
                                                lower.tail = FALSE)
        dat$pval_origin.outcome[index] <- "inferred"
      }
    }
  }
  if ("beta.outcome" %in% names(dat) & "se.outcome" %in% names(dat) &
      !"pval.outcome" %in% names(dat)) {
    message("Inferring p-values")
    dat$pval.outcome <- stats::pnorm(abs(dat$beta.outcome)/dat$se.outcome,
                                     lower.tail = FALSE) * 2
    dat$pval_origin.outcome <- "inferred"
  }
  if (ncase_col %in% names(dat)) {
    names(dat)[which(names(dat) == ncase_col)[1]] <- "ncase.outcome"
    if (!is.numeric(dat$ncase.outcome)) {
      warning(ncase_col, " column is not numeric")
      dat$ncase.outcome <- as.numeric(dat$ncase.outcome)
    }
  }
  if (ncontrol_col %in% names(dat)) {
    names(dat)[which(names(dat) == ncontrol_col)[1]] <- "ncontrol.outcome"
    if (!is.numeric(dat$ncontrol.outcome)) {
      warning(ncontrol_col, " column is not numeric")
      dat$ncontrol.outcome <- as.numeric(dat$ncontrol.outcome)
    }
  }
  if (samplesize_col %in% names(dat)) {
    names(dat)[which(names(dat) == samplesize_col)[1]] <- "samplesize.outcome"
    if (!is.numeric(dat$samplesize.outcome)) {
      warning(samplesize_col, " column is not numeric")
      dat$samplesize.outcome <- as.numeric(dat$samplesize.outcome)
    }
    if ("ncontrol.outcome" %in% names(dat) & "ncase.outcome" %in%
        names(dat)) {
      index <- is.na(dat$samplesize.outcome) & !is.na(dat$ncase.outcome) &
        !is.na(dat$ncontrol.outcome)
      if (any(index)) {
        message("Generating sample size from ncase and ncontrol")
        dat$samplesize.outcome[index] <- dat$ncase.outcome[index] +
          dat$ncontrol.outcome[index]
      }
    }
  }
  else if ("ncontrol.outcome" %in% names(dat) & "ncase.outcome" %in%
           names(dat)) {
    message("Generating sample size from ncase and ncontrol")
    dat$samplesize.outcome <- dat$ncase.outcome + dat$ncontrol.outcome
  }
  if (gene_col %in% names(dat)) {
    names(dat)[which(names(dat) == gene_col)[1]] <- "gene.outcome"
  }
  if (info_col %in% names(dat)) {
    names(dat)[which(names(dat) == info_col)[1]] <- "info.outcome"
  }
  if (z_col %in% names(dat)) {
    names(dat)[which(names(dat) == z_col)[1]] <- "z.outcome"
  }
  if (chr_col %in% names(dat)) {
    names(dat)[which(names(dat) == chr_col)[1]] <- "chr.outcome"
  }
  if (pos_col %in% names(dat)) {
    names(dat)[which(names(dat) == pos_col)[1]] <- "pos.outcome"
  }
  if (units_col %in% names(dat)) {
    names(dat)[which(names(dat) == units_col)[1]] <- "units.outcome"
    dat$units.outcome_dat <- as.character(dat$units.outcome)
    temp <- check_units(dat, type, "units.outcome")
    if (any(temp$ph)) {
      dat[[type]] <- paste0(dat[[type]], " (", dat$units.outcome,
                            ")")
    }
  }
  if (id_col %in% names(dat)) {
    names(dat)[which(names(dat) == id_col)[1]] <- "id.outcome"
    dat$id.outcome <- as.character(dat$id.outcome)
  }
  else {
    dat$id.outcome <- create_ids(dat[[type]])
  }
  if (any(dat$mr_keep.outcome)) {
    mrcols <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome")
    mrcols_present <- mrcols[mrcols %in% names(dat)]
    dat$mr_keep.outcome <- dat$mr_keep.outcome & apply(dat[,
                                                           mrcols_present], 1, function(x) !any(is.na(x)))
    if (any(!dat$mr_keep.outcome)) {
      warning("The following SNP(s) are missing required information for the MR tests and will be excluded\n",
              paste(subset(dat, !mr_keep.outcome)$SNP, collapse = "\n"))
    }
  }
  if (all(!dat$mr_keep.outcome)) {
    warning("None of the provided SNPs can be used for MR analysis, they are missing required information.")
  }
  for (col in c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome",
                "other_allele.outcome", "eaf.outcome")) {
    if (!col %in% names(dat)) {
      dat[[col]] <- NA
    }
  }
  names(dat) <- gsub("outcome", type, names(dat))
  rownames(dat) <- NULL
  return(dat)
}

#' Format Data Function
#'
#' This function formats a dataset for analysis, particularly for Mendelian Randomization (MR) studies.
#'
#' @param dat The dataset to be formatted.
#' @param type Character string specifying the type of data, defaults to "exposure".
#' @param snps Vector of SNP IDs to subset the data.
#' @param header Logical indicating if the dataset has column names, defaults to TRUE.
#' @param phenotype_col Character string specifying the column name for phenotype data, defaults to "Phenotype".
#' @param snp_col Character string specifying the column name for SNP IDs, defaults to "SNP".
#' @param beta_col Character string specifying the column name for beta coefficients, defaults to "beta".
#' @param se_col Character string specifying the column name for standard errors, defaults to "se".
#' @param eaf_col Character string specifying the column name for effect allele frequency, defaults to "eaf".
#' @param effect_allele_col Character string specifying the column name for effect allele, defaults to "effect_allele".
#' @param other_allele_col Character string specifying the column name for other allele, defaults to "other_allele".
#' @param pval_col Character string specifying the column name for p-values, defaults to "pval".
#' @param units_col Character string specifying the column name for units, defaults to "units".
#' @param ncase_col Character string specifying the column name for number of cases, defaults to "ncase".
#' @param ncontrol_col Character string specifying the column name for number of controls, defaults to "ncontrol".
#' @param samplesize_col Character string specifying the column name for sample size, defaults to "samplesize".
#' @param gene_col Character string specifying the column name for gene, defaults to "gene".
#' @param id_col Character string specifying the column name for IDs, defaults to "id".
#' @param min_pval Minimum p-value to replace non-finite or extreme p-values, defaults to 1e-200.
#' @param z_col Character string specifying the column name for z-scores, defaults to "z".
#' @param info_col Character string specifying the column name for info scores, defaults to "info".
#' @param chr_col Character string specifying the column name for chromosome, defaults to "chr".
#' @param pos_col Character string specifying the column name for position, defaults to "pos".
#' @param log_pval Logical indicating if p-values should be log-transformed, defaults to FALSE.
#' @return A formatted dataset ready for analysis.
#' @export
format_data2 <- function (dat, type = "exposure", snps = NULL, header = TRUE,
                          phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta",
                          se_col = "se", eaf_col = "eaf", effect_allele_col = "effect_allele",
                          other_allele_col = "other_allele", pval_col = "pval", units_col = "units",
                          ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize",
                          gene_col = "gene", id_col = "id", min_pval = 1e-200, z_col = "z",
                          info_col = "info", chr_col = "chr", pos_col = "pos", log_pval = FALSE)
{
  all_cols <- c(phenotype_col, snp_col, beta_col, se_col,
                eaf_col, effect_allele_col, other_allele_col, pval_col,
                units_col, ncase_col, ncontrol_col, samplesize_col,
                gene_col, id_col, z_col, info_col, chr_col, pos_col)
  i <- names(dat) %in% all_cols
  if (sum(i) == 0) {
    stop("None of the specified columns present")
  }
  dat <- dat[, i]
  if (!snp_col %in% names(dat)) {
    stop("SNP column not found")
  }
  names(dat)[names(dat) == snp_col] <- "SNP"
  snp_col <- "SNP"
  dat$SNP <- tolower(dat$SNP)
  dat$SNP <- gsub("[[:space:]]", "", dat$SNP)
  dat <- subset(dat, !is.na(SNP))
  if (!is.null(snps)) {
    dat <- subset(dat, SNP %in% snps)
  }
  if (!phenotype_col %in% names(dat)) {
    message("No phenotype name specified, defaulting to '", type, "'.")
    dat[[type]] <- type
  } else {
    dat[[type]] <- dat[[phenotype_col]]
  }
  if (log_pval) {
    dat$pval <- 10^-dat[[pval_col]]
  }

  # dat <- plyr::ddply(dat, type, function(x) {
  #   x <- plyr::mutate(x)
  #   dup <- duplicated(x$SNP)
  #   if (any(dup)) {
  #     warning("Duplicated SNPs present in exposure data for phenotype '",
  #             x[[type]][1], ". Just keeping the first instance:\n",
  #             paste(x$SNP[dup], collapse = "\n"))
  #     x <- x[!dup, ]
  #   }
  #   return(x)
  # })

  dat <- dat %>%
    group_by(!!sym(type)) %>%
    distinct(!!sym(snp_col), .keep_all = TRUE) %>%
    ungroup()

  mr_cols_required <- c(snp_col, beta_col, se_col, effect_allele_col)
  mr_cols_desired <- c(other_allele_col, eaf_col)
  if (!all(mr_cols_required %in% names(dat))) {
    warning("The following columns are not present and are required for MR analysis\n",
            paste(mr_cols_required[!mr_cols_required %in% names(dat)]),
            collapse = "\n")
    dat$mr_keep.outcome <- FALSE
  }
  else {
    dat$mr_keep.outcome <- TRUE
  }
  if (!all(mr_cols_desired %in% names(dat))) {
    warning("The following columns are not present but are helpful for harmonisation\n",
            paste(mr_cols_desired[!mr_cols_desired %in% names(dat)]),
            collapse = "\n")
  }
  i <- which(names(dat) == beta_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "beta.outcome"
    if (!is.numeric(dat$beta.outcome)) {
      warning("beta column is not numeric. Coercing...")
      dat$beta.outcome <- as.numeric(dat$beta.outcome)
    }
    index <- !is.finite(dat$beta.outcome)
    index[is.na(index)] <- TRUE
    dat$beta.outcome[index] <- NA
  }
  i <- which(names(dat) == se_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "se.outcome"
    if (!is.numeric(dat$se.outcome)) {
      warning("se column is not numeric. Coercing...")
      dat$se.outcome <- as.numeric(dat$se.outcome)
    }
    index <- !is.finite(dat$se.outcome) | dat$se.outcome <=
      0
    index[is.na(index)] <- TRUE
    dat$se.outcome[index] <- NA
  }
  i <- which(names(dat) == eaf_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "eaf.outcome"
    if (!is.numeric(dat$eaf.outcome)) {
      warning("eaf column is not numeric. Coercing...")
      dat$eaf.outcome <- as.numeric(dat$eaf.outcome)
    }
    index <- !is.finite(dat$eaf.outcome) | dat$eaf.outcome <=
      0 | dat$eaf.outcome >= 1
    index[is.na(index)] <- TRUE
    dat$eaf.outcome[index] <- NA
  }
  i <- which(names(dat) == effect_allele_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "effect_allele.outcome"
    if (is.logical(dat$effect_allele.outcome)) {
      dat$effect_allele.outcome <- substr(as.character(dat$effect_allele.outcome),
                                          1, 1)
    }
    if (!is.character(dat$effect_allele.outcome)) {
      warning("effect_allele column is not character data. Coercing...")
      dat$effect_allele.outcome <- as.character(dat$effect_allele.outcome)
    }
    dat$effect_allele.outcome <- toupper(dat$effect_allele.outcome)
    index <- !(grepl("^[ACTG]+$", dat$effect_allele.outcome) |
                 dat$effect_allele.outcome %in% c("D", "I"))
    index[is.na(index)] <- TRUE
    if (any(index)) {
      warning("effect_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded.")
      dat$effect_allele.outcome[index] <- NA
      dat$mr_keep.outcome[index] <- FALSE
    }
  }
  i <- which(names(dat) == other_allele_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "other_allele.outcome"
    if (is.logical(dat$other_allele.outcome)) {
      dat$other_allele.outcome <- substr(as.character(dat$other_allele.outcome),
                                         1, 1)
    }
    if (!is.character(dat$other_allele.outcome)) {
      warning("other_allele column is not character data. Coercing...")
      dat$other_allele.outcome <- as.character(dat$other_allele.outcome)
    }
    dat$other_allele.outcome <- toupper(dat$other_allele.outcome)
    index <- !(grepl("^[ACTG]+$", dat$other_allele.outcome) |
                 dat$other_allele.outcome %in% c("D", "I"))
    index[is.na(index)] <- TRUE
    if (any(index)) {
      warning("other_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded")
      dat$other_allele.outcome[index] <- NA
      dat$mr_keep.outcome[index] <- FALSE
    }
  }
  i <- which(names(dat) == pval_col)[1]
  if (!is.na(i)) {
    names(dat)[i] <- "pval.outcome"
    if (!is.numeric(dat$pval.outcome)) {
      warning("pval column is not numeric. Coercing...")
      dat$pval.outcome <- as.numeric(dat$pval.outcome)
    }
    index <- !is.finite(dat$pval.outcome) | dat$pval.outcome <
      0 | dat$pval.outcome > 1
    index[is.na(index)] <- TRUE
    dat$pval.outcome[index] <- NA
    index <- dat$pval.outcome < min_pval
    index[is.na(index)] <- FALSE
    dat$pval.outcome[index] <- min_pval
    dat$pval_origin.outcome <- "reported"
    if (any(is.na(dat$pval.outcome))) {
      if ("beta.outcome" %in% names(dat) & "se.outcome" %in%
          names(dat)) {
        index <- is.na(dat$pval.outcome)
        dat$pval.outcome[index] <- stats::pnorm(abs(dat$beta.outcome[index])/dat$se.outcome[index],
                                                lower.tail = FALSE)
        dat$pval_origin.outcome[index] <- "inferred"
      }
    }
  }
  if ("beta.outcome" %in% names(dat) & "se.outcome" %in% names(dat) &
      !"pval.outcome" %in% names(dat)) {
    message("Inferring p-values")
    dat$pval.outcome <- stats::pnorm(abs(dat$beta.outcome)/dat$se.outcome,
                                     lower.tail = FALSE) * 2
    dat$pval_origin.outcome <- "inferred"
  }
  if (ncase_col %in% names(dat)) {
    names(dat)[which(names(dat) == ncase_col)[1]] <- "ncase.outcome"
    if (!is.numeric(dat$ncase.outcome)) {
      warning(ncase_col, " column is not numeric")
      dat$ncase.outcome <- as.numeric(dat$ncase.outcome)
    }
  }
  if (ncontrol_col %in% names(dat)) {
    names(dat)[which(names(dat) == ncontrol_col)[1]] <- "ncontrol.outcome"
    if (!is.numeric(dat$ncontrol.outcome)) {
      warning(ncontrol_col, " column is not numeric")
      dat$ncontrol.outcome <- as.numeric(dat$ncontrol.outcome)
    }
  }
  if (samplesize_col %in% names(dat)) {
    names(dat)[which(names(dat) == samplesize_col)[1]] <- "samplesize.outcome"
    if (!is.numeric(dat$samplesize.outcome)) {
      warning(samplesize_col, " column is not numeric")
      dat$samplesize.outcome <- as.numeric(dat$samplesize.outcome)
    }
    if ("ncontrol.outcome" %in% names(dat) & "ncase.outcome" %in%
        names(dat)) {
      index <- is.na(dat$samplesize.outcome) & !is.na(dat$ncase.outcome) &
        !is.na(dat$ncontrol.outcome)
      if (any(index)) {
        message("Generating sample size from ncase and ncontrol")
        dat$samplesize.outcome[index] <- dat$ncase.outcome[index] +
          dat$ncontrol.outcome[index]
      }
    }
  }
  else if ("ncontrol.outcome" %in% names(dat) & "ncase.outcome" %in%
           names(dat)) {
    message("Generating sample size from ncase and ncontrol")
    dat$samplesize.outcome <- dat$ncase.outcome + dat$ncontrol.outcome
  }
  if (gene_col %in% names(dat)) {
    names(dat)[which(names(dat) == gene_col)[1]] <- "gene.outcome"
  }
  if (info_col %in% names(dat)) {
    names(dat)[which(names(dat) == info_col)[1]] <- "info.outcome"
  }
  if (z_col %in% names(dat)) {
    names(dat)[which(names(dat) == z_col)[1]] <- "z.outcome"
  }
  if (chr_col %in% names(dat)) {
    names(dat)[which(names(dat) == chr_col)[1]] <- "chr.outcome"
  }
  if (pos_col %in% names(dat)) {
    names(dat)[which(names(dat) == pos_col)[1]] <- "pos.outcome"
  }
  if (units_col %in% names(dat)) {
    names(dat)[which(names(dat) == units_col)[1]] <- "units.outcome"
    dat$units.outcome_dat <- as.character(dat$units.outcome)
    temp <- check_units(dat, type, "units.outcome")
    if (any(temp$ph)) {
      dat[[type]] <- paste0(dat[[type]], " (", dat$units.outcome,
                            ")")
    }
  }
  if (id_col %in% names(dat)) {
    names(dat)[which(names(dat) == id_col)[1]] <- "id.outcome"
    dat$id.outcome <- as.character(dat$id.outcome)
  }
  else {
    dat$id.outcome <- create_ids(dat[[type]])
  }
  if (any(dat$mr_keep.outcome)) {
    mrcols <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome")
    mrcols_present <- mrcols[mrcols %in% names(dat)]
    dat$mr_keep.outcome <- dat$mr_keep.outcome & apply(dat[,
                                                           mrcols_present], 1, function(x) !any(is.na(x)))
    if (any(!dat$mr_keep.outcome)) {
      warning("The following SNP(s) are missing required information for the MR tests and will be excluded\n",
              paste(subset(dat, !mr_keep.outcome)$SNP, collapse = "\n"))
    }
  }
  if (all(!dat$mr_keep.outcome)) {
    warning("None of the provided SNPs can be used for MR analysis, they are missing required information.")
  }
  for (col in c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome",
                "other_allele.outcome", "eaf.outcome")) {
    if (!col %in% names(dat)) {
      dat[[col]] <- NA
    }
  }
  names(dat) <- gsub("outcome", type, names(dat))
  rownames(dat) <- NULL
  return(dat)
}

#' Remove duplicate rows based on SNP column
#'
#' This function removes duplicate rows based on the "SNP" column, prioritizing rows with "cis" in the "cis_trans" column.
#' For a given SNP, it keeps only one row, prioritizing the row with "cis" in "cis_trans" and removing others.
#'
#' Example usage:
#' \code{#' for (i in seq_along(list_data_cis_trans)) {
#'   df <- list_data_cis_trans[[i]]
#'   df <- remove_duplicate_SNP(df)
#'   df <- remove_nearby_positions(df)
#'   list_data_cis_trans[[i]] <- df
#' }
#' }
#'
#' @param df A data frame containing SNP and cis_trans columns.
#' @return A modified data frame with duplicate rows removed based on SNP, prioritizing "cis" rows.
#' @export
remove_duplicate_SNP <- function(df) {
  cis_rows <- which(df$cis_trans == "cis")
  for (row_index in cis_rows) {
    SNP_value <- df$SNP[row_index]
    duplicate_rows <- which(df$SNP == SNP_value & df$cis_trans != "cis")
    if (length(duplicate_rows) > 0) {
      df <- df[-duplicate_rows, ]
    }
  }
  return(df)
}

#' Remove rows with nearby positions (excluding "cis" and "cis-trans_cis")
#'
#' This function removes rows with positions within ±1mb of a "cis" or "cis-trans_cis" row,
#' excluding the "cis" or "cis-trans_cis" row itself. It iterates through rows with "cis" or "cis-trans_cis" in "cis_trans",
#' checks for other rows with the same chromosome but positions within the threshold.
#' If such rows exist and their "cis_trans" is not "cis" or "cis-trans_cis", it removes those rows.
#'
#' Example usage:
#' \code{
#' for (i in seq_along(list_data_cis_trans)) {
#'   df <- list_data_cis_trans[[i]]
#'   df <- remove_duplicate_SNP(df)
#'   df <- remove_nearby_positions(df)
#'   list_data_cis_trans[[i]] <- df
#' }
#' }
#'
#' @param df A data frame containing chr.exposure, pos.exposure, and cis_trans columns.
#' @return A modified data frame with rows removed based on nearby positions (excluding "cis" and "cis-trans_cis").
#' @export
remove_nearby_positions <- function(df) {
  cis_rows <- which(df$cis_trans %in% c("cis", "cis-trans_cis"))
  for (row_index in cis_rows) {
    chr_exposure_value <- df$chr.exposure[row_index]
    pos_exposure_value <- df$pos.exposure[row_index]
    duplicate_chr_exposure <- which(df$chr.exposure == chr_exposure_value & !df$cis_trans %in% c("cis", "cis-trans_cis"))
    if (length(duplicate_chr_exposure) > 0) {
      rows_to_remove <- c()
      for (row_index2 in duplicate_chr_exposure) {
        # Check for missing values before accessing the cis_trans column
        if (!is.na(df$cis_trans[row_index2]) && !df$cis_trans[row_index2] %in% c("cis", "cis-trans_cis")) {
          pos_diff <- df$pos.exposure[row_index2] - pos_exposure_value
          if (abs(pos_diff) <= 1000000) {
            rows_to_remove <- c(rows_to_remove, row_index2)
          }
        }
      }
      # Remove rows outside the inner loop to prevent modifying the data frame while iterating
      df <- df[-rows_to_remove, ]
    }
  }
  return(df)
}

#' @title proxy_search: search for proxies when using read_outcome_data()
#' @description
#' This function gets proxy SNPs for MR analyses when using local outcome data.
#' It takes formatted exposure and outcome data from `TwoSampleMR`, searches for
#' missing SNPs in the outcome, identifies which of these missing-SNPs are
#' available in the provided reference panel, extracts all proxy-SNPs for the
#' missing-SNPs from the reference panel, returns a data frame with the top
#' proxy-SNP for each missing-SNP. You MUST have a local reference panel. Only
#' works with rsID.
#' @return a data frame
#' @param data_exposure exposure data frame
#' @param data_outcome outcome data frame
#' @param data_outcome_path file path used for `read_outcome_data()`
#' @param data_reference reference data; bim file
#' @param data_reference_path file path for your downloaded reference panel
#' @param tag_r2 r2 for proxy SNP; from `get_ld_proxies()`; default = 0.8
#' @param tag_kb window to look for proxy SNPs; from `get_ld_proxies()`; default = 5000
#' @param tag_nsnp from `get_ld_proxies()`; default = 5000
#' @param outcome_sep separator for your outcome GWAS
#' @param outcome_phenotype phenotype column name of your GWAS
#' @param outcome_SNP SNP column name of your GWAS
#' @param outcome_BETA BETA column name of your GWAS
#' @param outcome_SE SE column name of your GWAS
#' @param outcome_P P column name of your GWAS
#' @param outcome_EA EA column name of your GWAS
#' @param outcome_OA OA column name of your GWAS
#' @param outcome_EAF EAF column name of your GWAS
#' @param outcome_N N column name of your GWAS
#' @param outcome_ID ID column name of your GWAS
#' @param outcome_CHR CHR column name of your GWAS
#' @param outcome_POS POS column name of your GWAS
#' @export
proxy_search <- function(data_exposure, data_outcome, data_outcome_path, data_reference, data_reference_path,
                         tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                         outcome_sep, outcome_phenotype, outcome_SNP, outcome_BETA, outcome_SE, outcome_P,
                         outcome_EA, outcome_OA, outcome_EAF, outcome_N, outcome_ID, outcome_CHR, outcome_POS) {

  # Parameter Validation
  if (!file.exists(data_outcome_path) || !file.exists(paste0(data_reference))) {
    stop("Invalid file paths provided.")
  }

  if (!is.data.frame(data_exposure) || !is.data.frame(data_outcome)) {
    stop("data_exposure and data_outcome must be data frames.")
  }

  # exposure snps missing from outcome ====
  snps_missing <- setdiff(unique(as.factor(data_exposure$SNP)), unique(as.factor(data_outcome$SNP)))
  if (length(snps_missing) == 0) {
    message("No missing SNPs. Returning the original dataframe.")
    return(data_outcome)
  }

  # look-up missing SNPs in reference panel ====
  message(paste0("# 1. looking up ", length(unique(as.factor(snps_missing))), " missing-SNP(s) in the reference panel"))
  reference_header <- readLines(data_reference, n = 1) ## read the first row to identify the column containing "rs*"
  column_index <- grep("rs", strsplit(reference_header, "\t")[[1]]) ## get column index
  ## bash command to filter rows based on the identified column
  cmd <- sprintf("awk '$%d ~ /^(%s)$/' %s", column_index, paste(snps_missing, collapse = "|"), paste0(data_reference))
  ## use system to run the bash command and then read the result with fread
  reference <- data.table::fread(cmd = cmd, quote = "")
  snps_reference <- intersect(unique(as.factor(snps_missing)), unique(as.factor(reference$V2)))
  message(paste0("## ", length(unique(as.factor(snps_reference))), " of ", length(unique(as.factor(snps_missing))), " missing-SNP(s) available in the reference panel"))

  # find proxies for available SNPs ====
  message(paste0("# 2. extracting proxy-SNP(s) for the ", length(unique(as.factor(snps_reference))), " missing-SNP(s) from the reference panel"))
  gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
  proxies <- functions::get_ld_proxies(rsid = snps_missing,
                                       bfile = data_reference_path,
                                       searchspace = NULL,
                                       tag_kb = tag_kb,
                                       tag_nsnp = tag_nsnp,
                                       tag_r2 = tag_r2,
                                       threads = 1,
                                       out = tempfile())
  ## format proxy data: change column order and names, add proxy.outcome = TRUE
  proxies <- proxies %>%
    dplyr::select(target_snp.outcome = SNP_A,
                  proxy_snp.outcome = SNP_B,
                  target_a1.outcome = A1,
                  target_a2.outcome = A2,
                  proxy_a1.outcome = B1,
                  proxy_a2.outcome = B2,
                  R) %>%
    dplyr::mutate(proxy.outcome = TRUE,
                  SNP = proxy_snp.outcome) %>%
    dplyr::select(proxy.outcome, everything())
  message(paste0("## proxy-SNP(s) for ", length(unique(proxies$target_snp.outcome)), " missing-SNP(s) found; ", "proxy-SNP(s) for ", length(unique(as.factor(snps_reference))) - length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) not available (e.g., no proxy-SNP or r2 < provided)"))

  # extract proxies from outcome ====
  message(paste0("# 3. extracting proxy-SNP(s) from outcome"))
  proxy_snps <- proxies %>% # select unique proxy SNPs to extract
    dplyr::distinct(proxy_snp.outcome) %>%
    dplyr::pull(proxy_snp.outcome)
  data_outcome_proxies <- TwoSampleMR::read_outcome_data(filename = data_outcome_path,
                                                         snps = proxy_snps,
                                                         sep = outcome_sep,
                                                         phenotype_col = outcome_phenotype,
                                                         snp_col = outcome_SNP,
                                                         beta_col = outcome_BETA,
                                                         se_col = outcome_SE,
                                                         eaf_col = outcome_EAF,
                                                         effect_allele_col = outcome_EA,
                                                         other_allele_col = outcome_OA,
                                                         pval_col = outcome_P,
                                                         samplesize_col = outcome_N,
                                                         id_col = outcome_ID,
                                                         chr_col = outcome_CHR,
                                                         pos_col = outcome_POS)
  data_outcome_proxies <- dplyr::left_join(data_outcome_proxies, proxies, by = c("SNP" = "SNP"))
  message(paste0("## proxy-SNP(s) for ", length(unique(as.factor(data_outcome_proxies$target_snp.outcome))), " of ", length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) extracted"))

  # select proxy-SNP(s) with the highest R2 ====
  data_outcome_proxies <- data_outcome_proxies %>%
    dplyr::group_by(target_snp.outcome) %>%
    dplyr::filter(R == max(R)) %>%
    dplyr::slice(1) %>%
    dplyr::select(-R)

  ## Bind rows of data_outcome with data_outcome_proxies
  data_outcome <- dplyr::bind_rows(data_outcome, data_outcome_proxies)

  return(data_outcome)
}

#' finds proxy SNPs from an internal reference
#' split the SNP search into 10k blocks
#'
#' @param data_exposure Data frame containing exposure SNP data.
#' @param data_outcome Data frame containing outcome SNP data.
#' @param data_outcome_path Path to the outcome data file.
#' @param data_reference Data frame containing reference SNP data.
#' @param data_reference_path Path to the reference data file.
#' @param tag_r2 Threshold for LD (default: 0.8).
#' @param tag_kb Threshold for distance in kb (default: 5000).
#' @param tag_nsnp Threshold for number of SNPs (default: 5000).
#' @param outcome_sep Delimiter used in the outcome data file.
#' @param outcome_phenotype Column name for the outcome phenotype.
#' @param outcome_SNP Column name for the SNP in the outcome data.
#' @param outcome_BETA Column name for the beta value in the outcome data.
#' @param outcome_SE Column name for the standard error in the outcome data.
#' @param outcome_P Column name for the p-value in the outcome data.
#' @param outcome_EA Column name for the effect allele in the outcome data.
#' @param outcome_OA Column name for the other allele in the outcome data.
#' @param outcome_EAF Column name for the effect allele frequency in the outcome data.
#' @param outcome_N Column name for the sample size in the outcome data.
#' @param outcome_ID Column name for the ID in the outcome data.
#' @param outcome_CHR Column name for the chromosome in the outcome data.
#' @param outcome_POS Column name for the position in the outcome data.
#' @return Data frame with proxy SNPs added.
#' @export
proxy_search_split <- function(data_exposure, data_outcome, data_outcome_path, data_reference, data_reference_path,
                               tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                               outcome_sep, outcome_phenotype, outcome_SNP, outcome_BETA, outcome_SE, outcome_P,
                               outcome_EA, outcome_OA, outcome_EAF, outcome_N, outcome_ID, outcome_CHR, outcome_POS) {

  # Parameter Validation
  if (!file.exists(data_outcome_path) || !file.exists(paste0(data_reference))) {
    stop("Invalid file paths provided.")
  }

  if (!is.data.frame(data_exposure) || !is.data.frame(data_outcome)) {
    stop("data_exposure and data_outcome must be data frames.")
  }

  # Function to process a block of missing SNPs
  process_block <- function(snps_missing_block, data_outcome) {
    # look-up missing SNPs in reference panel ====
    message(paste0("# 1. looking up ", length(unique(as.factor(snps_missing_block))), " missing-SNP(s) in the reference panel"))
    reference_header <- readLines(data_reference, n = 1) ## read the first row to identify the column containing "rs*"
    column_index <- grep("rs", strsplit(reference_header, "\t")[[1]]) ## get column index
    ## bash command to filter rows based on the identified column
    cmd <- sprintf("awk '$%d ~ /^(%s)$/' %s", column_index, paste(snps_missing_block, collapse = "|"), paste0(data_reference))
    ## use system to run the bash command and then read the result with fread
    reference <- data.table::fread(cmd = cmd, quote = "")
    snps_reference <- intersect(unique(as.factor(snps_missing_block)), unique(as.factor(reference$V2)))
    message(paste0("## ", length(unique(as.factor(snps_reference))), " of ", length(unique(as.factor(snps_missing_block))), " missing-SNP(s) available in the reference panel"))

    # find proxies for available SNPs ====
    message(paste0("# 2. extracting proxy-SNP(s) for the ", length(unique(as.factor(snps_reference))), " missing-SNP(s) from the reference panel"))
    gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
    proxies <- functions::get_ld_proxies(rsid = snps_missing_block,
                                         bfile = data_reference_path,
                                         searchspace = NULL,
                                         tag_kb = tag_kb,
                                         tag_nsnp = tag_nsnp,
                                         tag_r2 = tag_r2,
                                         threads = 1,
                                         out = tempfile())
    ## format proxy data: change column order and names, add proxy.outcome = TRUE
    proxies <- proxies %>%
      dplyr::select(target_snp.outcome = SNP_A,
                    proxy_snp.outcome = SNP_B,
                    target_a1.outcome = A1,
                    target_a2.outcome = A2,
                    proxy_a1.outcome = B1,
                    proxy_a2.outcome = B2,
                    R) %>%
      dplyr::mutate(proxy.outcome = TRUE,
                    SNP = proxy_snp.outcome) %>%
      dplyr::select(proxy.outcome, everything())
    message(paste0("## proxy-SNP(s) for ", length(unique(proxies$target_snp.outcome)), " missing-SNP(s) found; ", "proxy-SNP(s) for ", length(unique(as.factor(snps_reference))) - length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) not available (e.g., no proxy-SNP or r2 < provided)"))

    # extract proxies from outcome ====
    message(paste0("# 3. extracting proxy-SNP(s) from outcome"))
    proxy_snps <- proxies %>% # select unique proxy SNPs to extract
      dplyr::distinct(proxy_snp.outcome) %>%
      dplyr::pull(proxy_snp.outcome)
    data_outcome_proxies <- TwoSampleMR::read_outcome_data(filename = data_outcome_path,
                                                           snps = proxy_snps,
                                                           sep = outcome_sep,
                                                           phenotype_col = outcome_phenotype,
                                                           snp_col = outcome_SNP,
                                                           beta_col = outcome_BETA,
                                                           se_col = outcome_SE,
                                                           eaf_col = outcome_EAF,
                                                           effect_allele_col = outcome_EA,
                                                           other_allele_col = outcome_OA,
                                                           pval_col = outcome_P,
                                                           samplesize_col = outcome_N,
                                                           id_col = outcome_ID,
                                                           chr_col = outcome_CHR,
                                                           pos_col = outcome_POS)
    data_outcome_proxies <- dplyr::left_join(data_outcome_proxies, proxies, by = c("SNP" = "SNP"))
    message(paste0("## proxy-SNP(s) for ", length(unique(as.factor(data_outcome_proxies$target_snp.outcome))), " of ", length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) extracted"))

    # select proxy-SNP(s) with the highest R2 ====
    data_outcome_proxies <- data_outcome_proxies %>%
      dplyr::group_by(target_snp.outcome) %>%
      dplyr::filter(R == max(R)) %>%
      dplyr::slice(1) %>%
      dplyr::select(-R)

    ## Bind rows of data_outcome with data_outcome_proxies
    data_outcome <- dplyr::bind_rows(data_outcome, data_outcome_proxies)

    return(data_outcome)
  }

  # exposure snps missing from outcome ====
  snps_missing <- setdiff(unique(as.factor(data_exposure$SNP)), unique(as.factor(data_outcome$SNP)))
  if (length(snps_missing) == 0) {
    message("No missing SNPs. Returning the original dataframe.")
    return(data_outcome)
  }

  # Check if splitting is necessary
  if (length(snps_missing) > 10000) {
    message("Splitting into blocks...")
    num_blocks <- ceiling(length(snps_missing) / 10000)
    for (i in 1:num_blocks) {
      start_index <- (i - 1) * 10000 + 1
      end_index <- min(i * 10000, length(snps_missing))
      snps_missing_block <- snps_missing[start_index:end_index]
      data_outcome <- process_block(snps_missing_block, data_outcome)
    }
  } else {
    # Process all missing SNPs at once
    data_outcome <- process_block(snps_missing, data_outcome)
  }

  return(data_outcome)
}


#' finds proxy SNPs from an internal reference (DECODE format)
#'
#' @param data_exposure Data frame containing exposure SNP data.
#' @param data_outcome Data frame containing outcome SNP data.
#' @param data_outcome_path Path to the outcome data file.
#' @param data_reference Data frame containing reference SNP data.
#' @param data_reference_path Path to the reference data file.
#' @param tag_r2 Threshold for LD (default: 0.8).
#' @param tag_kb Threshold for distance in kb (default: 5000).
#' @param tag_nsnp Threshold for number of SNPs (default: 5000).
#' @param outcome_sep Delimiter used in the outcome data file.
#' @return Data frame with proxy SNPs added.
#' @export
proxy_search_DECODE <- function (data_exposure, data_outcome, data_outcome_path, data_reference,
                                 data_reference_path, tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                                 outcome_sep)
{
  if (!file.exists(data_outcome_path) || !file.exists(paste0(data_reference))) {
    stop("Invalid file paths provided.")
  }
  if (!is.data.frame(data_exposure) || !is.data.frame(data_outcome)) {
    stop("data_exposure and data_outcome must be data frames.")
  }
  snps_missing <- setdiff(unique(as.factor(data_exposure$SNP)),
                          unique(as.factor(data_outcome$SNP)))
  if (length(snps_missing) == 0) {
    message("No missing SNPs. Returning the original dataframe.")
    return(data_outcome)
  }
  message(paste0("# 1. looking up ", length(unique(as.factor(snps_missing))),
                 " missing-SNP(s) in the reference panel"))
  reference_header <- readLines(data_reference, n = 1)
  column_index <- grep("rs", strsplit(reference_header, "\t")[[1]])
  cmd <- sprintf("awk '$%d ~ /^(%s)$/' %s", column_index,
                 paste(snps_missing, collapse = "|"), paste0(data_reference))
  reference <- data.table::fread(cmd = cmd, quote = "")
  snps_reference <- intersect(unique(as.factor(snps_missing)),
                              unique(as.factor(reference$V2)))
  message(paste0("## ", length(unique(as.factor(snps_reference))),
                 " of ", length(unique(as.factor(snps_missing))), " missing-SNP(s) available in the reference panel"))
  message(paste0("# 2. extracting proxy-SNP(s) for the ",
                 length(unique(as.factor(snps_reference))), " missing-SNP(s) from the reference panel"))
  gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
  proxies <- functions::get_ld_proxies(rsid = snps_missing,
                                       bfile = data_reference_path, searchspace = NULL, tag_kb = tag_kb,
                                       tag_nsnp = tag_nsnp, tag_r2 = tag_r2, threads = 1, out = tempfile())
  proxies <- proxies %>% dplyr::select(target_snp.outcome = SNP_A,
                                       proxy_snp.outcome = SNP_B, target_a1.outcome = A1, target_a2.outcome = A2,
                                       proxy_a1.outcome = B1, proxy_a2.outcome = B2, R) %>%
    dplyr::mutate(proxy.outcome = TRUE, SNP = proxy_snp.outcome) %>%
    dplyr::select(proxy.outcome, everything())
  message(paste0("## proxy-SNP(s) for ", length(unique(proxies$target_snp.outcome)),
                 " missing-SNP(s) found; ", "proxy-SNP(s) for ", length(unique(as.factor(snps_reference))) -
                   length(unique(as.factor(proxies$target_snp.outcome))),
                 " missing-SNP(s) not available (e.g., no proxy-SNP or r2 < provided)"))
  message(paste0("# 3. extracting proxy-SNP(s) from outcome"))
  proxy_snps <- proxies %>% dplyr::distinct(proxy_snp.outcome) %>%
    dplyr::pull(proxy_snp.outcome)
  data_outcome_proxies <- fread(data_outcome_path,
                                header = FALSE, sep = "\t",
                                col.names = c("chr.outcome", "pos.outcome", "SNPID", "SNP", "A1", "A2",
                                              "beta.outcome", "pval.outcome", "min_log10_pval", "se.outcome",
                                              "samplesize.outcome", "ImpMAF", "outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome"))
  data_outcome_proxies <- data_outcome_proxies %>%
    dplyr::filter(SNP %in% proxy_snps)
  data_outcome_proxies <- dplyr::left_join(data_outcome_proxies,
                                           proxies, by = c(SNP = "SNP"), relationship = "many-to-many")
  message(paste0("## proxy-SNP(s) for ", length(unique(as.factor(data_outcome_proxies$target_snp.outcome))),
                 " of ", length(unique(as.factor(proxies$target_snp.outcome))),
                 " missing-SNP(s) extracted"))
  data_outcome_proxies <- data_outcome_proxies %>% dplyr::group_by(target_snp.outcome) %>%
    dplyr::filter(R == max(R)) %>% dplyr::slice(1) %>% dplyr::select(-R)
  data_outcome_proxies <- format_data2(
    data_outcome_proxies,
    type = "outcome",
    snps = NULL,
    header = TRUE,
    phenotype_col = "outcome",
    id_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    chr_col = "chr.outcome",
    pos_col = "pos.outcome",
    samplesize_col = "samplesize.outcome",
    min_pval = 1e-200,
    log_pval = FALSE
  )
  data_outcome <- dplyr::bind_rows(data_outcome, data_outcome_proxies)
  return(data_outcome)
}

#' finds proxy SNPs from an internal reference (UKB format)
#'
#' @param data_exposure Data frame containing exposure SNP data.
#' @param data_outcome Data frame containing outcome SNP data.
#' @param data_outcome_path Path to the outcome data file.
#' @param data_reference Data frame containing reference SNP data.
#' @param data_reference_path Path to the reference data file.
#' @param tag_r2 Threshold for LD (default: 0.8).
#' @param tag_kb Threshold for distance in kb (default: 5000).
#' @param tag_nsnp Threshold for number of SNPs (default: 5000).
#' @param outcome_sep Delimiter used in the outcome data file.
#' @return Data frame with proxy SNPs added.
#' @export
proxy_search_UKB <- function (data_exposure, data_outcome, data_outcome_path, data_reference,
                              data_reference_path, tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000,
                              outcome_sep)
{
  if (!file.exists(data_outcome_path) || !file.exists(paste0(data_reference))) {
    stop("Invalid file paths provided.")
  }
  if (!is.data.frame(data_exposure) || !is.data.frame(data_outcome)) {
    stop("data_exposure and data_outcome must be data frames.")
  }
  snps_missing <- setdiff(unique(as.factor(data_exposure$SNP)),
                          unique(as.factor(data_outcome$SNP)))
  if (length(snps_missing) == 0) {
    message("No missing SNPs. Returning the original dataframe.")
    return(data_outcome)
  }
  message(paste0("# 1. looking up ", length(unique(as.factor(snps_missing))),
                 " missing-SNP(s) in the reference panel"))
  reference_header <- readLines(data_reference, n = 1)
  column_index <- grep("rs", strsplit(reference_header, "\t")[[1]])
  cmd <- sprintf("awk '$%d ~ /^(%s)$/' %s", column_index,
                 paste(snps_missing, collapse = "|"), paste0(data_reference))
  reference <- data.table::fread(cmd = cmd, quote = "")
  snps_reference <- intersect(unique(as.factor(snps_missing)),
                              unique(as.factor(reference$V2)))
  message(paste0("## ", length(unique(as.factor(snps_reference))),
                 " of ", length(unique(as.factor(snps_missing))), " missing-SNP(s) available in the reference panel"))
  message(paste0("# 2. extracting proxy-SNP(s) for the ",
                 length(unique(as.factor(snps_reference))), " missing-SNP(s) from the reference panel"))
  gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
  proxies <- functions::get_ld_proxies(rsid = snps_missing,
                                       bfile = data_reference_path, searchspace = NULL, tag_kb = tag_kb,
                                       tag_nsnp = tag_nsnp, tag_r2 = tag_r2, threads = 1, out = tempfile())
  proxies <- proxies %>% dplyr::select(target_snp.outcome = SNP_A,
                                       proxy_snp.outcome = SNP_B, target_a1.outcome = A1, target_a2.outcome = A2,
                                       proxy_a1.outcome = B1, proxy_a2.outcome = B2, R) %>%
    dplyr::mutate(proxy.outcome = TRUE, SNP = proxy_snp.outcome) %>%
    dplyr::select(proxy.outcome, everything())
  message(paste0("## proxy-SNP(s) for ", length(unique(proxies$target_snp.outcome)),
                 " missing-SNP(s) found; ", "proxy-SNP(s) for ", length(unique(as.factor(snps_reference))) -
                   length(unique(as.factor(proxies$target_snp.outcome))),
                 " missing-SNP(s) not available (e.g., no proxy-SNP or r2 < provided)"))
  message(paste0("# 3. extracting proxy-SNP(s) from outcome"))
  proxy_snps <- proxies %>% dplyr::distinct(proxy_snp.outcome) %>%
    dplyr::pull(proxy_snp.outcome)
  data_outcome_proxies <- fread(data_outcome_path,
                                header = FALSE,
                                col.names = c("ID", "REF", "ALT", "SNP", "POS19", "POS38"))
  data_outcome_proxies <- data_outcome_proxies %>%
    dplyr::filter(SNP %in% proxy_snps)
  data_outcome_proxies <- separate(data_outcome_proxies, ID, into = c("chr.outcome", "pos.outcome", "ID", "other_allele.outcome", "effect_allele.outcome", "eaf.outcome",
                                                                      "INFO", "samplesize.outcome", "TEST", "beta.outcome", "se.outcome", "CHISQ", "LOG10P", "EXTRA", "outcome", "ID2"), sep = " ")
  data_outcome_proxies <- dplyr::left_join(data_outcome_proxies,
                                           proxies, by = c(SNP = "SNP"), relationship = "many-to-many")
  message(paste0("## proxy-SNP(s) for ", length(unique(as.factor(data_outcome_proxies$target_snp.outcome))),
                 " of ", length(unique(as.factor(proxies$target_snp.outcome))),
                 " missing-SNP(s) extracted"))
  data_outcome_proxies <- data_outcome_proxies %>% dplyr::group_by(target_snp.outcome) %>%
    dplyr::filter(R == max(R)) %>% dplyr::slice(1) %>% dplyr::select(-R)

  data_outcome_proxies <- format_data2(
    data_outcome_proxies,
    type = "outcome",
    snps = NULL,
    header = TRUE,
    phenotype_col = "outcome",
    id_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    chr_col = "chr.outcome",
    pos_col = "pos.outcome",
    samplesize_col = "samplesize.outcome",
    min_pval = 1e-200,
    log_pval = FALSE
  )

  data_outcome <- dplyr::bind_rows(data_outcome, data_outcome_proxies)
  return(data_outcome)
}

#' calculate r2
#'
#' @param data Data frame containing SNP data.
#' @return Data frame with r2 values added.
#' @export
calculate_r2 <- function(data) {
  r2 <- (2 * (data$beta^2) * data$eaf.exposure * (1 - data$eaf.exposure)) /
    ((2 * (data$beta^2) * data$eaf.exposure * (1 - data$eaf.exposure)) +
       ((data$se.exposure^2) * (2 * data$samplesize.exposure) * data$eaf.exposure * (1 - data$eaf.exposure)))

  data$r2 <- r2
  return(data)
}

#' calculate fstat from r2
#'
#' @param data Data frame containing SNP data.
#' @param grouping_column Column to group by.
#' @return Data frame with fstat values added.
#' @export
calculate_fstat_from_r2 <- function(data, grouping_column) {
  # Group by the specified column and calculate fstat for each group
  result <- data %>%
    group_by({{ grouping_column }}) %>%
    mutate(k = n_distinct(SNP),  # Number of unique SNPs
           fstat_from_r2 = r2 * (samplesize.exposure - 1 - k) / ((1 - r2) * k)) %>%
    ungroup() %>%
    select(-k) %>%  # Remove the auxiliary column
    group_by({{ grouping_column }}) %>%
    summarize(mean_fstat_from_r2 = mean(fstat_from_r2, na.rm = TRUE)) %>%
    left_join(data, by = {{ grouping_column }})
  return(result)
}

#' calculate fstat from beta se
#'
#' @param data Data frame containing SNP data.
#' @return Data frame with fstat values added.
#' @export
calculate_fstat_from_beta_se <- function(data) {
  f_stats <- (data$beta.exposure / data$se.exposure)^2
  data$fstat_from_b_se <- f_stats
  # Calculate mean f_stats for each level of id.exposure
  fstat_mean_from_b_se <- data %>%
    group_by(id.exposure) %>%
    summarize(mean_f_stats = mean(f_stats, na.rm = TRUE)) %>%
    ungroup()
  # Merge the mean_f_stats column with the original data
  data <- left_join(data, mean_f_stats, by = "id.exposure")
  return(data)
}

#' add the mean sample size for missing SNP sample sizes
#'
#' @param data Data frame containing SNP data.
#' @param grouping_column Column to group by.
#' @param column_name Name of the column with missing values.
#' @return Data frame with missing values replaced by mean.
#' @export
replace_na_with_mean <- function(data, grouping_column, column_name) {
  mean_values <- data %>%
    dplyr::group_by({{ grouping_column }}) %>%
    dplyr::summarize(mean_value = mean({{ column_name }}, na.rm = TRUE)) %>%
    dplyr::filter(!is.infinite(mean_value))  # Filter out Inf values

  result <- data %>%
    dplyr::left_join(mean_values, by = {{ grouping_column }}) %>%
    dplyr::mutate({{ column_name }} := ifelse(is.na({{ column_name }}), mean_value, {{ column_name }})) %>%
    dplyr::select(-mean_value)

  return(result)
}

#' Fill Missing EAF Values in Data Frame
#'
#' This function checks if the specified column containing EAF values in the input dataframe contains
#' any missing values (NA). If missing values are found, it reads the necessary
#' reference data to fill in the missing values based on the provided allele
#' information.
#'
#' @param df The input dataframe containing the outcome data.
#' @param reference The file path to the reference data containing EAF information.
#' @param column_EAF The name of the column containing EAF values. If no column then this will be the name of the new column.
#' @return The input dataframe with missing EAF values filled.
#' @export
missing_EAF <- function(df, reference, column_EAF) {
  if (!(column_EAF %in% colnames(df)) || any(is.na(df[[column_EAF]]))) {
    if (!(column_EAF %in% colnames(df))) {
      df[[column_EAF]] <- NA
    }
    EAF <- fread(input = reference,
                 header = TRUE,
                 select = c("Predictor", "A1", "A2", "MAF"), # Select only necessary columns
                 data.table = FALSE) %>%
      dplyr::filter(Predictor %in% df$SNP) %>%
      dplyr::rename(EAF = MAF)

    df <- df %>%
      dplyr::left_join(EAF, by = c("SNP" = "Predictor")) %>%
      dplyr::mutate(!!column_EAF := ifelse(is.na(!!rlang::sym(column_EAF)),
                                    ifelse(effect_allele.outcome == A1 & other_allele.outcome == A2, EAF,
                                           ifelse(effect_allele.outcome == A2 & other_allele.outcome == A1, 1 - EAF, !!rlang::sym(column_EAF))),
                                    !!rlang::sym(column_EAF))) %>%
      dplyr::select(-A1, -A2, -EAF)
  }
  return(df)
}
