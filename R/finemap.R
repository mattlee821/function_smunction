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
#'   \item \code{label} The SNP identifier for the SNP with the highest PIP value within each credible set (NA for other SNPs).
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
  PIP <- data.frame(SNP = names(susieR_model$pip), pip_value = susieR_model$pip)

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
    dplyr::mutate(label = if_else(P == min(P), SNP, NA_character_)) %>%  # Identify SNP with smallest p-value
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
#'   \item label: SNP identifier of the SNP with the smallest p-value in each credible set.
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

#' Generate a Regional Association Plot with Recombination Rates and Genes
#'
#' This function creates a regional association plot showing -log10(p-values)
#' for SNPs, overlaid with recombination rates and gene annotations. Points can be
#' highlighted and labeled based on the SNP of interest.
#'
#' @param df Dataframe containing GWAS summary statistics.
#' @param rsid Column name for SNP identifiers in `df`.
#' @param chrom Column name for chromosome numbers in `df`.
#' @param pos Column name for SNP positions (in base pairs) in `df`.
#' @param p_value Column name for p-values in `df`.
#' @param label Optional vector of SNP IDs (rsIDs) to label in the plot. Default is NULL.
#' @param trait Column name for trait or phenotype in `df`. Default is NULL.
#' @param plot_pvalue_threshold Minimum p-value threshold for points to be plotted. Default is 0.1.
#' @param genome_build Genome build version for gene and recombination rate annotations. Default is "GRCh38".
#' @param population Population for recombination rate data (e.g., "EUR"). Default is "EUR".
#' @param plot_title Title of the plot. Default is NULL.
#' @param plot_subtitle Subtitle of the plot. Default is NULL.
#'
#' @return A combined `ggplot2` object with the regional association plot and gene annotations.
#' @export
#'
#' @examples
#' \dontrun{
#' gg_regionplot(df = gwas_data, rsid = "rsid", chrom = "CHR", pos = "POS",
#'               p_value = "P", label = c("rs12345"), plot_title = "Region Plot")
#' }
gg_regionplot <- function(df,
                          rsid,
                          chrom,
                          pos,
                          p_value,
                          label = NULL,
                          trait = NULL,
                          plot_pvalue_threshold = 0.1,
                          genome_build = "GRCh38",
                          population = "EUR",
                          plot_title = NULL,
                          plot_subtitle = NULL) {

  # Process and clean input dataframe ====
  df <- df %>%
    mutate(log10_pval = -log10(df[[p_value]])) %>%  # Compute -log10(p-value)
    dplyr::select(rsid = SNP, chromosome = CHR, position = pos, ref = EA, alt = OA, log10_pval, trait = phenotype) %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character)) %>%  # Convert factor columns to characters
    dplyr::mutate(ref = stringr::str_to_upper(ref), alt = stringr::str_to_upper(alt)) %>%  # Normalize allele case
    dplyr::group_by(trait, rsid) %>%
    dplyr::slice_max(log10_pval) %>%  # Keep the SNP with the highest log10(p-value) per trait
    dplyr::ungroup() %>%
    drop_na()

  # Extract recombination rate data ====
  ylim <- max(df$log10_pval, na.rm = TRUE) + 0.3 * max(df$log10_pval, na.rm = TRUE)  # Adjust y-axis limit
  recomb_df <- recomb_extract_locuszoom(
    chrom = unique(df$chromosome),
    start = min(df$position) - 1,
    end = max(df$position) + 1,
    genome_build = genome_build
  ) %>%
    select(position, recomb_rate)

  # Generate the regional association plot ====
  plot_region <- df %>%
    distinct(rsid, .keep_all = TRUE) %>%
    filter(log10_pval > -log10(plot_pvalue_threshold)) %>%
    mutate(
      size = ifelse(rsid %in% label, 6, 4),   # Larger size for labeled points
      color = ifelse(rsid %in% label, "purple", "black"),  # Purple color for labeled points
      alpha = ifelse(rsid %in% label, 1, 0.7),  # Higher opacity for labeled points
      shape = ifelse(rsid %in% label, 18, 16)  # Different shape for labeled points
    ) %>%
    ggplot(aes(
      x = position / 1e6,  # Convert position to Mb
      y = log10_pval)) +

    # Add recombination rate line (background layer)
    geom_line(
      data = recomb_df,
      mapping = aes(
        x = position / 1e6,
        y = recomb_rate),
      color = "lightblue",
      linewidth = 0.5
    ) +

    # Plot non-label points
    geom_point(data = subset(df, !rsid %in% label),
               aes(size = 4, color = "darkgrey", alpha = 0.7, shape = 16)) +

    # Plot label points
    geom_point(data = subset(df, rsid %in% label),
               aes(size = 6, color = "purple", alpha = 1, shape = 18)) +

    # Set scales for aesthetic mappings
    scale_size_identity() +
    scale_color_identity() +
    scale_alpha_identity() +
    scale_shape_identity() +

    # Add labels and theme
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Position (Mb)",
      y = "-log10(P)") +
    scale_x_continuous(labels = scales::comma) +
    cowplot::theme_cowplot() +

    # Add point labels
    ggrepel::geom_label_repel(
      data = subset(df, rsid %in% label),
      aes(label = rsid),
      size = 4,
      color = "black",
      fontface = "bold",
      fill = "white",
      min.segment.length = 0,
      box.padding = 1,
      alpha = 1,
      max.overlaps = 100,
      force = 10
    ) +

    # Add horizontal significance threshold line
    geom_hline(yintercept = -log10(5E-8), linetype = "dashed")

  # Add recombination rate as a secondary axis
  plot_region <- plot_region +
    scale_y_continuous(
      name = "-log10(P)",
      limits = c(0, ylim),
      sec.axis = sec_axis(
        ~ . * (100 / ylim),  # Scale secondary axis
        name = "Recombination rate (cM/Mb)"
      )
    ) +
    theme(axis.title.y.right = element_text(vjust = 1.5))

  # Generate the gene plot ====
  plot_gene <- gg_geneplot(
    chr = unique(df$chromosome),
    start = min(df$position) - 1,
    end = max(df$position) + 1,
    genome_build = genome_build
  ) +
    theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

  # Combine the region and gene plots ====
  plot_region <- patchwork::wrap_plots(list(
    plot_region +
      labs(x = "") +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(5.5, 5.5, 0, 5.5)
      ),
    plot_gene
  ),
  nrow = 2,
  heights = c(3, 1)
  )

  # Return the combined plot ====
  return(plot_region)
}
