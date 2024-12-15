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
#' @param labels Logical for whether to show text labels for `label` SNPs
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
                          labels = TRUE,
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

  # Conditionally add ggrepel labels ====
  if (labels) {
    plot_region <- plot_region +
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
      )
  }

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

#' @title Locus Plot for Genomic Regions
#' @description Creates a locus plot for visualizing genomic regions, highlighting SNPs, p-values, and optionally LD information.
#' @param df Dataframe containing the data to be plotted.
#' @param lead_snp Character, optional. The lead SNP to be highlighted. Defaults to the SNP with the lowest p-value if not provided.
#' @param rsid Column name for SNP IDs in `df`.
#' @param chrom Column name for chromosome information in `df`.
#' @param pos Column name for SNP position in `df`.
#' @param ref Column name for reference allele in `df`.
#' @param alt Column name for alternate allele in `df`.
#' @param effect Column name for effect size in `df`. Used to calculate p-values if `std_err` is provided.
#' @param std_err Column name for standard error in `df`. Used to calculate p-values if `effect` is provided.
#' @param p_value Column name for p-values in `df`. Used if `effect` and `std_err` are not provided.
#' @param trait Column name for traits in `df`, optional.
#' @param plot_pvalue_threshold Numeric, p-value threshold for plotting. Default is 0.1.
#' @param plot_subsample_prop Numeric, proportion of SNPs to subsample for plotting. Default is 0.25.
#' @param plot_distance Numeric, distance in base pairs to include around the lead SNP. Default is 500,000.
#' @param genome_build Character, genome build to use (e.g., "GRCh37"). Default is "GRCh37".
#' @param population Character, population for LD calculations. Default is "ALL".
#' @param plot_genes Logical, whether to include genes in the plot. Default is FALSE.
#' @param plot_recombination Logical, whether to include recombination rates in the plot. Default is FALSE.
#' @param plot_title Character, optional title for the plot.
#' @param plot_subtitle Character, optional subtitle for the plot.
#' @param label Character vector, optional SNP IDs to label in the plot.
#' @return A ggplot2 object representing the locus plot.
#' @examples
#' \dontrun{
#' gg_locusplot(df = my_data, rsid = rsid_col, chrom = chrom_col, pos = pos_col, ref = ref_col, alt = alt_col)
#' }
#' @export
gg_locusplot <- function(df,
                         lead_snp = NULL,
                         rsid = rsid,
                         chrom = chrom,
                         pos = pos,
                         ref = ref,
                         alt = alt,
                         effect = NULL,
                         std_err = NULL,
                         p_value = p_value,
                         trait = NULL,
                         plot_pvalue_threshold = 0.1,
                         plot_subsample_prop = 0.25,
                         plot_distance = 500000,
                         genome_build = "GRCh37",
                         population = "ALL",
                         plot_genes = FALSE,
                         plot_recombination = FALSE,
                         plot_title = NULL,
                         plot_subtitle = NULL,
                         label = NULL) {
  # Check input arguments to ensure they are of the correct type and within reasonable ranges
  checkmate::assert_data_frame(df)
  # checkmate::assert_string(lead_snp)
  checkmate::assert_numeric(plot_pvalue_threshold, upper = 1)
  checkmate::assert_numeric(plot_subsample_prop, lower = 0, upper = 1)
  checkmate::assert_numeric(plot_distance, lower = 0)
  checkmate::assert_logical(plot_genes)
  if (!is.null(label)) {
    checkmate::assert_character(label)
  }

  if(!rlang::quo_is_null(rlang::enquo(effect)) & !rlang::quo_is_null(rlang::enquo(std_err))) {
    checkmate::assert_numeric(df %>% pull({{ effect }}))
    checkmate::assert_numeric(df %>% pull({{ std_err }}))

    df <- df %>%
      rename(.effect = {{ effect }},
             .std_err = {{ std_err }}) %>%
      mutate(log10_pval = abs((pnorm(-abs(.effect/.std_err), log.p=TRUE) + log(2)) / log(10)))
  } else {
    df <- df %>%
      mutate(log10_pval = -log10({{ p_value }}))
  }

  if (rlang::quo_is_null(rlang::enquo(trait))) {
    df <- df %>%
      select(rsid = {{ rsid }}, chromosome = {{ chrom }}, position = {{ pos }}, ref = {{ ref }}, alt = {{ alt }}, log10_pval) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(ref = stringr::str_to_upper(ref), alt = stringr::str_to_upper(alt)) %>%
      group_by(rsid) %>%
      slice_max(log10_pval) %>%
      ungroup() %>%
      tidyr::drop_na()
  } else {
    df <- df %>%
      select(rsid = {{ rsid }}, chromosome = {{ chrom }}, position = {{ pos }}, ref = {{ ref }}, alt = {{ alt }}, log10_pval, trait = {{ trait }}) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(ref = stringr::str_to_upper(ref), alt = stringr::str_to_upper(alt)) %>%
      group_by(trait, rsid) %>%
      slice_max(log10_pval) %>%
      ungroup() %>%
      tidyr::drop_na()
  }

  # Create df containing information about lead SNP (by default, select SNP with lowest p-value, otherwise take user-supplied value)
  if (is.null(lead_snp)) {
    indep_snps <- df %>%
      slice_max(log10_pval, with_ties = FALSE, n = 1) %>%
      select(lead_rsid = rsid, lead_chromosome = chromosome, lead_position = position, lead_ref = ref, lead_alt = alt)

    cli::cli_alert_info("No lead_snp supplied. Defaulting to {indep_snps$lead_rsid} - {indep_snps$lead_chromosome}:{indep_snps$lead_position}, which has the lowest p-value in the region")
  } else if (!(lead_snp %in% df$rsid)) {
    # ensure lead_snp is in the supplied data; if not, use minimum p-value at locus
    indep_snps <- df %>%
      slice_max(log10_pval, with_ties = FALSE, n = 1) %>%
      select(lead_rsid = rsid, lead_chromosome = chromosome, lead_position = position, lead_ref = ref, lead_alt = alt)

    cli::cli_alert_info("Lead snp not present in supplied locus data. Defaulting to {indep_snps$lead_rsid} - {indep_snps$lead_chromosome}:{indep_snps$lead_position}, which has the lowest p-value in the region")
  } else {
    indep_snps <- df %>%
      select(lead_rsid = rsid, lead_chromosome = chromosome, lead_position = position, lead_ref = ref, lead_alt = alt) %>%
      filter(lead_rsid == lead_snp) %>%
      distinct(lead_rsid, .keep_all = TRUE)
  }

  # Create dataframe of variants within the region size specified by the user
  suppressMessages(locus_snps <- df %>%
                     filter(rsid %in% indep_snps$lead_rsid) %>%
                     select(chromosome, position, lead_rsid = rsid) %>%
                     purrr::pmap_dfr(function(chromosome_filter = first, position_filter = second, lead_rsid = third) {
                       df %>%
                         filter(chromosome == chromosome_filter & between(position, position_filter - plot_distance / 2, position_filter + plot_distance / 2)) %>%
                         mutate(lead_rsid = lead_rsid) %>%
                         left_join(indep_snps)
                     }))

  # Extract LD and format colors
  possibly_ld_extract_locuszoom <- purrr::possibly(locusplotr::ld_extract_locuszoom, otherwise = NULL)

  ld_extracted <- possibly_ld_extract_locuszoom(chrom = indep_snps$lead_chromosome, pos = indep_snps$lead_position, ref = indep_snps$lead_ref, alt = indep_snps$lead_alt, start = min(locus_snps$position), stop = max(locus_snps$position), genome_build = genome_build, population = population)

  # Create dataframe with variants at locus, LD information, color codes, and labels in preparation for plotting
  if (!(is.null(ld_extracted))) {
    # Join locus df with LD information
    locus_snps_ld <- ld_extracted %>%
      select(chromosome = chromosome2, position = position2, variant2, correlation) %>%
      mutate(chromosome = as.numeric(chromosome), position = as.numeric(position)) %>%
      tidyr::separate(variant2, into = c("chr_pos", "ref_alt"), sep = "_") %>%
      tidyr::separate(ref_alt, into = c("ref", "alt"), sep = "/") %>%
      right_join(locus_snps, by = c("chromosome" = "chromosome", "position" = "position"), relationship = "many-to-many") %>%
      filter((ref.x == ref.y & alt.x == alt.y) | (ref.x == alt.y & alt.x == ref.y)) %>%
      select(-ends_with(".y"), -chr_pos) %>%
      rename_with(~ stringr::str_replace(.x, ".x", ""), .cols = ends_with(".x"))

    # Create color codes and labels
    locus_snps_ld <- locus_snps_ld %>%
      mutate(color_code = as.character(cut(as.numeric(correlation), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("blue4", "skyblue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
      mutate(legend_label = as.character(cut(as.numeric(correlation), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1"), include.lowest = TRUE))) %>%
      mutate(lead = rsid == lead_rsid) %>%
      mutate(label = case_when(
        rsid == lead_rsid ~ lead_rsid,
        rsid %in% label ~ rsid,
        TRUE ~ NA_character_
      )) %>%
      mutate(color_code = case_when(
        rsid == lead_rsid ~ "purple",
        TRUE ~ color_code
      )) %>%
      mutate(color_code = forcats::fct_expand(color_code, "purple", "red", "orange", "darkgreen", "skyblue", "blue4")) %>%
      mutate(color_code = forcats::fct_relevel(color_code, "purple", "red", "orange", "darkgreen", "skyblue", "blue4")) %>%
      mutate(legend_label = case_when(
        rsid == lead_rsid ~ "Ref",
        TRUE ~ legend_label
      )) %>%
      mutate(legend_label = forcats::fct_expand(legend_label, "Ref", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2")) %>%
      mutate(legend_label = forcats::fct_relevel(legend_label, "Ref", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"))
  } else {
    # Deal with scenario where lead variant is not present in LD database
    cli::cli_alert_info("No linkage disequilibrium information found")
    locus_snps_ld <- locus_snps %>%
      mutate(correlation = NA_integer_) %>%
      mutate(lead = rsid == lead_rsid) %>%
      mutate(label = case_when(
        rsid == lead_rsid ~ lead_rsid,
        TRUE ~ NA_character_
      )) %>%
      mutate(color_code = case_when(
        rsid == lead_rsid ~ "purple",
        TRUE ~ "grey50"
      )) %>%
      mutate(legend_label = case_when(
        rsid == lead_rsid ~ "Ref",
        TRUE ~ "Other"
      ))
  }

  # group locus by trait if necessary
  if (!rlang::quo_is_null(rlang::enquo(trait))) {
    locus_snps_ld <- locus_snps_ld %>%
      group_by(.data = ., trait)
  }

  locus_snps_ld_label <- locus_snps_ld %>% filter(!is.na(label))

  # Make plot (sample non-significant p-values to reduce overplotting)
  regional_assoc_plot <- locus_snps_ld %>%
    distinct(rsid, .keep_all = TRUE) %>%
    filter(log10_pval > -log10(plot_pvalue_threshold) | correlation > 0.2 | legend_label == "Ref") %>% # improve overplotting
    bind_rows(locus_snps_ld %>%
                filter(log10_pval <= -log10(plot_pvalue_threshold) & correlation < 0.2 & legend_label != "Ref") %>%
                slice_sample(prop = plot_subsample_prop)) %>%
    arrange(color_code, log10_pval) %>%
    ggplot(aes(position, log10_pval)) +
    geom_point(aes(fill = factor(color_code), size = lead, alpha = lead, shape = lead)) +
    ggrepel::geom_label_repel(data = locus_snps_ld_label, aes(label = label),
                              size = 4,
                              color = "black",
                              fontface = "bold",
                              fill = "white",
                              min.segment.length = 0,
                              box.padding = 1,
                              alpha = 1,
                              nudge_y = 4
    ) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
    scale_fill_identity(parse(text = "r^2"), guide = "legend", labels = levels(forcats::fct_drop(locus_snps_ld$legend_label)), na.translate = FALSE) +
    scale_size_manual(values = c(3, 5), guide = "none") +
    scale_shape_manual(values = c(21, 23), guide = "none") +
    scale_alpha_manual(values = c(0.8, 1), guide = "none") +
    scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 6),
                               position = "inside")) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = glue::glue("Position on Chromosome {unique(indep_snps$lead_chromosome)} (Mb)"),
      y = "-log<sub>10</sub>(P-value)"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      # legend.text = element_text(size = 10),
      # legend.title = element_text(size = 10, hjust = 0.5),
      # legend.justification.inside = c("right", "top"),
      # legend.position.inside = c(0.99, 0.99),
      # legend.spacing.y = unit(0, "pt"),
      strip.text = element_text(color = "black"),
      strip.text.x = element_blank(),
      axis.title.y = ggtext::element_markdown()
    )

  if (plot_recombination) {
    cli::cli_alert_info("Extracting recombination rates for the region {indep_snps$lead_chromosome}:{indep_snps$lead_position - plot_distance/2}-{indep_snps$lead_position + plot_distance/2}")
    ylim <- max(pull(locus_snps_ld, log10_pval), na.rm = TRUE) +
      0.3 * max(pull(locus_snps_ld, log10_pval), na.rm = TRUE)

    recomb_df <- recomb_extract_locuszoom(chrom = indep_snps$lead_chromosome, start = indep_snps$lead_position - plot_distance / 2, end = indep_snps$lead_position + plot_distance / 2, genome_build = genome_build) %>%
      select(position, recomb_rate)
    # return(recomb_df)
    suppressMessages(
      regional_assoc_plot <- regional_assoc_plot +
        geom_line(data = recomb_df, mapping = aes(x = position, y = recomb_rate), color = "lightblue", linewidth = 0.5) +
        scale_y_continuous(
          name = "-log<sub>10</sub>(P-value)",
          limits = c(0, ylim),
          sec.axis = sec_axis(
            ~. * (100 / ylim),
            name = "Recombination rate (cM/Mb)"
            # labels = scales::percent_format()
          )
        ) +
        theme(axis.title.y.right = element_text(vjust = 1.5))
    )

    regional_assoc_plot <- gginnards::move_layers(regional_assoc_plot, "GeomLine", "bottom")

  }

  # Add plot of genes if requested by user
  if (plot_genes) {
    cli::cli_alert_info("Extracting genes for the region {indep_snps$lead_chromosome}:{indep_snps$lead_position - plot_distance/2}-{indep_snps$lead_position + plot_distance/2}")
    geneplot <- gg_geneplot(chr = indep_snps$lead_chromosome, start = indep_snps$lead_position - plot_distance / 2, end = indep_snps$lead_position + plot_distance / 2, genome_build = genome_build) +
      theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

    suppressWarnings(suppressMessages(regional_assoc_plot <- patchwork::wrap_plots(list(
      regional_assoc_plot +
        labs(x = "") +
        xlim(indep_snps$lead_position - plot_distance / 2, indep_snps$lead_position + plot_distance / 2) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 0, 5.5)
        ),
      geneplot
    ), nrow = 2, heights = c(3, 1))))
  }

  return(regional_assoc_plot)
}
