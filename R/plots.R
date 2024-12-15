# Function to plot regional association with LD
#' Create a regional association plot
#' The code is adapted from the [locusplotr package](https://github.com/mglev1n/locusplotr).
#'
#' Returns a ggplot object containing a regional association plot (-log10(p-value) as a function of chromosomal position, with variants colored by linkage disequilibrium to reference variant).
#' This function allows the user to integrate genome wide association study (GWAS) summary statistics for a locus of interest with linkage disequilibrium information (obtained using the University of Michigan LocusZoom API <https://portaldev.sph.umich.edu/>) for that locus to create a regional association plot.
#'
#' @param df Dataframe containing columns with rsid, chromosome, position, reference/effect allele, alternate/non-effect allele, and p-value for all variants within the range of interest
#' @param lead_snp A character vector containing a lead variant of interest. When NULL (default), the variant with the lowest p-value will be selected as the lead variant.
#' @param rsid Rsid column
#' @param chrom Chromosome column
#' @param pos Position column
#' @param ref Reference/effect allele column
#' @param alt Alternate/non-effect allele column
#' @param effect Effect size column (on beta or log-odds scale)
#' @param std_err Standard error column
#' @param p_value P-value column
#' @param plot_pvalue_threshold Threshold for plotting p-value on regional association plot (default = 0.1) - reducing the number of points decreases file size and improves performance
#' @param plot_subsample_prop Proportion of points above p-value threshold to plot (default = 0.25; range = 0-1) - reducing the number of points decreases file size and improves performance
#' @param plot_distance Integer corresponding to the size of the locus that should be plotted
#' @param genome_build Character - one of "GRCh37" or "GRCh38"
#' @param population Character - one of "ALL", "AFR", "AMR", "EAS", "EUR", "SAS" referring to the reference population of interest for obtaining linkage disequilibrium information (default = "ALL")
#' @param plot_genes Logical - Include a plot of genes/transcripts within the region of interest beneath the regional association plot (default = FALSE)
#' @param plot_recombination Logical - Include a secondary y-axis of recombination rate within the region of interest
#' @param plot_title A character string corresponding to plot title (default = NULL)
#' @param plot_subtitle A character string corresponding to plot subtitle (default = NULL)
#' @param trait (optional) Column containing the name of the trait
#' @param label Labels for SNPs to highlight (optional).
#'
#' @return A ggplot object containing a regional association plot for the locus of interest
#'
#' @examples
#' \dontrun{
#' # Basic regional association plot
#' gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value)
#'
#' # Use "plot_genes = TRUE" to add a plot of genes within the region beneath the regional association plot
#' gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = TRUE)
#' }
#'
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
                         population = "EUR",
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

#' Create a Region Plot with Genetic Data and Recombination Rates
#'
#' This function generates a region plot with genetic data for a given chromosome and
#' position range, and overlays it with recombination rates. Significant genetic
#' variants are highlighted with labels, and a gene plot is included below the main plot.
#'
#' @param df A data frame containing genetic data. It must include the following columns:
#'   - `SNP`: SNP identifier (rsid).
#'   - `CHR`: Chromosome number.
#'   - `pos`: Position of the SNP (e.g., in base pairs).
#'   - `EA`: Effect allele (ref).
#'   - `OA`: Other allele (alt).
#'   - `phenotype`: The trait or phenotype associated with the SNP.
#'   - `p_value`: P-value for the association.
#' @param rsid A character vector of SNP identifiers to label on the plot.
#' @param chrom The chromosome number to be used for plotting.
#' @param pos The column name that contains the SNP positions.
#' @param p_value The column name containing the p-values for each SNP.
#' @param label A character vector of SNPs to label on the plot. Default is `NULL`.
#' @param labels Logical for whether to show text labels for `label` SNPs
#' @param trait The trait associated with the SNP. Default is `NULL`.
#' @param plot_pvalue_threshold The p-value threshold for plotting. SNPs with p-values greater than this threshold are excluded from the plot. Default is `0.1`.
#' @param genome_build The genome build used for the data (e.g., "GRCh38"). Default is `"GRCh38"`.
#' @param population The population for which the data was analyzed. Default is `"EUR"`.
#' @param plot_title The title of the plot. Default is `NULL`.
#' @param plot_subtitle The subtitle of the plot. Default is `NULL`.
#'
#' @return A `ggplot2` object containing the combined plot of region data with genetic variants, recombination rates, and gene annotations.
#'
#' @examples
#' # Example usage
#' \dontrun{
#' gg_regionplot(df = my_data,
#'               rsid = c("rs12345", "rs67890"),
#'               chrom = 1,
#'               pos = "position",
#'               p_value = "p_value",
#'               label = c("rs12345"),
#'               plot_title = "Region Plot Example")
#'}
#'
#' @export
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
