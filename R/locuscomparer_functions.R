#' Create a merged data frame for locus comparison.
#' This function creates and merges data frames from exposure and outcome data,
#' adds labels, and joins linkage disequilibrium (LD) data.
#'
#' @import dplyr
#' @import reshape2
#' @param data_coloc_exposure (data.frame) Data frame containing exposure data with columns `snp`, `pval`, `position`, and `LD`.
#' @param data_coloc_outcome (data.frame) Data frame containing outcome data with columns `snp` and `pval`.
#' @param SNP_causal_exposure (string) The causal SNP identifier for the exposure data.
#' @return A merged data frame ready for locus comparison plotting.
#' @export
locuscomparer_make_df <- function(data_coloc_exposure, data_coloc_outcome, SNP_causal_exposure) {
  # Create data frames
  df1 <- data.frame(
    rsid = data_coloc_exposure$snp,
    pval = data_coloc_exposure$pval,
    logp = -log10(data_coloc_exposure$pval),
    pos = data_coloc_exposure$position
  )

  df2 <- data.frame(
    rsid = data_coloc_outcome$snp,
    pval = data_coloc_outcome$pval,
    logp = -log10(data_coloc_outcome$pval)
  )

  # Merge data frames
  df <- merge(df1, df2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
  df <- locuscomparer_add_label(df, SNP_causal_exposure)

  # Replace Inf values with 10% larger than the maximum finite value
  max_logp <- max(df$logp1[is.finite(df$logp1)])
  V1 <- max_logp * 1.1
  df$logp1[is.infinite(df$logp1)] <- V1
  max_logp <- max(df$logp2[is.finite(df$logp2)])
  V1 <- max_logp * 1.1
  df$logp2[is.infinite(df$logp2)] <- V1

  # Convert LD to df
  df_ld <- reshape2::melt(data_coloc_exposure$LD)
  colnames(df_ld) <- c("SNP_A", "SNP_B", "R2")
  df_ld <- df_ld %>%
    dplyr::filter(SNP_A == SNP_causal_exposure)

  # Join LD data
  df <- df %>%
    dplyr::left_join(df_ld, by = c("rsid" = "SNP_B"))

  # make VARS
  colour = locuscomparer::assign_color(rsid = df$rsid, snp = SNP_causal_exposure, ld = df_ld)
  shape = ifelse(df$rsid == SNP_causal_exposure, 23, 21)
  names(shape) = df$rsid
  size = ifelse(df$rsid == SNP_causal_exposure, 3, 2)
  names(size) = df$rsid

  return(list(df = df, df_ld = df_ld, colour = colour, shape = shape, size = size))
}

#' Add Labels to df Data Frame
#'
#' This function adds labels to a df data frame based on the presence of specified SNPs.
#'
#' @param df (data.frame) A data frame containing at least the column `rsid`.
#' @param snp (character vector) A vector of SNP rsIDs to be labeled.
#' @return (data.frame) The input data frame with an additional column `label` that contains the rsIDs
#' of the specified SNPs and empty strings for others.
#' @export
locuscomparer_add_label <- function(df, snp) {
  df$label <- ifelse(df$rsid %in% snp, df$rsid, '')
  return(df)
}

#' Generate a combined plot with two LocusZoom plots and a LocusCompare plot using coloc data
#'
#' This function creates a combined plot with two LocusZoom plots and a LocusCompare plot. Each LocusZoom plot represents an association study.
#'
#' @import ggplot2
#' @import cowplot
#' @param data_coloc_exposure (data.frame) Data frame containing the exposure data with the following columns:
#' \itemize{
#'   \item \code{snp} (character) - SNP rsID.
#'   \item \code{pval} (numeric) - p-value for the association.
#'   \item \code{chr} (character) - Chromosome name.
#'   \item \code{position} (integer) - Position on the chromosome.
#'   \item \code{LD} (matrix) - LD matrix for the SNPs.
#' }
#' @param data_coloc_outcome (data.frame) Data frame containing the outcome data with the following columns:
#' \itemize{
#'   \item \code{snp} (character) - SNP rsID.
#'   \item \code{pval} (numeric) - p-value for the association.
#' }
#' @param SNP_causal_exposure (character vector) A vector of SNP rsIDs to be labeled.
#' @param trait1_title (character) The title for the first trait (exposure). Default is "exposure".
#' @param trait2_title (character) The title for the second trait (outcome). Default is "outcome".
#' @return (ggplot) Combined plot of the two LocusZoom plots and the LocusCompare plot.
#' @export
locuscomparer_make_plot <- function(data_coloc_exposure, data_coloc_outcome, SNP_causal_exposure, trait1_title = "exposure", trait2_title = "outcome") {

  # make data
  list <- locuscomparer_make_df(data_coloc_exposure = data_coloc_exposure, data_coloc_outcome = data_coloc_outcome, SNP_causal_exposure = SNP_causal_exposure)
  df <- list[["df"]]
  df_ld <- list[["df_ld"]]
  colour <- list[["colour"]]
  shape <- list[["shape"]]
  size <- list[["size"]]

  # make the combined plot
  p <- locuscomparer_make_combined_plot(df = df,
                                        trait1_title = trait1_title, trait2_title = trait2_title,
                                        ld = df_ld,
                                        colour = colour,
                                        shape = shape,
                                        size = size)

  return(p)
}

#' Generate a combined plot with two LocusZoom plots and a LocusCompare plot using coloc data.
#'
#' This function creates a combined plot featuring a LocusCompare plot alongside two LocusZoom plots. Each LocusZoom plot represents an association study for a specific trait.
#' The LocusCompare plot compares the -log10(p-values) of two traits, while the LocusZoom plots show the detailed locus zoom views for the two traits.
#'
#' @param df (data.frame) An input data.frame containing the following columns:
#' \itemize{
#'   \item \code{rsid} (character) - SNP rsID.
#'   \item \code{pval1} (numeric) - p-value for study 1.
#'   \item \code{logp1} (numeric) - -log10(p-value) for study 1.
#'   \item \code{logp2} (numeric) - -log10(p-value) for study 2.
#'   \item \code{chr} (character) - Chromosome name (without 'chr').
#'   \item \code{pos} (integer) - Position on the chromosome.
#'   \item \code{label} (character) - Optional label for SNPs.
#' }
#' See the example for `get_lead_snp()` for details on how to generate this data.frame.
#'
#' @param trait1_title (string) Title for the x-axis of the LocusCompare plot, representing the first trait (exposure). Default is `"exposure"`.
#' @param trait2_title (string) Title for the y-axis of the LocusCompare plot, representing the second trait (outcome). Default is `"outcome"`.
#' @param ld (data.frame) The output from `retrieve_LD()`. This data.frame should contain columns:
#' \itemize{
#'   \item \code{SNP_A} (character) - SNP rsID for the reference SNP.
#'   \item \code{SNP_B} (character) - SNP rsID for the comparison SNP.
#'   \item \code{R2} (numeric) - LD (R^2) value between the reference and comparison SNPs.
#' }
#' @param colour (named vector) A vector specifying the colors for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = "#FF0000", "rs2" = "#0000FF", "rs3" = "#00FF00")`.
#' @param shape (named vector) A vector specifying the shapes for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = 21, "rs2" = 22, "rs3" = 23)`.
#' @param size (named vector) A vector specifying the sizes for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = 2, "rs2" = 3, "rs3" = 4)`.
#' @param snp (string, optional) SNP rsID to be highlighted. If NULL, the function will select the lead SNP. Default: NULL.
#' @param combine (boolean, optional) Should the three plots be combined into one plot? If FALSE, a list of three plots will be returned. Default: TRUE.
#' @param legend (boolean, optional) Should the legend be shown in the LocusCompare plot? Default: TRUE.
#' @param legend_position (string, optional) Position of the legend in the LocusCompare plot. Choose from 'bottomright', 'topright', or 'topleft'. Default: 'bottomright'.
#'
#' @return (ggplot) A combined plot with a LocusCompare plot and two LocusZoom plots, or a list of the three individual plots if `combine` is set to FALSE.
#' @export
locuscomparer_make_combined_plot = function (df, trait1_title = "exposure", trait2_title = "outcome",
                                             ld, colour = colour, shape = shape, size = size,
                                             snp = NULL,
                                             combine = TRUE, legend = TRUE,
                                             legend_position = c('bottomright','topright','topleft')) {

  p1 <- locuscomparer_make_scatterplot(df, trait1_title = trait1_title, trait2_title = trait2_title, colour = colour,
                        shape = shape, size = size, legend = legend, legend_position = legend_position)

  df1 <- df[,c('rsid', 'logp1', 'pos', 'label')]
  colnames(df1)[which(colnames(df1) == 'logp1')] = 'logp'
  p2 <- locuscomparer_make_locuszoom(df = df1, trait_title = trait1_title, colour = colour, shape = shape, size = size)

  df2 = df[,c('rsid', 'logp2', 'pos', 'label')]
  colnames(df2)[which(colnames(df2) == 'logp2')] = 'logp'
  p3 <- locuscomparer_make_locuszoom(df = df2, trait_title = trait2_title, colour = colour, shape = shape, size = size)

  if (combine) {
    p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    p4 = cowplot::plot_grid(p2, p3, align = "v", nrow = 2, rel_heights=c(0.8,1))
    p5 = cowplot::plot_grid(p1, p4)
    return(p5)
  }
  else {
    return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
  }
}

#' Create a Scatter Plot Comparing Two Traits
#'
#' This function generates a scatter plot comparing the -log10(p-values) of two different traits. The plot includes points representing SNPs, with customizable aesthetics such as color, shape, and size based on the SNPs. Additionally, it can display a legend with configurable position.
#'
#' @import ggplot2
#' @import cowplot
#'
#' @param df (data.frame) A data frame containing the following columns:
#' \itemize{
#'   \item \code{logp1} (numeric) - -log10(p-value) for the first trait.
#'   \item \code{logp2} (numeric) - -log10(p-value) for the second trait.
#'   \item \code{rsid} (character) - SNP rsID.
#'   \item \code{label} (character) - Optional label for SNPs, used for annotation.
#'   \item Additional columns for plotting aesthetics: `colour`, `shape`, `size` can be specified for different SNPs.
#' }
#' The `logp1` and `logp2` columns represent the transformed p-values for the two traits being compared.
#'
#' @param trait1_title (string) The title for the x-axis of the scatter plot representing the first trait. Default is `"Trait 1"`.
#' @param trait2_title (string) The title for the y-axis of the scatter plot representing the second trait. Default is `"Trait 2"`.
#' @param colour (named vector) A vector specifying the colors for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = "#FF0000", "rs2" = "#0000FF", "rs3" = "#00FF00")`.
#' @param shape (named vector) A vector specifying the shapes for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = 21, "rs2" = 22, "rs3" = 23)`.
#' @param size (named vector) A vector specifying the sizes for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = 2, "rs2" = 3, "rs3" = 4)`.
#' @param legend (boolean) Whether to include a legend in the plot. If \code{TRUE}, a legend will be displayed showing the colors, shapes, and sizes associated with different SNPs. Default is \code{TRUE}.
#' @param legend_position (string) Position of the legend on the plot. Choose from 'bottomright', 'topright', or 'topleft'. Default is 'bottomright'.
#'
#' @return A \code{ggplot} object representing the scatter plot comparing two traits.
#' @export
locuscomparer_make_scatterplot = function (df,
                                           trait1_title,
                                           trait2_title,
                                           colour, shape, size,
                                           legend = TRUE, legend_position = c('bottomright','topright','topleft')) {

  x_max <- max(df$logp1, na.rm = TRUE) * 1.1
  y_max <- max(df$logp2, na.rm = TRUE) * 1.1

  p <- ggplot(df, aes(logp1, logp2)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid),
               alpha = 0.8) +
    geom_point(data = df[df$label != "",],
               aes(x = logp1, y = logp2,
                   fill = rsid, size = rsid, shape = rsid)) +
    xlab(bquote(.(trait1_title) ~ -log[10] * '(P)')) +
    ylab(bquote(.(trait2_title) ~ -log[10] * '(P)')) +
    scale_fill_manual(values = colour, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    ggrepel::geom_text_repel(aes(label = label)) +
    scale_x_continuous(limits = c(0, x_max)) +
    scale_y_continuous(limits = c(0, y_max)) +
    theme_classic()

  if (legend == TRUE) {
    legend_position = match.arg(legend_position)
    if (legend_position == 'bottomright'){
      legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
    } else if (legend_position == 'topright'){
      legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
    } else {
      legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
    }

    p = cowplot::ggdraw(p) +
      geom_rect(data = legend_box,
                aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                color = "black",
                fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
      cowplot::draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 20) +
      cowplot::draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 20) +
      cowplot::draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 20) +
      cowplot::draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 20)
  p
    }

  return(p)
}

#' Create a LocusZoom Plot for SNP Associations
#'
#' This function generates a LocusZoom-style plot, which visualizes the association between SNPs and a trait along a chromosome. The plot displays SNP positions and their significance (-log10(p-values)), with customizable aesthetics for different SNPs. For more details on LocusZoom plots, see \url{http://locuszoom.org/}.
#'
#' @import ggplot2
#' @import cowplot
#'
#' @param df (data.frame) A data frame with the following columns:
#' \itemize{
#'   \item \code{pos} (numeric) - Position of the SNP on the chromosome.
#'   \item \code{logp} (numeric) - -log10(p-value) for the SNP-trait association.
#'   \item \code{rsid} (character) - SNP rsID.
#'   \item \code{label} (character) - Optional label for SNPs, used for annotation. Only SNPs with non-empty labels will be annotated.
#'   \item Additional columns for aesthetics, such as `colour`, `shape`, `size` that map to `rsid` values for visualization.
#' }
#' Example data frame:
#' \preformatted{
#' df <- data.frame(
#'   pos = c(1, 2, 3, 4, 5),
#'   logp = c(3.0, 2.5, 4.0, 1.5, 2.0),
#'   rsid = c("rs1", "rs2", "rs3", "rs4", "rs5"),
#'   label = c("", "", "Lead SNP", "", "Another SNP")
#' )
#' }
#'
#' @param ylab (string) The label for the y-axis of the plot. This represents the measure of association, typically "-log10(p-value)". Default is \code{"-log10(p)"}.
#' @param trait_title (string) The title for the plot. This represents the name of the trait or association being visualized. Default is \code{"exposure"}.
#' @param colour (named vector) A vector specifying the colors for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = "#FF0000", "rs2" = "#0000FF", "rs3" = "#00FF00")`.
#' @param shape (named vector) A vector specifying the shapes for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = 21, "rs2" = 22, "rs3" = 23)`.
#' @param size (named vector) A vector specifying the sizes for different SNPs. The names of the vector should correspond to `rsid` values.
#' Example: `c("rs1" = 2, "rs2" = 3, "rs3" = 4)`.
#'
#' @return A \code{ggplot} object representing the LocusZoom plot showing the SNP positions and their -log10(p-values) for the trait.
#' @export
locuscomparer_make_locuszoom <- function(df, ylab = "-log10(p)", trait_title = "exposure", colour, shape, size) {

  # Define minimum and maximum position for the x-axis
  x_min <- min(df$pos)
  x_max <- max(df$pos)

  # Create ggplot object
  p <- ggplot(df, aes(x = pos, y = logp)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(data = df[df$label != '', ],
               aes(x = pos, y = logp, fill = rsid, size = rsid, shape = rsid)) +
    scale_fill_manual(values = colour, guide = 'none') +
    scale_shape_manual(values = shape, guide = 'none') +
    scale_size_manual(values = size, guide = 'none') +
    scale_x_continuous(
      breaks = c(x_min, x_max),  # Set breaks to min and max values
      name = "Position",
      limits = c(x_min, x_max)  # Ensure limits include the data range
    ) +
    ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf) +
    labs(x = "Position", y = ylab, title = trait_title) +
    theme_classic() +
    theme(
      plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"),
      axis.text.x = element_text(size = 12),  # Adjust text size for x-axis labels
      axis.title.x = element_text(size = 14, face = "bold"),  # Adjust size and style for x-axis title
      axis.title.y = element_text(size = 14, face = "bold")  # Adjust size and style for y-axis title
    )

  return(p)
}

