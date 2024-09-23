#' Plots sensitivity analysis for coloc results
#'
#' @param obj Output from `coloc.detail` or `coloc.process`.
#' @param rule Decision rule for acceptable posterior probabilities.
#' @param trait1_title Title for the first trait (default: "trait 1").
#' @param trait2_title Title for the second trait (default: "trait 2").
#' @param dataset1 Optional dataset for the first trait.
#' @param dataset2 Optional dataset for the second trait.
#' @param npoints Number of points for evaluating prior values. Default is `100`.
#' @param row Row number for coloc summary if multiple rows are present. Default is `1`.
#' @param data_check_trait1 Data for trait 1 to be checked.
#' @param data_check_trait2 Data for trait 2 to be checked.
#' @param suppress_messages Logical, whether to suppress messages (default is FALSE).
#' @return List of prior matrix, posterior matrix, and a pass/fail indicator (returned invisibly).
#' @export
coloc_sensitivity <- function(obj, rule = "H4 > 0.8", trait1_title = "trait 1", trait2_title = "trait 2",
                              dataset1 = NULL, dataset2 = NULL,
                              npoints = 100, suppress_messages = FALSE,
                              row = 1, data_check_trait1 = NULL, data_check_trait2 = NULL) {
  # check ====
  stopifnot("list" %in% class(obj))
  stopifnot("priors" %in% names(obj))
  stopifnot("summary" %in% names(obj))
  if (rule == "")
    stop("please supply a rule to define colocalisation, eg 'H4 > thr' where thr is some probability of H4 that you accept as colocalisation")
  rule.init <- rule
  rule <- gsub("(H.)", "PP.\\1.abf", rule, perl = TRUE)

  ## massage results object
  results <- obj$results
  ## multiple signals?
  multiple <- FALSE
  if (data.table::is.data.table(obj$summary)) { # we're not in coloc.abf anymore
    if (!(row %in% 1:nrow(obj$summary)))
      stop("row must be between 1 and ", nrow(obj$summary))
    pp <- unlist(c(obj$summary[row, grep("PP|nsnp", names(obj$summary)), with = FALSE]))
    if (paste0("SNP.PP.H4.row", row) %in% names(results)) {
      multiple <- TRUE
      results[["SNP.PP.H4"]] <- results[[paste0("SNP.PP.H4.row", row)]]
    }
    if (paste0("z.df1.row", row) %in% names(results)) { # might be passed here or in separate dataset objects
      results[["z.df1"]] <- results[[paste0("z.df1.row", row)]]
      results[["z.df2"]] <- results[[paste0("z.df2.row", row)]]
    } else {
      pp <- unlist(c(obj$summary[row, grep("PP|nsnp", names(obj$summary)), with = FALSE]))
    }
  } else {
    pp <- obj$summary
  }
  ## need to add z score from datasets?
  if (!is.null(dataset1) && !is.null(dataset2)) {
    df1 <- with(dataset1, data.table(snp = snp, position = position, z.df1 = beta / sqrt(varbeta)))
    df2 <- with(dataset2, data.table(snp = snp, position = position, z.df2 = beta / sqrt(varbeta)))
    df <- merge(df1, df2, by = c("snp", "position"), all = TRUE)
    results <- merge(results, df, by = "snp")
  }

  p12 <- obj$priors["p12"]
  p1 <- obj$priors["p1"]
  p2 <- obj$priors["p2"]
  check <- function(pp) { with(as.list(pp), eval(parse(text = rule))) }
  pass.init <- check(pp)
  if (!suppress_messages) {
    message("Results ", if (check(pp)) { "pass" } else { "fail" }, " decision rule ", rule.init)
  }

  testp12 <- 10^seq(log10(p1 * p2), log10(min(p1, p1)), length.out = npoints)
  testH <- coloc_prior.snp2hyp(pp["nsnps"], p12 = testp12, p1 = p1, p2 = p2)
  testpp <- as.data.frame(coloc_prior.adjust(summ = pp, newp12 = testp12, p1 = p1, p2 = p2, p12 = p12))
  colnames(testpp) <- gsub("(H.)", "PP.\\1.abf", colnames(testpp), perl = TRUE)
  pass <- check(testpp)
  w <- which(pass)

  # plots ====
  # Data check plot ====
  # plot_distribution <- cowplot::plot_grid(
  #   coloc_plot_dataset(d = data_check_trait1, label = trait1_title),
  #   coloc_plot_dataset(d = data_check_trait2, label = trait2_title),
  #   ncol = 1, nrow = 2
  # )
  plot_alignment <- cowplot::plot_grid(
    coloc_check_alignment(D = data_check_trait1),
    coloc_check_alignment(D = data_check_trait2),
    ncol = 1, nrow = 2
  )

  # locuszoom plots ====
  plot_locuscompare <- locuscomparer_make_plot(data_coloc_exposure = data_check_trait1,
                                               data_coloc_outcome = data_check_trait2,
                                               SNP_causal_exposure = SNP_causal_exposure$snp,
                                               trait1_title = trait1_title, trait2_title = trait2_title)

  # Manhattan plots ====
  if ("z.df1" %in% colnames(results) && "z.df2" %in% colnames(results)) {
    plot_manhattan <- cowplot::plot_grid(
      coloc_manh.plot(results, 1),
      coloc_manh.plot(results, 2),
      ncol = 1, nrow = 2)
  } else {
    warning("Please supply dataset1 and dataset2 to view Manhattan plots.")
  }

  # Probability plots ====
  m <- list(testH, as.matrix(testpp))
  ti <- list("Prior probabilities", "Posterior probabilities")
  plot_probabilities <- list()
  for (i in 1:2) {
    # Set max y axis value
    ym <- if (i == 1) max(m[[i]][, -1]) else max(m[[i]])

    # Data
    df <- data.frame(
      testp12 = testp12,
      m = m[[i]],
      pass = pass
    ) %>%
      tidyr::pivot_longer(
        cols = if (i == 1) c("m.H0", "m.H1", "m.H2", "m.H3", "m.H4")
        else c("m.PP.H0.abf", "m.PP.H1.abf", "m.PP.H2.abf", "m.PP.H3.abf", "m.PP.H4.abf"),
        names_to = "name", values_to = "value"
      ) %>%
      dplyr::mutate(name = forcats::fct_recode(name,
                               H0 = if (i == 1) "m.H0" else "m.PP.H0.abf",
                               H1 = if (i == 1) "m.H1" else "m.PP.H1.abf",
                               H2 = if (i == 1) "m.H2" else "m.PP.H2.abf",
                               H3 = if (i == 1) "m.H3" else "m.PP.H3.abf",
                               H4 = if (i == 1) "m.H4" else "m.PP.H4.abf"))

    # Check data
    df <- df %>%
      dplyr::filter(!is.na(testp12) & !is.na(value)) %>%
      dplyr::filter(testp12 >= 1e-8 & testp12 <= 1e-4)
    p12 <- ifelse(p12 >= 1e-8 & p12 <= 1e-4, p12, NA)

    # Initialize the plot
    plot_probabilities[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x = testp12, y = value, colour = name)) +
      ggplot2::scale_x_log10(breaks = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4)) +
      ggplot2::scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      ggplot2::labs(x = "p12", y = "Prob",
           title = ti[[i]],
           subtitle = paste("shaded region:", rule.init)) +
      ggplot2::geom_vline(xintercept = p12, linetype = "dashed", color = "gray") +
      ggplot2::annotate("text", x = p12, y = 0.5, label = "results", angle = 90, color = "gray40") +
      cowplot::theme_cowplot() +
      ggplot2::theme(legend.position = "right", legend.title = ggplot2::element_blank())

    # Add rectangle if pass condition is met
    if (any(df$pass)) {
      # Calculate x-axis limits for the rectangle
      xleft <- min(df$testp12[df$pass])
      xright <- max(df$testp12[df$pass])

      plot_probabilities[[i]] <- plot_probabilities[[i]] +
        ggplot2::geom_rect(ggplot2::aes(xmin = xleft, xmax = xright, ymin = 0, ymax = 1),
                  fill = "#E2FFBF", colour = NA, alpha = 0.005)
    }

    # Add the points on top of the rectangle
    plot_probabilities[[i]] <- plot_probabilities[[i]] +
      ggplot2::geom_point()
  }
  # combined
  plot_probabilities <- cowplot::plot_grid(plot_probabilities[[1]],
                                           plot_probabilities[[2]],
                                           ncol = 1, nrow = 2)

  # combined plot ====
  plot_row1 <- cowplot::plot_grid(
    # plot_distribution,  # A
    plot_alignment,     # B
    plot_manhattan,     # C
    labels = c("A", "B"),
    ncol = 2,
    rel_widths = c(1, 1)
  )
  plot_row2 <- cowplot::plot_grid(
    plot_locuscompare,  # D
    plot_probabilities, # E
    labels = c("C", "D"),
    ncol = 2,
    rel_widths = c(1.3, 0.7)
  )

  # Combine the two rows into one plot
  plot_sensitivity <- cowplot::plot_grid(
    plot_row1,
    plot_row2,
    ncol = 1,
    rel_heights = c(1, 1)
  )

  plot_sensitivity

  return(plot_sensitivity)
  invisible(cbind(testpp, p12 = testp12, pass = pass))
}

#' Convert SNP Prior to Hypotheses
#'
#' Converts SNP prior probabilities to hypotheses probabilities for colocalisation.
#'
#' @param nsnp Number of SNPs.
#' @param p12 Probability p12.
#' @param p1 Probability p1.
#' @param p2 Probability p2.
#'
#' @return Matrix of hypotheses probabilities.
#' @export
coloc_prior.snp2hyp <- function(nsnp, p12 = 1e-6, p1 = 1e-4, p2 = 1e-4) {
  if (any(p12 < p1 * p2) || any(p12 > p1) || any(p12 > p2))
    return(NULL)
  tmp <- cbind(nsnp * p1,
               nsnp * p2,
               nsnp * (nsnp - 1) * p1 * p2,
               nsnp * p12)
  tmp <- cbind(1 - rowSums(tmp), tmp)
  ## if(nrow(tmp)==1) {
  ##     tmp <- c(tmp)
  ##     names(tmp) <- paste0("H",0:4)
  ## } else
  colnames(tmp) <- paste0("H", 0:4)
  tmp
}

#' Adjust Colocalisation Prior Probabilities
#'
#' Adjusts colocalisation prior probabilities based on new p12 values.
#'
#' @param summ Summary vector or list containing necessary data.
#' @param newp12 New p12 values for adjustment.
#' @param p1 Probability p1.
#' @param p2 Probability p2.
#' @param p12 Probability p12.
#'
#' @return Adjusted posterior probabilities matrix.
#' @export
coloc_prior.adjust <- function(summ, newp12, p1 = 1e-4, p2 = 1e-4, p12 = 1e-6) {
  if (is.list(summ) && "summary" %in% names(summ))
    summ <- summ$summary
  if (!identical(names(summ), c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")))
    stop("not a coloc summary vector")
  ## back calculate likelihoods
  f <- function(p12)
    coloc_prior.snp2hyp(summ["nsnps"], p12 = p12, p1 = p1, p2 = p2)
  pr1 <- f(newp12)
  pr0 <- matrix(f(p12), nrow = nrow(pr1), ncol = ncol(pr1), byrow = TRUE)
  ## make everything conformable
  ## if(is.matrix(summ) && nrow(summ)==1) summ <- as.vector(summ)
  ## if(nrow(pr1)==1) pr1 <- as.vector(pr1)
  ## if(nrow(pr1)>1) pr1 <- t(pr1)
  newpp <- matrix(summ[-1], nrow = nrow(pr1), ncol = ncol(pr1), byrow = TRUE) * pr1 / pr0 # prop to, not equal to
  newpp / rowSums(newpp)
}

#' Plot Manhattan Plot for Colocalisation
#'
#' Generates a Manhattan plot based on colocalisation data.
#'
#' @param df Data frame containing colocalisation data.
#' @param wh Which trait to plot (1 or 2).
#' @param position Chromosome positions or SNP numbers.
#'
#' @export
coloc_manh.plot <- function(df, wh,
                            position = if ("position" %in% names(df)) {
                              df$position
                            } else {
                              1:nrow(df)
                            }) {
  znm <- if (wh == 1) { "z.df1" } else { "z.df2" }

  # Calculate -log10(p) values
  logp <- - (stats::pnorm(-abs(df[[znm]]), log.p = TRUE) + log(2)) / log(10)

  # Color palette for plotting
  Pal <- grDevices::colorRampPalette(c('white', 'blue'))
  Col <- Pal(100)[ceiling(100 * df$SNP.PP.H4)]

  # Create ggplot object
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = position, y = logp, fill = SNP.PP.H4)) +
    ggplot2::geom_point(shape = 21, size = 3, color = "grey", ggplot2::aes(fill = SNP.PP.H4)) +
    ggplot2::scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1), breaks = seq(0, 1)) +
    ggplot2::labs(x = if ("position" %in% names(df)) "position" else "SNP number",
         y = "-log10(p)",
         title = paste("Distribution / H4"),
         fill = "H4") +
    ggplot2::scale_x_continuous(
      breaks = c(min(df$position), max(df$position)),  # Set breaks to min and max values
      labels = c(min(df$position), max(df$position))  # Set labels to min and max values
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "right",
          legend.title = ggplot2::element_text(hjust = 0, vjust = 1))  # Center-align the legend title
  return(p)
}

#' Check Alignment for Colocalisation
#'
#' This function checks the alignment for colocalisation based on the product of Z scores and LD.
#'
#' @param D Dataset containing beta, varbeta, and LD.
#' @param thr Threshold for alignment check.
#' @param do_plot Logical indicating whether to generate a plot.
#'
#' @return Mean of values greater than zero, indicating alignment quality.
#'
#' @export
coloc_check_alignment <- function(D, thr = 0.2, do_plot = TRUE) {
  coloc::check_dataset(D)
  bprod <- outer(D$beta / sqrt(D$varbeta), D$beta / sqrt(D$varbeta), "*")
  tmp <- (bprod / D$LD)[abs(D$LD) > 0.2]

  if (do_plot) {
    # Calculate % positive values
    percent_positive <- 100 * mean(tmp > 0)

    # Create histogram with ggplot2
    p <- ggplot2::ggplot(data = data.frame(tmp = tmp), ggplot2::aes(x = tmp)) +
      ggplot2::geom_histogram(bins = 30, color = "black") +
      ggplot2::labs(x = "Ratio: product of Z scores to LD",
           caption = "most values should be positive;\n symmetry = potentially poor alignment",
           y = "Count",
           title = paste0("Alignment: ", round(percent_positive,2), "% positive")) +
      ggplot2::geom_vline(xintercept = 0, color = "red") +
      cowplot::theme_cowplot()
    return(p)
  }

  return(mean(tmp > 0))
}

#' @title plot a coloc dataset
#' @param d a coloc dataset
#' @param label title
#' @param susie_obj optional, the output of a call to runsusie()
#' @param highlight_list optional, a list of character vectors. any snp in the
#'   character vector will be highlighted, using a different colour for each
#'   list.
#' @param alty default is to plot a standard manhattan. If you wish to plot a
#'   different y value, pass it here. You may also want to change ylab to
#'   describe what you are plotting.
#' @param ylab label for y axis, default is -log10(p) and assumes you are
#'   plotting a manhattan
#' @param show_legend optional, show the legend or not. default is TRUE
#' @param color optional, specify the colours to use for each credible set when
#'   susie_obj is supplied. Default is shamelessly copied from
#'   susieR::susie_plot() so that colours will match
#' @param ... other arguments passed to the base graphics plot() function
#' @author Chris Wallace
#' @export
coloc_plot_dataset <- function(d, label = "distribution", susie_obj = NULL, highlight_list = NULL, alty = NULL,
                               ylab = "-log10(p)", show_legend = TRUE,
                               color = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                                         "gold1", "skyblue2", "#FB9A99", "palegreen2",
                                         "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
                                         "maroon", "orchid1", "deeppink1", "blue1",
                                         "steelblue4", "darkturquoise", "green1",
                                         "yellow4", "yellow3", "darkorange4", "brown"), ...) {

  # Convert list to data frame
  df <- data.frame(
    position = d$position,
    beta = d$beta,
    varbeta = d$varbeta,
    MAF = d$MAF,
    type = d$type,
    N = d$N,
    snp = d$snp,
    LD = as.vector(d$LD)  # Flatten LD matrix to vector
  )

  # Calculate y values
  if (is.null(alty)) {
    df$y <- -(stats::pnorm(-abs(df$beta) / sqrt(df$varbeta), log.p = TRUE) + log(2)) / log(10)
  } else {
    df$y <- alty
  }

  # Create ggplot object
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = position, y = y)) +
    ggplot2::geom_point(shape = 21, size = 3, color = "white", fill = "grey") +
    ggplot2::labs(x = "Position", y = ylab, title = label) +
    ggplot2::scale_x_continuous(
      breaks = c(min(df$position), max(df$position)),  # Set breaks to min and max values
      labels = c(min(df$position), max(df$position))  # Set labels to min and max values
    ) +
    cowplot::theme_cowplot()
  # Add highlight points if highlight_list is provided
  if (!is.null(highlight_list)) {
    p <- p +
      ggplot2::geom_point(data = df[df$snp %in% unlist(highlight_list), ], ggplot2::aes(color = factor(snp))) +
      ggplot2::scale_color_manual(values = color[1:length(highlight_list)]) +
      ggplot2::guides(color = ggplot2::guide_legend(title = "Highlight", keywidth = 1, keyheight = 1))
  }
  return(p)
}

