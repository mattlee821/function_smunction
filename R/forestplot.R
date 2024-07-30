# modification to ggforestplot forestplot() function to increase the useable space for plotting points

#' @title forestplot function
#' @param df A data frame with the data to plot. It must contain at least three
#' variables, a character column with the names to be displayed on the y-axis
#' (see parameter \code{name}), a numeric column with the value (or the log of the
#' value) to display (see parameter \code{estimate}) and a numeric value with
#' the corresponding standard errors (see parameter \code{se}). It may contain
#' additional columns, e.g. the corresponding p-values (see parameter \code{pvalue})
#' in which case, in conjuction with the threshold given in \code{psignif}, the
#' non-significant results will be displayed as hollow points. Other variables
#' may be used as aesthetics to define the colour and the shape of the points
#' to be plotted.
#' @param name the variable in \code{df} that contains the y-axis
#' names. This argument is automatically \link[rlang:quotation]{quoted} and
#' \link[rlang:eval_tidy]{evaluated} in the context of the \code{df} data frame.
#' See Note.
#' @param estimate the variable in \code{df} that contains the values (or log of
#' values) to be displayed. This argument is automatically
#' \link[rlang:quotation]{quoted} and \link[rlang:eval_tidy]{evaluated} in the
#' context of the \code{df} data frame.
#' See Note.
#' @param se the variable in the \code{df} data frame that contains the standard
#' error values. This argument is automatically \link[rlang:quotation]{quoted}
#' and \link[rlang:eval_tidy]{evaluated} in the context of the \code{df} data
#' frame. See Note.
#' @param pvalue the variable in \code{df} that contains the
#' p-values. Defaults to NULL. When explicitly defined, in conjuction with
#' the p-value threshold provided in the \code{psignif}, the non-significant
#' entries will be drawn as hollow points. This argument is automatically
#' \link[rlang:quotation]{quoted} and \link[rlang:eval_tidy]{evaluated} in the
#' context of the \code{df} data frame. See Note.
#' @param colour the variable in \code{df} by which to colour the different
#' groups of points. This argument is automatically \link[rlang:quotation]{quoted} and
#' \link[rlang:eval_tidy]{evaluated} in the context of the \code{df} data frame.
#' See Note.
#' @param shape the variable in \code{df} by which to shape the different groups of
#' points. This argument is automatically \link[rlang:quotation]{quoted} and
#' \link[rlang:eval_tidy]{evaluated} in the context of the \code{df} data frame.
#' See Note.
#' @param logodds logical (defaults to FALSE) specifying whether the \code{estimate}
#' parameter should be treated as log odds/hazards ratio (TRUE) or not (FALSE). When
#' \code{logodds} = TRUE the estimates and corresponding confidence intervals will be
#' exponentiated and a log scale will be used for the x-axis.
#' @param psignif numeric, defaults to 0.05. The p-value threshold
#' for statistical significance. Entries with larger than \code{psignif} will be
#' drawn with a hollow point.
#' @param ci A number between 0 and 1 (defaults to 0.95) indicating the type of
#' confidence interval to be drawn.
#' @param ... \code{ggplot2} graphical parameters such as \code{title},
#' \code{ylab}, \code{xlab}, \code{xtickbreaks} etc. to be passed along.
#' @return A \code{ggplot} object.
#' @note  See \code{vignette(programming, package = "dplyr")} for an
#' introduction to non-standard evaluation.
#' @author Maria Kalimeri, Ilari Scheinin, Vilma Jagerroos
#' @export
#' @importFrom stats qnorm
#' @importFrom rlang := !! enquo quo_is_null
#' @importFrom grDevices axisTicks
#' @importFrom scales trans_breaks
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_shape_manual labs coord_cartesian ylim geom_vline guide_legend
#' @importFrom ggforestplot theme_forest scale_colour_ng_d scale_fill_ng_d geom_stripes geom_effect
#' @importFrom ggstance position_dodgev

forestplot <- function (df, name = name, estimate = estimate, se = se, pvalue = NULL,
                           colour = NULL, shape = NULL, logodds = FALSE, psignif = 0.05,
                           ci = 0.95, ...)
{
  stopifnot(is.data.frame(df))
  stopifnot(is.logical(logodds))
  name <- enquo(name)
  estimate <- enquo(estimate)
  se <- enquo(se)
  pvalue <- enquo(pvalue)
  colour <- enquo(colour)
  shape <- enquo(shape)
  args <- list(...)
  const <- stats::qnorm(1 - (1 - ci)/2)
  df <- df %>% dplyr::mutate(`:=`(!!name, factor(!!name, levels = !!name %>%
                                                   unique() %>% rev(), ordered = TRUE)), .xmin = !!estimate -
                               const * !!se, .xmax = !!estimate + const * !!se, .filled = TRUE,
                             .label = sprintf("%.2f", !!estimate))
  if (logodds) {
    df <- df %>% dplyr::mutate(.xmin = exp(.data$.xmin), .xmax = exp(.data$.xmax),
                        `:=`(!!estimate, exp(!!estimate)))
  }
  if (!quo_is_null(pvalue)) {
    df <- df %>% dplyr::mutate(.filled = !!pvalue < !!psignif)
  }
  g <- ggplot2::ggplot(df, aes(x = !!estimate, y = !!name))
  if (logodds) {
    if ("xtickbreaks" %in% names(args)) {
      g <- g + ggplot2::scale_x_continuous(trans = "log10", breaks = args$xtickbreaks)
    }
    else {
      g <- g + ggplot2::scale_x_continuous(trans = "log10", breaks = scales::log_breaks(n = 7))
    }
  }
  g <- g +
    ggforestplot::theme_forest() +
    ggforestplot::scale_colour_ng_d() +
    ggforestplot::scale_fill_ng_d() +
    ggforestplot::geom_stripes() +
    ggplot2::geom_vline(xintercept = ifelse(test = logodds, yes = 1, no = 0), linetype = "solid", size = 0.4, colour = "black")
  g <- g +
    ggforestplot::geom_effect(aes(xmin = .data$.xmin, xmax = .data$.xmax, colour = !!colour, shape = !!shape, filled = .data$.filled),
                       position = ggstance::position_dodgev(height = 0.9)) +
    ggplot2::scale_shape_manual(values = c(21L, 22L, 23L, 24L, 25L)) +
    ggplot2::guides(colour = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE))
  if ("title" %in% names(args)) {
    g <- g + ggplot2::labs(title = args$title)
  }
  if ("subtitle" %in% names(args)) {
    g <- g + ggplot2::labs(subtitle = args$subtitle)
  }
  if ("caption" %in% names(args)) {
    g <- g + ggplot2::labs(caption = args$caption)
  }
  if ("xlab" %in% names(args)) {
    g <- g + ggplot2::labs(x = args$xlab)
  }
  if (!"ylab" %in% names(args)) {
    args$ylab <- ""
  }
  g <- g + ggplot2::labs(y = args$ylab)
  if ("xlim" %in% names(args)) {
    g <- g + ggplot2::coord_cartesian(xlim = args$xlim)
  }
  if ("ylim" %in% names(args)) {
    g <- g + ggplot2::ylim(args$ylim)
  }
  if ("xtickbreaks" %in% names(args) & !logodds) {
    g <- g + ggplot2::scale_x_continuous(breaks = args$xtickbreaks)
  }
  g
}
