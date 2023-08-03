#' @title qqplot
#' @param df A data frame AND p-value column i.e., data$P
#' @param ci A number between 0 and 1 (defaults to 0.95) indicating the type of confidence interval to be drawn. Use NA to remove.

#' @param point_size the size of the points - defaults to 2
#' @param point_colour the colour of the points - defaults to black
#' @param point_alpha the alpha/transparency of the points - defaults to 0.8
#' @param title a ggtitle placed top left of the plot
#' @export
#' @import ggplot2
#' @import cowplot

qqplot <- function(df,
                      ci = 0.95,
                      point_size = 2,
                      point_colour = "black",
                      point_alpha = 0.8,
                      title = NA
                      ) {
  n  <- length(df)
  df <- data.frame(
    observed = -log10(sort(df)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), alpha = point_alpha, size = point_size, colour = point_colour) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    ggtitle(title) +
    xlab(log10Pe) +
    ylab(log10Po) +

    # Customise the theme:
    theme_cowplot() +
    theme(
      axis.text.x = element_text(hjust = 0.5,
                                 vjust = 1),
      axis.text.y = element_text(hjust = 1,
                                 vjust = 0.5)
    )
}
