#' Create a blank ggplot2 plot
#'
#' This function creates a blank ggplot2 plot with no data or aesthetics. It can be used as a placeholder
#' or to customize axes and labels without displaying any data points.
#'
#' @return A ggplot2 plot object
#' @param message Optional. A character string to be displayed as a message within the plot. Defaults to NULL (no message).
#' @export
blank_plot <- function(message)
{
  ggplot2::ggplot(data.frame(a=0,b=0,n=message)) +
    ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) +
    ggplot2::labs(x=NULL,y=NULL) +
    ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}
