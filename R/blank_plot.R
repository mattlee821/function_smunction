#' Create a blank ggplot2 plot
#'
#' This function creates a blank ggplot2 plot with no data or aesthetics
#'
#' @return A ggplot2 plot object
#' @param message a message that no plot was made
#' @export
blank_plot <- function(message)
{
  requireNamespace("ggplot2", quietly=TRUE)
  ggplot2::ggplot(data.frame(a=0,b=0,n=message)) +
    ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) +
    ggplot2::labs(x=NULL,y=NULL) +
    ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}
