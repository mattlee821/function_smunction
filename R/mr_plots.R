#' Funnel plot
#'
#' Create funnel plot from single SNP analyses.
#'
#' @param singlesnp_results from mr_singlesnp().
#'
#' @export
#' @return List of plots
mr_funnel_plot <- function(singlesnp_results)
{
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        functions::blank_plot("Insufficient number of SNPs")
      )
    }
    am <- grep("All", d$SNP, value=TRUE)
    d$SNP <- gsub("All - ", "", d$SNP)
    am <- gsub("All - ", "", am)
    ggplot2::ggplot(subset(d, ! SNP %in% am), ggplot2::aes(y = 1/se, x=b)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(data=subset(d, SNP %in% am), ggplot2::aes(xintercept=b, colour = SNP)) +
      # ggplot2::scale_colour_brewer(type="qual") +
      ggplot2::scale_colour_manual(values = c("#a6cee3",
                                              "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                                              "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                                              "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(y=expression(1/SE[IV]), x=expression(beta[IV]), colour="MR Method") +
      ggplot2::ggtitle(paste0("exposure = ", d$id.exposure[1], "; outcome = ", d$id.outcome[1])) +
      ggplot2::theme(legend.position="top", legend.direction="vertical")
  })
  res
}

#' modification of TwoSampleMR mr_scatter_plot() function that adds exposure/outcome labels
#'
#' Requires dev version of ggplot2
#'
#' @param mr_results Output from mr().
#' @param dat Output from harmonise_data().
#' @export
#' @return List of plots
#'
mr_scatter_plot <- function (mr_results, dat)
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"),
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, y = beta.outcome)) +
                           ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), colour = "grey", width = 0) +
                           ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
                           ggplot2::geom_point(ggplot2::aes(text = paste("SNP:", SNP))) +
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method), show.legend = TRUE) +
                           ggplot2::labs(colour = "MR Test", x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on", d$outcome[1]), title = paste(d$id.exposure[1], "on", d$outcome[1])) +
                           ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal") +
                           ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))
                       })
  mrres
}


