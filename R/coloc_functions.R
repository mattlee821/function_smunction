##' Plots sensitivity analysis for coloc results
##'
##' @param obj Output from `coloc.detail` or `coloc.process`.
##' @param rule Decision rule for acceptable posterior probabilities.
##' @param trait1_title Title for the first trait.
##' @param trait2_title Title for the second trait.
##' @param dataset1 Optional dataset for the first trait.
##' @param dataset2 Optional dataset for the second trait.
##' @param npoints Number of points for evaluating prior values. Default is `100`.
##' @param doplot Whether to draw the plot. Default is `TRUE`.
##' @param plot.manhattans Whether to show Manhattan plots. Default is `TRUE`.
##' @param preserve.par Whether to preserve the current graphics device parameters. Default is `FALSE`.
##' @param row Row number for coloc summary if multiple rows are present. Default is `1`.
##' @return List of prior matrix, posterior matrix, and a pass/fail indicator (returned invisibly).
##' @export
##'
##' @import data.table
##' @import dplyr
my_sensitivity <- function(obj, rule = "", trait1_title, trait2_title,
                           dataset1 = NULL, dataset2 = NULL,
                           npoints = 100, doplot = TRUE, plot.manhattans = TRUE,
                           preserve.par = FALSE,
                           row = 1) {
  stopifnot("list" %in% class(obj))
  stopifnot("priors" %in% names(obj))
  stopifnot("summary" %in% names(obj))
  if (rule == "")
    stop("please supply a rule to define colocalisation, eg 'H4 > thr' where thr is some probability of H4 that you accept as colocalisation")
  rule.init <- rule
  rule <- gsub("(H.)", "PP.\\1.abf", rule, perl = TRUE)

  ## massage results object
  results = obj$results
  ## multiple signals?
  multiple = FALSE
  if (data.table::is.data.table(obj$summary)) { # we're not in coloc.abf anymore
    if (!(row %in% 1:nrow(obj$summary)))
      stop("row must be between 1 and ", nrow(obj$summary))
    pp <- unlist(c(obj$summary[row, grep("PP|nsnp", names(obj$summary)), with = FALSE]))
    if (paste0("SNP.PP.H4.row", row) %in% names(results)) {
      multiple = TRUE
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
    df1 <- with(dataset1, data.table::data.table(snp = snp, position = position, z.df1 = beta / sqrt(varbeta)))
    df2 <- with(dataset2, data.table::data.table(snp = snp, position = position, z.df2 = beta / sqrt(varbeta)))
    df <- merge(df1, df2, by = c("snp", "position"), all = TRUE)
    results <- merge(results, df, by = "snp")
  }

  p12 <- obj$priors["p12"]
  p1 <- obj$priors["p1"]
  p2 <- obj$priors["p2"]
  check <- function(pp) { with(as.list(pp), eval(parse(text = rule))) }
  pass.init <- check(pp)
  message("Results ", if (check(pp)) { "pass" } else { "fail" }, " decision rule ", rule.init)

  testp12 <- 10^seq(log10(p1 * p2), log10(min(p1, p1)), length.out = npoints)
  testH <- prior.snp2hyp(pp["nsnps"], p12 = testp12, p1 = p1, p2 = p2)
  testpp <- as.data.frame(prior.adjust(summ = pp, newp12 = testp12, p1 = p1, p2 = p2, p12 = p12))
  colnames(testpp) <- gsub("(H.)", "PP.\\1.abf", colnames(testpp), perl = TRUE)
  pass <- check(testpp)
  w <- which(pass)

  if (doplot) {
    H <- as.character(0:4)
    op <- par('mfcol', 'mar', 'mgp')
    if (!preserve.par) {
      on.exit(par(op))
    }
    if (plot.manhattans) {
      layout(matrix(1:3, ncol = 1), heights = c(2, 1, 1))
      par(mar = c(0.1, 4.1, 4.1, 2.1))
    } else {
      layout(matrix(1:2, ncol = 1), heights = c(2, 1))
    }
    matplot(log10(testp12), testpp[, paste0("PP.H", H, ".abf")], type = "l",
            col = c("black", "red", "green", "blue", "purple"),
            ylab = "posterior", xlab = "", lwd = 2)
    mtext(expression(p[12]), 1, line = 2)
    legend("topright", lwd = 2, col = c("black", "red", "green", "blue", "purple"), legend = paste0("H", H))
    abline(v = log10(p12), col = "grey")
    title(main = paste(trait1_title, " and ", trait2_title, sep = ""))
    u <- par('usr')
    rect(u[1], u[3], log10(testp12)[w[1]], u[4], col = "grey", border = NA)
    rect(log10(testp12)[w[length(w)]], u[3], u[2], u[4], col = "grey", border = NA)
    if (plot.manhattans) {
      par(mar = c(2.1, 4.1, 0.1, 2.1))
      manh.plot(results, wh = 1)
      title(xlab = "Chromosome position")
      par(mar = c(4.1, 4.1, 0.1, 2.1))
      manh.plot(results, wh = 2)
      title(xlab = "Chromosome position")
    }
  }
  invisible(list(prior = testH, posterior = testpp, pass = pass))
}


##' Performs colocalization analysis with varying prior probabilities and saves the results
##'
##' @param coloc_data_exposure Colocalization data for the exposure trait.
##' @param coloc_data_outcome Colocalization data for the outcome trait.
##' @param label Label for the analysis.
##' @param outcome Outcome trait.
##' @param FILE_NAME Name of the file to save results.
##' @param window Window size for the analysis.
perform_colocalization <- function(coloc_data_exposure, coloc_data_outcome, label, outcome, FILE_NAME, window) {
  if (is.character(coloc_data_exposure) && file.exists(coloc_data_exposure)) {
    coloc_data_exposure <- readRDS(coloc_data_exposure)
  }
  if (is.character(coloc_data_outcome) && file.exists(coloc_data_outcome)) {
    coloc_data_outcome <- readRDS(coloc_data_outcome)
  }

  new.df <- data.frame()
  counter <- 1
  for (p12 in c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) {
    try({
      my.res <- coloc::coloc.abf(dataset1 = coloc_data_exposure, dataset2 = coloc_data_outcome, p1 = 1e-4, p2 = 1e-4, p12 = p12)
      my.res$priors["p12"] <- p12
      my.summary <- coloc.process(my.res)
      sensitivity.plot <- my_sensitivity(my.res, rule = "H4 > 0.5", trait1_title = label, trait2_title = outcome, plot.manhattans = FALSE, doplot = TRUE)

      if (!dir.exists("output")) {
        dir.create("output")
      }

      write.csv(my.summary, paste0("output/", FILE_NAME, "_p12_", p12, "_summary.csv"), row.names = FALSE)
      pdf(paste0("output/", FILE_NAME, "_p12_", p12, "_sensitivity.pdf"))
      print(sensitivity.plot)
      dev.off()

      if (counter == 1) {
        new.df <- my.summary
      } else {
        new.df <- rbind(new.df, my.summary)
      }

      counter <- counter + 1
    })
  }

  write.csv(new.df, paste0("output/", FILE_NAME, "_all_summaries.csv"), row.names = FALSE)

  invisible(new.df)
}
