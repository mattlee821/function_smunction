rm(list=ls())
set.seed(821)

# environment ====
remotes::install_github("chr1swallace/coloc@main",build_vignettes=TRUE)
library(coloc)
remotes::install_github("mattlee821/functions",build_vignettes=TRUE)
library(functions)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

# data ====
data(coloc_test_data)
data_coloc_exposure <- coloc_test_data$D1[c("beta","varbeta","snp","position","type","sdY", "LD", "MAF")]
data_coloc_exposure$N <- 10000
data_coloc_exposure$chr <- 1
data_coloc_exposure$pval <- runif(n = length(data_coloc_exposure$snp), min = 5e-200, max = 0.1)
str(data_coloc_exposure)

data_coloc_outcome <- coloc_test_data$D2[c("beta","varbeta","snp","position", "LD", "MAF")]
data_coloc_outcome$type <- "cc"
data_coloc_outcome$N <- 10000
data_coloc_outcome$chr <- 1
data_coloc_outcome$pval <- runif(n = length(data_coloc_outcome$snp), min = 5e-200, max = 0.1)
str(data_coloc_outcome)

# VARS ====
priors <- list(
  list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
  list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-6),
  list(p1 = 1e-4, p2 = 1e-5, p12 = 1e-6),
  list(p1 = 1e-5, p2 = 1e-4, p12 = 1e-6),
  list(p1 = 1e-5, p2 = 1e-5, p12 = 1e-7)
)
label_priors <- c(
  "p1=1e-4;p2=1e-4;p12=1e-5",
  "p1=1e-4;p2=1e-4;p12=1e-6",
  "p1=1e-4;p2=1e-5;p12=1e-6",
  "p1=1e-5;p2=1e-4;p12=1e-6",
  "p1=1e-5;p2=1e-5;p12=1e-7"
)

window <- list(
  "1mb",
  "1mb",
  "1mb",
  "1mb",
  "1mb"
)

SNP <- "s105"

# data check ====
coloc::check_dataset(d = data_coloc_exposure, suffix = 1, warn.minp=5e-8)
coloc::check_dataset(d = data_coloc_outcome, suffix = 2, warn.minp=5e-8)
## plot dataset
cowplot::plot_grid(
  coloc_plot_dataset(d = data_coloc_exposure, label = "exposure"),
  coloc_plot_dataset(d = data_coloc_outcome, label = "outcome"),
  coloc_check_alignment(D = data_coloc_exposure),
  coloc_check_alignment(D = data_coloc_outcome),
  ncol = 2)

# finemap check ====
SNP_causal_exposure <- coloc::finemap.abf(dataset = data_coloc_exposure) %>%
  filter(SNP.PP == max(SNP.PP)) %>%
  select(snp, SNP.PP)

SNP_causal_outcome <- coloc::finemap.abf(dataset = data_coloc_outcome) %>%
  filter(SNP.PP == max(SNP.PP)) %>%
  select(snp, SNP.PP)

# coloc.abf ====
coloc_results <- lapply(priors, function(params) {
  coloc::coloc.abf(
    dataset1 = data_coloc_exposure,
    dataset2 = data_coloc_outcome,
    p1 = params$p1,
    p2 = params$p2,
    p12 = params$p12
  )
})

# results table ====
table_coloc <- data.frame()
# Loop through the coloc_results list
for (j in seq_along(coloc_results)) {
  results_coloc <- coloc_results[[j]]
  results <- data.frame(
    exposure = "exposure",
    exposure_sex = "exposure_sex",
    exposure_study = "exposure_study",
    exposure_data = "exposure_data",
    exposure_population = "exposure_population",
    outcome = "outcome",
    outcome_sex = "outcome_sex",
    outcome_study = "outcome_study",
    outcome_data = "outcome_data",
    outcome_population = "outcome_population",
    SNP = SNP,
    window = window[[j]],
    finemap_snp_exposure = SNP_causal_exposure$snp,
    finemap_snp_exposure_PP = SNP_causal_exposure$SNP.PP,
    finemap_snp_outcome = SNP_causal_outcome$snp,
    finemap_snp_outcome_PP = SNP_causal_outcome$SNP.PP,
    nsnps = results_coloc$summary[["nsnps"]],
    h0 = results_coloc$summary[["PP.H0.abf"]],
    h1 = results_coloc$summary[["PP.H1.abf"]],
    h2 = results_coloc$summary[["PP.H2.abf"]],
    h3 = results_coloc$summary[["PP.H3.abf"]],
    h4 = results_coloc$summary[["PP.H4.abf"]],
    prior_p1 = results_coloc$priors[["p1"]],
    prior_p2 = results_coloc$priors[["p2"]],
    prior_p12 = results_coloc$priors[["p12"]],
    priors_label = label_priors[j]
  )
  # Bind the current results to the accumulated table
  table_coloc <- bind_rows(table_coloc, results)
}

# sensitivity  ====
plot <- coloc_sensitivity(
  obj = coloc_results[[1]],
  rule = "H4 > 0.8", npoints = 100, row = 1, suppress_messages = TRUE,
  trait1_title = "exposure", trait2_title = "outcome",
  dataset1 = NULL, dataset2 = NULL,
  data_check_trait1 = data_coloc_exposure, data_check_trait2 = data_coloc_outcome
)
plot
