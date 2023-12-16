#' @title proxy_search: search for proxies when using read_outcome_data()
#' @return a data frame
#' @param data_exposure exposure data frame
#' @param data_outcome outcome data frame
#' @param data_outcome_path file path used for read_outcome_data()
#' @param data_reference reference data; bim file (i.e., paste(data_reference_path, ".bim"))
#' @param data_reference_path file path for your downloaded reference panel
#' @param tag_r2 r2 for proxy SNP; from get_ld_proxies(); default = 0.8
#' @param tag_kb window to look for proxy SNPs; from get_ld_proxies(); default = 5000
#' @param tag_nsnp from get_ld_proxies(); default = 5000

proxy_search <- function(data_exposure, data_outcome, data_outcome_path, data_reference, data_reference_path,
                         tag_r2 = 0.8, tag_kb = 5000, tag_nsnp = 5000) {

  # Parameter Validation
  if (!file.exists(data_outcome_path) || !file.exists(paste0(data_reference))) {
    stop("Invalid file paths provided.")
  }

  if (!is.data.frame(data_exposure) || !is.data.frame(data_outcome)) {
    stop("data_exposure and data_outcome must be data frames.")
  }

  # exposure snps missing from outcome ====
  snps_missing <- setdiff(unique(as.factor(data_exposure$SNP)), unique(as.factor(data_outcome$SNP)))

  # look-up missing SNPs in reference panel ====
  message(paste0("# 1. looking up ", length(unique(as.factor(snps_missing))), " missing-SNP(s) in the reference panel"))
  reference_header <- readLines(data_reference, n = 1) ## read the first row to identify the column containing "rs*"
  column_index <- grep("rs", strsplit(reference_header, "\t")[[1]]) ## get column index
  ## bash command to filter rows based on the identified column
  cmd <- sprintf("awk '$%d ~ /^(%s)$/' %s", column_index, paste(snps_missing, collapse = "|"), paste0(data_reference))
  ## use system to run the bash command and then read the result with fread
  reference <- data.table::fread(cmd = cmd, quote = "")
  snps_reference <- intersect(unique(as.factor(snps_missing)), unique(as.factor(reference$V2)))
  message(paste0("## ", length(unique(as.factor(snps_reference))), " of ", length(unique(as.factor(snps_missing))), " missing-SNP(s) available in the reference panel"))

  # find proxies for available SNPs ====
  message(paste0("# 2. extracting proxy-SNP(s) for the ", length(unique(as.factor(snps_reference))), " missing-SNP(s) from the reference panel"))
  gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
  proxies <- functions::get_ld_proxies(rsid = snps_missing,
                                       bfile = data_reference_path,
                                       searchspace = NULL,
                                       tag_kb = tag_kb,
                                       tag_nsnp = tag_nsnp,
                                       tag_r2 = R2,
                                       threads = 1,
                                       out = tempfile())
  ## format proxy data: change column order and names, add proxy.outcome = TRUE
  proxies <- proxies %>%
    select(target_snp.outcome = SNP_A,
           proxy_snp.outcome = SNP_B,
           target_a1.outcome = A1,
           target_a2.outcome = A2,
           proxy_a1.outcome = B1,
           proxy_a2.outcome = B2,
           R) %>%
    mutate(proxy.outcome = TRUE,
           SNP = proxy_snp.outcome) %>%
    select(proxy.outcome, everything())
  message(paste0("## proxy-SNP(s) for ", length(unique(proxies$target_snp.outcome)), " missing-SNP(s) found; ", "proxy-SNP(s) for ", length(unique(as.factor(snps_reference))) - length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) not available (e.g., no proxy-SNP or r2 < provided)"))

  # extract proxies from outcome ====
  message(paste0("# 3. extracting proxy-SNP(s) from outcome"))
  proxy_snps <- proxies %>% # select unique proxy SNPs to extract
    distinct(proxy_snp.outcome) %>%
    pull(proxy_snp.outcome)
  data_outcome_proxies <- TwoSampleMR::read_outcome_data(filename = data_outcome_path,
                                                         snps = proxy_snps,
                                                         sep = "\t",
                                                         phenotype_col = "phenotype",
                                                         snp_col = "SNP",
                                                         beta_col = "BETA",
                                                         se_col = "SE",
                                                         eaf_col = "EAF",
                                                         effect_allele_col = "EA",
                                                         other_allele_col = "OA",
                                                         pval_col = "P",
                                                         samplesize_col = "N",
                                                         id_col = "phenotype",
                                                         chr_col = "CHR",
                                                         pos_col = "POS")
  data_outcome_proxies <- left_join(data_outcome_proxies, proxies, by = c("SNP" = "SNP"))
  message(paste0("## proxy-SNP(s) for ", length(unique(as.factor(data_outcome_proxies$target_snp.outcome))), " of ", length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) extracted"))

  # select proxy-SNP(s) with the highest R2 ====
  data_outcome_proxies <- data_outcome_proxies %>%
    group_by(target_snp.outcome) %>%
    filter(R == max(R)) %>%
    slice(1) %>%
    select(-R)

  ## Bind rows of data_outcome with data_outcome_proxies
  data_outcome <- bind_rows(data_outcome, data_outcome_proxies)

  return(data_outcome)
}
