#' @title proxy_search: search for proxies when using read_outcome_data()
#' @return a data frame
#' @param data_exposure exposure data frame
#' @param data_outcome outcome data frame
#' @param data_outcome_path file path used for read_outcome_data()
#' @param data_reference_path file path for your downloaded reference panel

proxy_search <- function(data_exposure, data_outcome, data_outcome_path, data_reference_path) {

  # exposure snps missing from outcome ====
  snps_missing <- setdiff(unique(as.factor(data_exposure$SNP)), unique(as.factor(data_outcome$SNP)))

  # look-up missing SNPs in reference panel ====
  print(paste0("# 1. looking up ", length(unique(as.factor(snps_missing))), " missing-SNP(s) in reference panel"))
  reference <- data.table::fread(paste0(data_reference_path,".bim"))[V2 %in% snps_missing]
  snps_reference <- intersect(unique(as.factor(snps_missing)), unique(as.factor(reference$V2)))
  print(paste0("## ", length(unique(as.factor(snps_reference))), " of ", length(unique(as.factor(snps_missing))), " missing-SNP(s) available in reference panel"))

  # find proxies for available SNPs ====
  print(paste0("# 2. extracting proxy-SNP(s) for the ", length(unique(as.factor(snps_reference))), " missing-SNP(s) from the reference panel"))
  gwasvcf::set_plink(genetics.binaRies::get_plink_binary())
  proxies <- gwasvcf::get_ld_proxies(rsid = snps_missing,
                                     bfile = data_reference_path,
                                     searchspace = NULL,
                                     tag_kb = 5000,
                                     tag_nsnp = 5000,
                                     tag_r2 = 0.8,
                                     threads = 1,
                                     out = tempfile())
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
  print(paste0("## proxy-SNP(s) for ", length(unique(proxies$target_snp.outcome)), " missing-SNP(s) found; ", "proxy-SNP(s) for ", length(unique(as.factor(snps_reference))) - length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) not available (e.g., no proxy-SNP or r2 < provided)"))

  # extract proxies from outcome ====
  print(paste0("# 3. extracting proxy-SNP(s) from outcome"))
  data_outcome_proxies <- TwoSampleMR::read_outcome_data(filename = data_outcome_path,
                                            snps = proxies$proxy_snp.outcome,
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
  print(paste0("## proxy-SNP(s) for ", length(unique(as.factor(data_outcome_proxies$target_snp.outcome))), " of ", length(unique(as.factor(proxies$target_snp.outcome))), " missing-SNP(s) extracted"))

  # select proxy-SNP(s) with highest R2 ====
  data_outcome_proxies <- data_outcome_proxies %>%
    group_by(target_snp.outcome) %>%
    filter(R == max(R)) %>%
    slice(1) %>%
    select(-R)

  ## Bind rows of data_outcome with data_outcome_proxies
  data_outcome <- bind_rows(data_outcome, data_outcome_proxies)

  return(data_outcome)
}
