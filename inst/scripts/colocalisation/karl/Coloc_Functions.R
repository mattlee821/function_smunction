### Load in the libraries we need 
library(dplyr)
library(tidyverse)
library(LDlinkR)
library(ggrepel)
library(ggsci)
library(coloc)
library(locuscomparer)
library(TwoSampleMR)
library(phenoscanner)
library(ieugwasr)


#source("~/Documents/Misc_Functions/MR_Functions.R")


####################### Contributing Functions

# thresh = 0.8
find_proxies <- function(snp,LD, thresh = 0.8){  
  list_proxies <- LD %>% select(.,all_of(snp)) %>% 
    filter(. > thresh) %>% filter(row.names(.) %in% outcome_dat$SNP) %>%
    arrange(desc(.)) %>% dplyr::slice(1) %>%
    mutate(SNP = colnames(.), Proxy = row.names(.)) %>%
    select(-1)
  return(list_proxies)
}

get_proxy_out_dat_LD <- function(exp_dat, outcome_dat,LD, thresh = 0.8){
  ## Find out which SNPs are in the lung dat and which need proxies 
  exp_dat_no_Proxy <- outcome_dat[outcome_dat$SNP %in% exp_dat$SNP,]
  
  ## Lets look for proxies for those SNPs that we dont have in our data
  SNPs_to_be_proxied <- exp_dat[!exp_dat$SNP %in% exp_dat_no_Proxy$SNP,]
  
  proxy_snps <- data.frame()
  for (snp in SNPs_to_be_proxied$SNP) {
    proxy_snps <- rbind(find_proxies(snp = snp,LD, thresh = thresh),proxy_snps)
  }
  
  
  filtered_proxies <- data.frame()
  for (snp in unique(proxy_snps$Proxy)) {
    chosen_dup <- exp_dat[exp_dat$SNP %in% proxy_snps$SNP[proxy_snps$Proxy == snp],] %>% 
      dplyr::slice(which.min(pval.exposure)) %>% select(SNP)
    filtered_proxies <- rbind(proxy_snps[proxy_snps$SNP == chosen_dup$SNP,],filtered_proxies) 
  }
  
  snp_for_proxy <- outcome_dat[outcome_dat$SNP %in% filtered_proxies$Proxy,]
  
  colnames(filtered_proxies) <- c("Proxy", "SNP")
  snp_for_proxy <- merge(snp_for_proxy, filtered_proxies)
  
  formated_out_dat_proxy <- format_proxies(snp_for_proxy)
  
  ## combine with proxy SNPs
  out_dat_All <- rbind(test_formating, exp_dat_no_Proxy)
  
  
  Harm_dat <- harmonise_data(exp_dat, out_dat_All)
  
  return(Harm_dat)
}

  
### This function returns the LD matrix
get_LD_mat_ldlink <- function(exp_dat, token ="4134a0a348dc" ){
  ### we want to get the LD matrix for the SNPs in the region
  LD_Full <- LDlinkR::LDmatrix(exp_dat[,"SNP"], pop = "CEU", r2d = "r2", token = "4134a0a348dc", file = FALSE)
  
  ## format the LD matrix
  LD <- LD_Full[,-1]
  LD <- LD[,!colnames(LD) %in% names(which(colSums(is.na(as.matrix(LD))) > dim(LD)-1)) ]
  LD <- LD[complete.cases(LD),]
  rownames(LD) <- colnames(LD)
  return(list(LD_Anal = LD, LD_Plot = LD_Full))
}


get_coloc_all <- function(Harm_dat, LD, N_exp, N_out, prop_cases){
  ### check that the N and prop and eaf are present
  
  ## Make sure there are no duplicates
  Harm_dat <- Harm_dat[!duplicated(Harm_dat$SNP),]
  
  Harm_dat <- Harm_dat %>% filter(SNP %in% colnames(LD))
  
  ## do some formatting
  LD <- data.frame(columnNameILike = row.names(LD), LD)
  LD <- LD[,-1]
  ## keep only the SNPs in the LD matrix 
  LD <- LD %>% dplyr::select(Harm_dat$SNP) %>% filter(row.names(.) %in% colnames(.))
  LD <- LD[order(row.names(LD)),] %>% select(row.names(.))
  ## format the input dataframes
  ### Exposure
  D1=list(beta=Harm_dat$beta.exposure,
          varbeta=Harm_dat$se.exposure^2,
          snp=Harm_dat$SNP, 
          MAF=Harm_dat$eaf.outcome, 
          LD = as.matrix(LD),
          position = Harm_dat$pos.exposure,
          N=N_exp, 
          sdY=1,
          type="quant")
  
  ### Outcome
  D2=list(beta=Harm_dat$beta.outcome,
          varbeta=Harm_dat$se.outcome^2,
          snp=Harm_dat$SNP, 
          MAF=Harm_dat$eaf.outcome, 
          LD = as.matrix(LD),
          position = Harm_dat$pos.exposure,
          s=prop_cases, # This is the case proportion. Make sure to specify if using a case-control GWAS (get this information from the paper
          N=N_out, #number of participants
          type="cc")
  
  ### Run the conditional iterative method
  colo_res_Cond_Iter <- coloc::coloc.signals(
    D1,
    D2,
    MAF = NULL,
    LD = NULL,
    method = c("cond"),
    mode = c("iterative"),
    p1 = 1e-03,
    p2 = 1e-04,
    p12 = 1e-05,
    maxhits = 3,
    r2thr = 0.01,
    pthr = 1e-06
  )
  colo_res_Cond_Iter$summary$method <- "Cond_Iter"
  
  ### Run the single method
  colo_res_simple <- coloc::coloc.signals(
    D1,
    D2,
    MAF = NULL,
    LD = NULL,
    method = c("single"),
    p1 = 1e-03,
    p2 = 1e-04,
    p12 = 1e-05,
    maxhits = 3,
    r2thr = 0.01,
    pthr = 1e-06
  )
  colo_res_simple$summary$method <- "Single"
  
  ## rbind the two methods
  comb_coloc <- rbind(colo_res_simple$summary,colo_res_Cond_Iter$summary )
  
  ## re-order columns 
  comb_coloc <- comb_coloc[,c("nsnps", "hit1", "hit2",
                              "PP.H0.abf","PP.H1.abf",
                              "PP.H2.abf","PP.H3.abf",
                              "PP.H4.abf", "method")]
  #######################
  #######################
  ## run susie
  #######################
  ## Do the finemapping
  # S1=tryCatch(expr = 
  #   {temp <- coloc::runsusie(D1,repeat_until_convergence = F, 
  #                    maxit = 1000)}, 
  #       error = function(e) { 
  #       message("Susie Didn't Converge for Exposure")
  #       NULL
  #   })
  # 
  # S2=tryCatch(expr = 
  #   {temp <- coloc::runsusie(D2,repeat_until_convergence = F, 
  #                    maxit = 1000)}, 
  #       error = function(e) { 
  #       message("Susie Didn't Converge for Outcome")
  #       NULL
  #   })
  # if(!is.null(S1) & !is.null(S2)){
  #   ## Run Susie Coloc
  #   susie.res=coloc::coloc.susie(S1,S2,
  #                              p1 = 1e-03,
  #                              p2 = 1e-04,
  #                              p12 = 1e-05)
  # } else {
  #     susie.res <- NULL
  #   }
  # if(!is.null(susie.res$summary)){
  #   ### Some Results formatting
  #   ## select the columns we want
  #   susie.res$summary <- susie.res$summary[, 1:8]
  #   ## create a methods variabe
  #   susie.res$summary$method <- "Susie"
  #   ########################
  #   ########################
  #   ## Combine with existing results
  #   comb_coloc <- rbind(comb_coloc,susie.res$summary )
  # } else {
  #   print("Susie found nothing")
  # }
  
  ## Give the columns meaningful names for the exposure and the outcome
  comb_coloc$Exposure <- unique(Harm_dat$exposure)
  comb_coloc$Outcome <- unique(Harm_dat$outcome)
  
  ## create a nice summary dataframe 
  summary_return <- comb_coloc %>% group_by(Outcome, Exposure, method) %>% filter(PP.H4.abf == max(PP.H4.abf)) %>% 
    select(PP.H4.abf, method, Exposure, Outcome) %>%
    pivot_wider(names_from = method, values_from = PP.H4.abf)
  
  return(list(Summary = summary_return, Results = comb_coloc))
  
}


Make_ZZ_Plot <- function(LD_Mat, lead_SNP, Harm_dat, 
                         exposure_name =  base::unique(Harm_dat[,"exposure"]),
                         outcome_name = base::unique(Harm_dat[,"outcome"])) {

  # if(!lead_SNP %in% colnames(LD_Mat)){
  #   print("Lead SNP not in LD Matrix")
  #   lead_SNP <- Harm_dat %>% filter(SNP %in% colnames(LD_Mat)) %>%
  #     arrange(desc(.)) %>% dplyr::slice(1) %>% select(SNP)
  #   }
  
  ## We want to extract the column that is the lead SNP for plotting
  LD_TEMP <- LD_Mat[,c("RS_number", lead_SNP)]
  
  ## We want to merge through the LD vector for plottign
  temp_dat_Format <- merge(Harm_dat, LD_TEMP, by.x ="SNP" , by.y = "RS_number", all.x = T)
  
  ## Annoyingly long ifelse to create color factor variable
  temp_dat_Format[,lead_SNP] <- ifelse(is.na(temp_dat_Format[,lead_SNP]), 0, temp_dat_Format[,lead_SNP])
  temp_dat_Format$LD <- ifelse(temp_dat_Format[,lead_SNP]  > 0 &temp_dat_Format[,lead_SNP]  <=.2, "LD < 0.2",
                               ifelse(temp_dat_Format[,lead_SNP]  > 0.2 & temp_dat_Format[,lead_SNP]  <=.4, "0.2 > LD < 0.4",
                                      ifelse(temp_dat_Format[,lead_SNP]  > 0.4 & temp_dat_Format[,lead_SNP] <=.6, "0.4 > LD < 0.6",
                                             ifelse(temp_dat_Format[,lead_SNP]  > 0.6 & temp_dat_Format[,lead_SNP]  <=.8, "0.6 > LD < 0.8",
                                                    ifelse(temp_dat_Format[,lead_SNP]  > 0.8, "LD > 0.8", "No LD")))))
  temp_dat_Format$LD <- ifelse(is.na(temp_dat_Format$LD), "No LD", temp_dat_Format$LD)
  temp_dat_Format$LD <- ifelse(temp_dat_Format$SNP == lead_SNP, "Lead SNP", temp_dat_Format$LD)
  
  ## Create Z scores for plotting
  temp_dat_Format$Z_exp <- temp_dat_Format$beta.exposure/temp_dat_Format$se.exposure
  temp_dat_Format$Z_out <- temp_dat_Format$beta.outcome/temp_dat_Format$se.outcome
  
  
  pos <- position_jitter(width = 0.5, seed = 1)
  # Reorder following the value of another column:
  Plot <- temp_dat_Format %>% mutate(LD = fct_reorder(LD, get(lead_SNP))) %>%
    ggplot(aes(Z_exp, Z_out, color = LD)) + geom_point(size = 2) + 
    theme_bw() + xlab(paste(exposure_name," Z-score", sep = "")) +
    ylab(paste(outcome_name," Z-score", sep = "")) + 
    ggtitle(paste("Z-Z Locus Plot for: ", exposure_name, " and ", outcome_name, sep = "")) +
    theme(axis.text = element_text(hjust = 1,  size =20),
          plot.title = element_text(hjust = 0.5,size=22,face="bold"),
          axis.title=element_text(size=25,face="bold"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20, face = "bold")) + 
    geom_label_repel(size = 6,
                     data=temp_dat_Format %>% filter(SNP==lead_SNP), # Filter data first
                     aes(label=SNP),show.legend = FALSE
    )+ 
    labs(color= 'LD with Lead SNP') +  
    scale_color_manual(values=c("No LD" = "#D3D3D3",
                                "LD < 0.2"="#00468BB2",
                                "0.2 > LD < 0.4"="#0099B4B2",
                                "0.4 > LD < 0.6"="#20854EB2",
                                "0.6 > LD < 0.8"= "#DF8F44B2",
                                "LD > 0.8"="#ED0000B2",
                                "Lead SNP" = "#cb99fe"
    ))
  return(Plot)
}


Make_Locus_Compare <- function(Harm_dat,lead_SNP, 
                               exposure_name =  unique(Harm_dat[,"exposure"]),
                               outcome_name = unique(Harm_dat[,"outcome"])){
  ## Make  LOCUS ZOOM plot 
  ## we want to recalculate the p values
  Harm_dat$pval.exposure <- pnorm(abs(Harm_dat$beta.exposure)/Harm_dat$se.exposure, lower.tail = FALSE) * 2
  Harm_dat$pval.outcome <- pnorm(abs(Harm_dat$beta.outcome)/Harm_dat$se.outcome, lower.tail = FALSE) * 2
  exp_gwas_lc <- Harm_dat[,c("SNP", "pval.exposure")]
  colnames(exp_gwas_lc) <- c("rsid", "pval")
  
  out_gwas_lc <- Harm_dat[,c("SNP", "pval.outcome")]
  colnames(out_gwas_lc) <- c("rsid", "pval")
  
  locus_plot <- tryCatch(expr =  
             {locuscompare(in_fn1 = exp_gwas_lc, in_fn2 = out_gwas_lc, title = exposure_name, 
                           title2 = outcome_name)},
           error = function(e) {
             message("Locus Plot Didnt Work...")
             ggplot(Harm_dat, aes(beta.exposure, beta.outcome)) + geom_point()+
               ggtitle("Replacement for Locus Plot where SNPs 'not on same chromosome'")
           })

  return(locus_plot)
  
}



Do_Coloc_Plots <- function(Harm_dat, LD, N_exp, 
                           N_out, prop_cases, 
                           lead_SNP){
  ## Do the coloc
  coloc_res <- get_coloc_all(Harm_dat ,LD[[1]], N_exp, N_out, prop_cases)
  ## get the zz plot
  ZZ_Plot <- Make_ZZ_Plot(LD[[2]], lead_SNP, Harm_dat = Harm_dat)
  ## get the locus plot
  locus_Plot <- Make_Locus_Compare(Harm_dat, lead_SNP)
  
  return(list(Harmonised_Dataframe = Harm_dat, Coloc_Results_Summary = coloc_res[[1]],
              Coloc_Results =  coloc_res[[2]], Z_Z_Plot=ZZ_Plot, LocusComparePlot = locus_Plot))
  #return(list(Harmonised_Dataframe = Harm_dat, Coloc_Results_Summary = coloc_res[[1]],
  #       Coloc_Results =  coloc_res[[2]]))
}

Coloc_LD_Exp_Format <- function(exp_dat){
  #if(dim(exp_dat)[1] > 1000) {
  #  print(paste0("Chr is:",unique(exp_dat$chr.exposure)," there are ", dim(exp_dat)[1], " SNPs. We will clump at 0.999 to get it down and exclude exact signal duplications" ))
  #  exp_dat <- clump_data(exp_dat, clump_r2 = 0.9999)
  #  print(paste0("After clumping there are ",dim(exp_dat)[1], " SNPs" ))
  #}
  LD_Return <- get_ld_matrix(exp_dat$SNP,plink_loc = plink_loc, bfile_loc = bfile_loc)
  ## restrict only to SNPs with LD values
  exp_dat <- exp_dat[exp_dat$SNP %in% colnames(LD_Return[[1]]), ]
  exp_dat <- exp_dat <- exp_dat[order(colnames(LD_Return[[1]])),]
  return(list(LD=LD_Return, exp_dat = exp_dat))
}

Coloc_LD_Exp_Format_LD <- function(exp_dat){
  #if(dim(exp_dat)[1] > 1000) {
  #  print(paste0("Chr is:",unique(exp_dat$chr.exposure)," there are ", dim(exp_dat)[1], " SNPs. We will clump at 0.999 to get it down and exclude exact signal duplications" ))
  #  exp_dat <- clump_data(exp_dat, clump_r2 = 0.9999)
  #  print(paste0("After clumping there are ",dim(exp_dat)[1], " SNPs" ))
  #}
  LD_Return <- get_LD_mat_ldlink(exp_dat)
  ## restrict only to SNPs with LD values
  exp_dat <- exp_dat[exp_dat$SNP %in% colnames(LD_Return[[1]]), ]
  exp_dat <- exp_dat <- exp_dat[order(colnames(LD_Return[[1]])),]
  return(list(LD=LD_Return, exp_dat = exp_dat))
}

bfile_loc <- "QC_1000G_P3"
plink_loc <- "./plink"
get_ld_matrix <- function(rsid_list, plink_loc, bfile_loc, with_alleles = F){
  setwd("/Users/karlsm/Documents/Genetics/1000G_EUR/")
  LD_Full <- ieugwasr::ld_matrix_local(rsid_list, bfile = bfile_loc, 
                                       plink_bin= plink_loc , with_alleles = F)
  
  ## format the LD matrix
  LD <- LD_Full[,!colnames(LD_Full) %in% names(which(colSums(is.na(as.matrix(LD_Full))) > dim(LD_Full)-1)) ]
  LD <- LD[complete.cases(LD),]
  rownames(LD) <- colnames(LD)
  LD_Full <- as.data.frame(LD_Full)
  LD_Full$RS_number <- rownames(LD)
  return(list(LD_Anal = LD, LD_Plot = LD_Full))
  
}


Master_Coloc <- function(exp_dat= NA, exposure_id = NA, 
                         outcome_dat = NA, outcome_id = NA, 
                         lead_SNP, window = 75000, token ="4134a0a348dc", 
                         N_exp, N_out, prop_cases){
  
  if(is.na(exposure_id)){
    print("Internal exposure dataset")
  } else if(!is.na(exposure_id)){ 
    print("Fetchting external exposure dataset")
    exp_dat <- get_openGWAS_protein(lead_SNP, exposure_id, window = window)
  } else {
    stop("Neither External Exposure Data nor OpenGWAS ID provided")
  }
  ## keep only variants we can query in LDLinkR
  exp_dat <- exp_dat %>% filter(str_detect(SNP,"^rs"))
  
  if(dim(exp_dat[is.na(exp_dat$eaf.exposure),])[1] > 0) {
    print("Fetching EAF...")
    exp_dat <- EAF_SNP(exp_dat)
  }
  dim_exp_pre <- dim(exp_dat)[1]
  exp_dat <- exp_dat %>% filter(!is.na(eaf.exposure))
  dim_exp_post <- dim(exp_dat)[1]
  print(paste0(abs(dim_exp_post-dim_exp_pre), " SNPs Excluded because no EAF"))
  
  LD_and_Exp_Dat <- Coloc_LD_Exp_Format_LD(exp_dat)
  
  if(is.na(outcome_id)){
    print("Internal Outcome dataset")
    Harm_dat <- harmonise_data(LD_and_Exp_Dat$exp_dat, outcome_dat)
  } else if(!is.na(outcome_id)){ 
    print("Fetchting external outcome dataset")
    outcome_Open_GWAS <- extract_outcome_data(LD_and_Exp_Dat$exp_dat$SNP,outcome_id )
    Harm_dat <- harmonise_data(exp_dat, outcome_Open_GWAS)
  } else {
    stop("Neither External Outcome Data nor OpenGWAS ID provided")
  }
  
  if(sum(is.na(Harm_dat$other_allele.exposure))) {

    Harm_dat$other_allele.exposure <- Harm_dat$other_allele.outcome
  }
  if(sum(is.na(Harm_dat$other_allele.outcome))) {
    Harm_dat$other_allele.outcome <- Harm_dat$other_allele.exposure
  }

  
  if(dim(Harm_dat[is.na(Harm_dat$eaf.outcome),])[1] < 
     dim(Harm_dat[is.na(Harm_dat$eaf.exposure),])[1]) {
    Harm_dat$eaf.exposure <- Harm_dat$eaf.outcome
  } else {
    Harm_dat$eaf.outcome<- Harm_dat$eaf.exposure
  }
  
  Harm_dat <- Harm_dat[!is.na(Harm_dat$eaf.outcome),]
  
  if(!lead_SNP %in% Harm_dat$SNP){
    print("Lead SNP not in Outcome Data. Finding Proxy for Plots")
    
    out_dat_Proxy <- tryCatch(expr =  
                             {find_proxies(snp = lead_SNP,LD = LD_and_Exp_Dat$LD[[2]])},
                             error = function(e) {
                               message("we couldnt find the lead SNP in the LD Mat")
                               data.frame()
                             })
    if(dim(out_dat_Proxy)[1] < 1){
      print("Lead SNP doesnt have good proxy for plot. Taking label for min(p exposure)")
      out_snp <- Harm_dat %>% filter(SNP %in% colnames(LD_and_Exp_Dat$LD[[2]])) %>%
        mutate(temp_z = beta.exposure/se.exposure) %>%
        arrange(desc(abs(temp_z))) %>% dplyr::slice(1) %>% select(SNP)
    } else {
      out_snp <- out_dat_Proxy$Proxy
    }
    
    print(paste("Lead SNP replacement is: ",  out_snp, sep = ""))
    lead_SNP <- out_snp
  }
  
  lead_SNP <- as.character(lead_SNP)
  results_to_return <- Do_Coloc_Plots(Harm_dat, 
                                      LD_and_Exp_Dat$LD, N_exp, 
                                      N_out, prop_cases, lead_SNP)
  
  names(results_to_return) <-paste(names(results_to_return),"_", lead_SNP, sep= "")
  return(results_to_return)
}







