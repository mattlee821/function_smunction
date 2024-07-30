## Functions 
EAF_SNP <- function(in_dat){
  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  library(rsnps)

  if (dim(in_dat)[1] > 100) {
    print("Chunking data...")
    temp <- data.frame()
    i=2
    split_dat <-chunk(in_dat, i)
    dims_dat <- data.frame()
    for(i in 1:length(split_dat)){
      temp_dim <- dim(split_dat[[i]])[1]
      dims_dat <- rbind(dims_dat, temp_dim)
    }
    max(dims_dat)
    while(max(dims_dat) >100){
    i = i + 1
      split_dat <-chunk(in_dat, i)
      dims_dat <- data.frame()
      for(i in 1:length(split_dat)){
        temp_dim <- dim(split_dat[[i]])[1]
        dims_dat <- rbind(dims_dat, temp_dim)
      }
    }
    
    
    

    for (d in split_dat) {
      
    temp_interim <- phenoscanner::phenoscanner(d$SNP)$snp
    temp <- rbind(temp, temp_interim)
    
    }
    if(dim(temp)[1]>0) {
      temp$a1 <- tolower(temp$a1)
      temp$a2<- tolower(temp$a2)
      temp$eur <- as.numeric(as.character(temp$eur))
      
      ret_dat <- merge(in_dat, temp[, c("snp", "a1", "a2", "eur")], by.x = "SNP", by.y = "snp", all.x = T)
      
      ret_dat$eaf.exposure <- ifelse(tolower(ret_dat$effect_allele.exposure) == ret_dat$a1,ret_dat$eur,
                                     ifelse(tolower(ret_dat$effect_allele.exposure) == ret_dat$a2,abs(ret_dat$eur-1),NA))
      ret_dat <- ret_dat[, 1:13]
    }
    else {
      temp_snp_2 <- ncbi_snp_query(in_dat$SNP)
      
      ret_dat <- merge(in_dat, temp_snp_2[, c("query", "minor", "ref_seq", "maf")], by.x = "SNP", by.y = "query", all.x = T)
      
      ret_dat$eaf.exposure <- ifelse(tolower(ret_dat$effect_allele.exposure) == ret_dat$minor,ret_dat$maf,
                                     ifelse(tolower(ret_dat$effect_allele.exposure) == ret_dat$ref_seq,abs(ret_dat$maf-1),NA))
      
      ret_dat <- ret_dat[, 1:13]
    }
  } else {
  temp <- phenoscanner::phenoscanner(d$SNP)$snp
  
  if(dim(temp)[1]>0) {
    temp$a1 <- tolower(temp$a1)
    temp$a2<- tolower(temp$a2)
    temp$eur <- as.numeric(as.character(temp$eur))
    
    ret_dat <- merge(in_dat, temp[, c("snp", "a1", "a2", "eur")], by.x = "SNP", by.y = "snp", all.x = T)
    
    ret_dat$eaf.exposure <- ifelse(tolower(ret_dat$effect_allele.exposure) == ret_dat$a1,ret_dat$eur,
                                   ifelse(tolower(ret_dat$effect_allele.exposure) == ret_dat$a2,abs(ret_dat$eur-1),NA))
    ret_dat <- ret_dat[, 1:13]
  }
  else {
    temp_snp_2 <- ncbi_snp_query(in_dat$SNP)
    
    ret_dat <- merge(in_dat, temp_snp_2[, c("query", "minor", "ref_seq", "maf")], by.x = "SNP", by.y = "query", all.x = T)
    
    ret_dat$eaf.exposure <- ifelse(tolower(ret_dat$effect_allele.exposure) == ret_dat$minor,ret_dat$maf,
                                   ifelse(tolower(ret_dat$effect_allele.exposure) == ret_dat$ref_seq,abs(ret_dat$maf-1),NA))
    
    ret_dat <- ret_dat[, 1:13]
  }
  }
  
  return(ret_dat)
}





## get proxies
get_proxy_snps <- function(snp, r2_thresh = 0.8, GWAS_dat) {
  # Check if package is installed, and if not then grab it.
  if("LDlinkR" %in% rownames(installed.packages()) == FALSE) {install.packages("LDlinkR")}
  # Load package into environment
  library("LDlinkR")
  library(dplyr)
  library(data.table)
  
  # Grab r2 dataframe
  r2_dat <- tryCatch(expr = {LDproxy(snp,pop = "CEU", r2d = "r2", 
                                     token = "4134a0a348dc")},
                     error = function(e) { 
                       message("we had an error")
                       replace_data <- LDproxy("rs10993994",pop = "CEU", r2d = "r2", token = "4134a0a348dc")
                       replace_data[-c(1:dim(replace_data)[1]),]
                     }
  )
  if(grepl("error",r2_dat[1,1],fixed=T)){
    print("SNP is not a biallelic variant")
  }else {
  # Remove the target SNP
  r2_dat <- r2_dat[r2_dat$RS_Number != snp,]
  # remove snps that are not above the ld threshold.
  r2_dat <- r2_dat[r2_dat$R2> r2_thresh,]
  
  # remove non-rsids
  r2_dat <- r2_dat[r2_dat$RS_Number %like% "rs",]
  # Check if we indeed have any proxies
  if(dim(r2_dat)[1] ==0) {
    print("No Proxies above specified threshold")
  }
  
  # We want to rank the dataframe by the strength of
  # the r2 that a given proxy has with our target SNP
  else {
    print(paste("We have", dim(r2_dat)[1], "potential proxy SNPs", sep = " "))
    # Grab order
    order.r2<- order(-r2_dat$R2,r2_dat$RS_Number)
    # Create empty vector 
    r2_dat$rank <- NA
    # Insert ranking variable
    r2_dat$rank[order.r2] <- 1:nrow(r2_dat)
    # Check which, if any of our proxies are in the 
    temp_dat <- GWAS_dat[GWAS_dat$SNP %in% r2_dat$RS_Number,]
    # Check if our proxies have any betas in the target GWAS
    if(dim(temp_dat)[1] ==0){
      print(paste("There are no available betas for any proxy SNPs of", snp, sep = " "))
    }
    else{
      temp_dat <- merge(temp_dat, r2_dat[,c("RS_Number", "rank")], by.x = "SNP", by.y = "RS_Number", all.x == T)
      temp_dat <- temp_dat[which.min(temp_dat$rank),]
      temp_dat$Proxy = snp
      print(paste(temp_dat$SNP,"is your proxy for", snp,", have fun", sep = " "))
      return(temp_dat)
    }
  }
  }
}



format_proxies <- function(input_dat) {
  library(LDlinkR)
  return_dat <- data.frame()
  for (i in seq(1, length(input_dat$SNP), 1)) {
    print(input_dat$SNP[i])
    dat <- input_dat[i,]
    ## Grab out the ldmatrix
    ld_mat <- tryCatch({LDpair(dat$SNP, dat$Proxy,  token = "4134a0a348dc")}, error=function(e){})
  
    ## grab out the allele frequency of the effect allele for the proxy variant
    if(dat$effect_allele.outcome ==  ld_mat$var1_a1){
      effect.allele.maf <-  ld_mat$var1_a1_freq
    } else {
      effect.allele.maf <-  ld_mat$var1_a2_freq
    }
    ## assign the correct allele for the proxy variant based on A1 MAF and the SNP that we have in the exposure data
    if(abs(effect.allele.maf-ld_mat$var2_a1_freq) < abs(effect.allele.maf-ld_mat$var2_a2_freq)) {
      dat$effect_allele.outcome <- ld_mat$var2_a1
      dat$other_allele.outcome <- ld_mat$var2_a2
      
    } else {
      dat$effect_allele.outcome <- ld_mat$var2_a2
      dat$other_allele.outcome <- ld_mat$var2_a1
    }
    ## grab the eaf 
    dat$eaf.outcome <- effect.allele.maf
    ## remove the extra columns
    #dat <- dat[,-c(13:14)]
    ## rbind it with the return_dat
    return_dat <- rbind(return_dat, dat)
    
  }
  return(return_dat)
  
}


get_rsid_pos_simp <- function(SNP, window_size = 75000) {
  
  ## Get the chr and position
  temp <- tryCatch(expr = {ieugwasr::variants_rsid(SNP)[c("chr", "pos")]},
                       error = function(e) { 
                         message("we had an error")
                         replace_data <- phenoscanner::phenoscanner(SNP)$snps[c("chr", "pos_hg19")]
                       }
                       
                       
                       
  )
  
  list_ret <- list(Chr = paste("chr",temp[,1], sep = ""), Position = as.numeric(temp[,2]))
  return(list_ret)
}







get_rsid_pos <- function(SNP, window_size = 75000) {
  
  ## Get the chr and position
  temp <- tryCatch(expr = {ieugwasr::variants_rsid(SNP)[c("chr", "pos")]},
                   error = function(e) { 
                     message("we had an error")
                     replace_data <- phenoscanner::phenoscanner(SNP)$snps[c("chr", "pos_hg19")]
                   }
                   
                   
                   
  )
  
  #BiocManager::install("rtracklayer")
  library(rtracklayer)
  session <- browserSession("UCSC")
  genome(session) <- "hg19"
  #trackNames(session) ## list the track names
  temp_chr <- paste("chr",temp[,1], sep = "")
  temp_pos <- temp[,2]
  ## choose the Conservation track for a portion of mm9 chr1
  query <- ucscTableQuery(session, "dbSnp153Composite",
                          GRangesForUCSCGenome("hg19", temp_chr,
                                               IRanges(as.numeric(temp_pos)-window_size, 
                                                       as.numeric(temp_pos)+window_size)))
  ## get the phastCons30way track
  tableName(query) <-  "dbSnp153Common"
  SNPs <- getTable(query)
  list_ret <- list(Chr = paste("chr",temp[,1], sep = ""), Position = as.numeric(temp[,2]))
  list_ret <- list(SNP_Coord = SNPs,Chr = paste("chr",temp[,1], sep = ""), Position = as.numeric(temp[,2]))
  return(list_ret)
}


format_chr_pos <- function(temp_dat, var_name = "MarkerName" , delim = ":", window = window, temp_coord_dat = temp_coord) {
  ## check if the var_name is in the data frame 
  if(!var_name %in% colnames(temp_dat)) {
    print("Chr:Pos variable not in dataframe")
  } else {
    
    ## generate the chromosome
    temp_dat$Chr <- data.frame(str_split_fixed(temp_dat[, get(var_name)], delim, 3))[,1]
    ## format chr variable
    temp_dat <- temp_dat[temp_dat$Chr ==gsub("chr", "", temp_coord_dat["Chr"]),]
    ## generate the position variable and make it numeric
    temp_dat$Position <- data.frame(str_split_fixed(temp_dat[, get(var_name)], delim, 3))[,2]
    temp_dat$Position  <- as.numeric(as.character(temp_dat$Position))
    
    ## Extract window
    temp_dat <- temp_dat[ between(Position, abs(as.numeric(temp_coord_dat["Position"])-window), 
                                  abs(as.numeric(temp_coord_dat["Position"])+window)),]
    
    temp_dat <- merge(temp_dat, temp_coord_dat$SNP_Coord, by.x = "Position", by.y = "chromEnd")
    
    
    temp_dat_af <- ieugwasr::afl2_rsid(temp_dat$name)
    names(temp_dat_af)[names(temp_dat_af) == 'ref'] <- 'REF'
    temp_dat <- merge(temp_dat, temp_dat_af, by.x = "name", by.y = "id")
    
    #temp_dat$eaf <- ifelse(temp_dat$Allele1 == tolower(temp_dat$REF), as.numeric(temp_dat$AF.EUR), 
    #                      1-as.numeric(temp_dat$AF.EUR))
    ## return the dataframe
    return(temp_dat)
  }
}


Get_Coloc <- function(file_addr, window, SNPs, exposure_name,outcome_name, outcome_dat_id) {
  
  library(stringr)
  library(TwoSampleMR)
  library(LDlinkR)
  library(coloc)
  library(ggplot2)
  library(ggsci)
  library(ggrepel)
  library(dplyr)
  library(forcats)
  library(locuscomparer)
  
  ## read in the files
  temp_dat_full <- fread(file_addr)
  
  return_list <- list()
  for(lead_SNP in SNPs ) {
    print(c(lead_SNP,exposure_name))
    ## Load in the coord
    temp_coord <- get_rsid_pos(lead_SNP, window)
    
    
    
    
    ## format the data and prune to the window that we want
    ## temp_dat <- format_chr_pos(temp_dat_full, "MarkerName", ":", window,
    ## temp_coord)
    
    temp_dat <- temp_dat_full[temp_dat_full$chr == gsub("chr", replacement = "", temp_coord$Chr) &
                                between(pos, abs(as.numeric(temp_coord["Position"])-window), 
                                        abs(as.numeric(temp_coord["Position"])+window)),]
    
    ## convert to the twosample format
    temp_dat_Format <- format_data(temp_dat, 
                                   snp_col = "SNP" ,
                                   beta_col = "Effect",
                                   se_col = "StdErr",
                                   effect_allele_col = "Allele1" ,
                                   other_allele_col =  "Allele2",
                                   pval_col = "P-value"
    )
    ## give things a name
    temp_dat_Format$exposure <- exposure_name
    
    
    
    temp_dat_Format <- temp_dat_Format[temp_dat_Format$effect_allele.exposure %in% c("A", "C", "G", "T") &
                                         temp_dat_Format$other_allele.exposure %in% c("A", "C", "G", "T"), ]
    ## Grab out the outcome data
    outcome_dat <- TwoSampleMR::extract_outcome_data(temp_dat_Format$SNP, outcome_dat_id)
    
    ## Harmonise the data
    Harm_dat <- harmonise_data(temp_dat_Format,outcome_dat)
    
    ## Grab out the LD
    LD <- LDmatrix(Harm_dat$SNP, pop = "CEU", r2d = "r2", token = "4134a0a348dc", file = FALSE)
    
    ## We want to extract the column that is the lead SNP for plotting
    LD_TEMP <- LD[,c("RS_number", lead_SNP)]
    
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
      theme_bw() + xlab(paste(exposure_name," Z score", sep = "")) +
      ylab(paste(outcome_name," Z score", sep = "")) + 
      #ggtitle("Cis eQTL in Normal Lung Tissue") +
      theme(axis.text = element_text(hjust = 1,  size =20),
            plot.title = element_text(hjust = 0.5,size=22,face="bold"),
            axis.title=element_text(size=25,face="bold"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20, face = "bold")) + 
      geom_label_repel(size = 6,
                       data=temp_dat_Format %>% filter(SNP==lead_SNP), # Filter data first
                       aes(label=SNP),show.legend = FALSE
      )+ 
      labs(color= 'LD with Lead Cis SNP') +  
      scale_color_manual(values=c("No LD" = "#D3D3D3",
                                  "LD < 0.2"="#00468BB2",
                                  "0.2 > LD < 0.4"="#0099B4B2",
                                  "0.4 > LD < 0.6"="#20854EB2",
                                  "0.6 > LD < 0.8"= "#DF8F44B2",
                                  "LD > 0.8"="#ED0000B2",
                                  "Lead SNP" = "#cb99fe"
      ))
    
    
    ## format the LD matrix
    LD <- LD[,-1]
    LD <- LD[,!colnames(LD) %in% names(which(colSums(is.na(as.matrix(LD))) > dim(LD)-1)) ]
    LD <- LD[complete.cases(LD),]
    rownames(LD) <- colnames(LD)
    ## restrict only to SNPs with LD values
    Harm_dat <- Harm_dat[Harm_dat$SNP %in% colnames(LD), ]
    Harm_dat <- Harm_dat <- Harm_dat[order(colnames(LD)),]
    
    ## Make sure that we have eaf for all SNPs
    Harm_dat <- Harm_dat[!is.na(Harm_dat$eaf.outcome),]
    
    
    ## Grab the position from the temp_dat_Format above
    
    #Harm_dat <- merge(Harm_dat, temp_dat[,c("name", "pos")], by.x = "SNP", by.y = "name", all.x = T)
    Harm_dat <- Harm_dat[!duplicated(Harm_dat$SNP),]
    D1=list(beta=Harm_dat$beta.exposure,
            varbeta=Harm_dat$se.exposure^2,
            snp=Harm_dat$SNP, 
            MAF=Harm_dat$eaf.outcome, 
            LD = as.matrix(LD),
            position = Harm_dat$pos.exposure,
            N=40000, 
            sdY=1,
            type="quant")
    
    
    D2=list(beta=Harm_dat$beta.outcome,
            varbeta=Harm_dat$se.outcome^2,
            snp=Harm_dat$SNP, 
            MAF=Harm_dat$eaf.outcome, 
            LD = as.matrix(LD),
            position = Harm_dat$pos.exposure,
            s=0.20, # This is the case proportion. Make sure to specify if using a case-control GWAS (get this information from the paper
            N=244230, #number of participants
            type="cc")
    
    
    colo_res <- coloc.abf(dataset1 = D1, dataset2 = D2)
    
    
    ## Make  LOCUS ZOOM plot 
    
    ## we want to recalculate the p values
    Harm_dat$pval.exposure <- pnorm(abs(Harm_dat$beta.exposure)/Harm_dat$se.exposure, lower.tail = FALSE) * 2
    Harm_dat$pval.outcome <- pnorm(abs(Harm_dat$beta.outcome)/Harm_dat$se.outcome, lower.tail = FALSE) * 2
    exp_gwas_lc <- Harm_dat[,c("SNP", "pval.exposure")]
    colnames(exp_gwas_lc) <- c("rsid", "pval")
    
    out_gwas_lc <- Harm_dat[,c("SNP", "pval.outcome")]
    colnames(out_gwas_lc) <- c("rsid", "pval")
    
    locus_plot <- locuscompare(in_fn1 = exp_gwas_lc, in_fn2 = out_gwas_lc, title = exposure_name, 
                               title2 = outcome_name, snp = lead_SNP)
    
    ## save list to be returned
    return_dat <- list(Plot, colo_res, locus_plot)
    names(return_dat) <- c(paste(exposure_name,snp, 
                                 "Z_Z_Plot", sep = "_"), 
                           paste(exposure_name,snp, 
                                 "coloc", sep = "_"), 
                           paste(exposure_name,snp, 
                                 "Locus_Plot", sep = "_"))
    return_list <- c(return_list, return_dat)
  }
  
  return(return_list)
  
}


get_openGWAS_protein <- function(lead_SNP,protein_id, window = 75000 ) {
  ## get the SNP chromosome and position
  print("Retrieving SNP Position and Chromosome")
  snp_deets <- get_rsid_pos_simp(lead_SNP, window)
  
  ## Get SNP chr
  lead_SNP_chr <- gsub(pattern = "chr",replacement = "",x= snp_deets$Chr)
  
  ## Get SNP position
  lead_SNP_position <- snp_deets$Position
  
  ## Extract the SNPs associations with exposure protein for a cis window
  print("Grabbing out the associations within the specified window")
  exposure_window <- ieugwasr::associations(variants=paste(lead_SNP_chr, ":", 
                                                           lead_SNP_position-window, "-",
                                                           lead_SNP_position+window, sep = "" ), 
                                            id=c(protein_id), proxies=0)
  ## Get name of protein
  exposure_name <- unique(exposure_window$trait)
  ## convert to the twosample format
  temp_Format <- TwoSampleMR::format_data(exposure_window, 
                                          snp_col = "rsid" ,
                                          beta_col = "beta",
                                          se_col = "se",
                                          effect_allele_col = "ea" ,
                                          other_allele_col =  "nea",
                                          eaf_col = "eaf",
                                          pval_col = "p",
                                          pos_col = "position",
                                          chr_col = "chr",
                                          phenotype_col = "trait"
  )
  
  return(temp_Format)
  
}

coloc_multi_method <- function(Harm_dat, lead_SNP) {
  library(coloc)
  ## Grab out the LD
  LD <- LDlinkR::LDmatrix(Harm_dat$SNP, pop = "CEU", r2d = "r2", token = "4134a0a348dc", file = FALSE)
  
  ## We want to extract the column that is the lead SNP for plotting
  #LD_TEMP <- LD[,c("RS_number", lead_SNP)]
  
  
  ## Create Z scores for plotting
  Harm_dat$Z_exp <- Harm_dat$beta.exposure/Harm_dat$se.exposure
  Harm_dat$Z_out <- Harm_dat$beta.outcome/Harm_dat$se.outcome
  
  ## format the LD matrix
  LD <- LD[,-1]
  LD <- LD[,!colnames(LD) %in% names(which(colSums(is.na(as.matrix(LD))) > dim(LD)-1)) ]
  LD <- LD[complete.cases(LD),]
  rownames(LD) <- colnames(LD)
  ## restrict only to SNPs with LD values
  Harm_dat <- Harm_dat[Harm_dat$SNP %in% colnames(LD), ]
  Harm_dat <- Harm_dat <- Harm_dat[order(colnames(LD)),]
  
  ## Make sure that we have eaf for all SNPs
  Harm_dat <- Harm_dat[!is.na(Harm_dat$eaf.exposure),]
  
  
  ## Grab the position from the temp_dat_Format above
  
  #Harm_dat <- merge(Harm_dat, temp_dat[,c("name", "pos")], by.x = "SNP", by.y = "name", all.x = T)
  Harm_dat <- Harm_dat[!duplicated(Harm_dat$SNP),]
  D1=list(beta=Harm_dat$beta.exposure,
          varbeta=Harm_dat$se.exposure^2,
          snp=Harm_dat$SNP, 
          MAF=Harm_dat$eaf.exposure, 
          LD = as.matrix(LD),
          position = Harm_dat$pos.exposure,
          N=40000, 
          sdY=1,
          type="quant")
  
  D2=list(beta=Harm_dat$beta.outcome,
          varbeta=Harm_dat$se.outcome^2,
          snp=Harm_dat$SNP, 
          MAF=Harm_dat$eaf.exposure, 
          LD = as.matrix(LD),
          position = Harm_dat$pos.exposure,
          s=0.20, # This is the case proportion. Make sure to specify if using a case-control GWAS (get this information from the paper
          N=244230, #number of participants
          type="cc")
  
  
  colo_res_Cond_Iter <- coloc.signals(
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
  
  colo_res_mask_Iter <- coloc.signals(
    D1,
    D2,
    MAF = NULL,
    LD = NULL,
    method = c("mask"),
    mode = c("iterative"),
    p1 = 1e-03,
    p2 = 1e-04,
    p12 = 1e-05,
    maxhits = 3,
    r2thr = 0.01,
    pthr = 1e-06
  )
  colo_res_mask_Iter$summary$method <- "Mask_Iter"
  
  colo_res_mask_allbutone <-coloc.signals(
    D1,
    D2,
    MAF = NULL,
    LD = NULL,
    method = c("mask"),
    mode = c("allbutone"),
    p1 = 1e-03,
    p2 = 1e-04,
    p12 = 1e-05,
    maxhits = 3,
    r2thr = 0.01,
    pthr = 1e-06
  )
  colo_res_mask_allbutone$summary$method <- "Mask_All"
  
  
  colo_res <- coloc.abf(dataset1 = D1, 
                        dataset2 = D2,
                        p1 = 1e-03,
                        p2 = 1e-04,
                        p12 = 1e-05)
  
  coloc_res_format <- data.frame(hit2 = NA, hit1 = NA,
                                 nsnps = colo_res$summary[1],
                                 PP.H0.abf = colo_res$summary[2],
                                 PP.H1.abf = colo_res$summary[3],
                                 PP.H2.abf = colo_res$summary[4],
                                 PP.H3.abf = colo_res$summary[5],
                                 PP.H4.abf = colo_res$summary[6],
                                 best1 = NA, best2 = NA, best4 = NA,
                                 hit1.margz = NA, hit2.margz = NA,
                                 method = "Simple")
  coloc_res_All <- rbind(colo_res_Cond_Iter$summary,
                         colo_res_mask_Iter$summary,
                         colo_res_mask_allbutone$summary,coloc_res_format)
  coloc_res_All$protein <- unique(Harm_dat$exposure)
  coloc_res_All$outcome <- unique(Harm_dat$outcome)
  return(coloc_res_All)
}

coloc_public_prot_outcome <- function(SNP, exposure, outcome, window) {
  ## grab out the exposure information
  exposure_dat <- get_openGWAS_protein(SNP, exposure, window = window)
  ## grab out the outcome information
  outcome_dat <-TwoSampleMR::extract_outcome_data(exposure_dat$SNP, outcome)
  
  ## try to harmonise
  Harm_Dat <- TwoSampleMR::harmonise_data(
    exposure_dat,
    outcome_dat
  )
  
  temp_coloc_res <- coloc_multi_method(Harm_Dat, SNP)
  
}



get_peitzner_SNPs <- function(SNP, window, protein){
  library(dplyr)
  ## Get the SNP position
  temp_pos_chr <-get_rsid_pos_simp(SNP)
  
  ## We need a bash command
  command <- paste('tabix http://omicscience.org/apps/covidpgwas/data/all.grch37.tabix.gz ', gsub("chr","",temp_pos_chr[[1]]),":",
                   temp_pos_chr[[2]]-window,"-",temp_pos_chr[[2]]+window, sep = "")
  
  ## read in pQTL data for region of interest
  snp_dat_temp <- read.table(text = system(command,intern=TRUE))
  ## get column names
  colnames(snp_dat_temp) <- c("chr","pos","rsid","MarkerName","Allele1","Allele2","Freq1","FreqSE","MinFreq","MaxFreq",
                              "Somamer","Effect","StdErr","Pvalue","Direction","HetISq","HetChiSq","HetDf","HetPVal","TotalSampleSize")
  ## grab the somamer labels
  merge_file <- read.table("https://omicscience.org/apps/covidpgwas/data/somamer_label.txt", sep = "\t")
  ## name file nicely
  colnames(merge_file) <- c("Somamer","Somer_ID_2", "Short_Name", "Long_Name" )
  
  ## rename 
  merge_file$Short_Name <- recode_factor(merge_file$Short_Name, `gp130, soluble` = "IL-6RB")
  ## merge in the ids
  snp_dat_temp <- right_join(merge_file, snp_dat_temp)
  ## keep only the protein we are interested in
  return_dat <- snp_dat_temp %>% filter(Short_Name == protein)
  ## keep only those with rsids
  return_dat <- return_dat[return_dat$rsid %like% "rs",]
  
  
  return_dat <- TwoSampleMR::format_data(return_dat, type = "exposure",snp_col  = "rsid", 
                                         effect_allele_col ="Allele1", other_allele_col = "Allele2",
                                         eaf_col = "Freq1", beta_col = "Effect", se_col = "StdErr", 
                                         pval_col = "Pvalue", pos_col = "pos", chr_col = "chr", 
                                         phenotype_col = "Short_Name", samplesize_col = "TotalSampleSize")
  
  
  return(return_dat)
  
}
