sink("Coloc.txt")

library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(coloc)
library(tidyverse)
library(plinkbinr)
library(biomaRt)
setwd("/home/sw20203/MiniProject2/genexpression/Coloc")

distal_samplesize<-98715
proximal_samplesize<-57515
distal_samplesize<-55978
rectal_samplesize<-57249
left_samplesize<-70103
colon_samplesize<-71835

#Get a list of analyses needed to perform
df<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/034/working/results/MiniProject2/MR/MR_Results.csv")

df<-df[which(df$bh_p<0.05),]

#eQTLGen
eQTL<-df[which(df$type=="eQTLGen"),]

#distal
distal_eQTL<-eQTL[which(eQTL$id.outcome=="distal"),]
exposure_dat<-extract_instruments(outcomes=paste("eqtl-a-",distal_eQTL$id.exposure,sep=""), p1=1e-06)

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
gene<-t(data.frame(unique(lapply(exposure_dat$id.exposure,gsub,pattern="eqtl-a-",replacement=""))))
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","start_position","chromosome_name"), values=gene, mart=ensembl,useCache=FALSE)
genes$cisstart <- genes$start_position-1000000
genes$cisend <- genes$start_position+1000000

exposure_dat$ensembl_gene_id<-as.character(lapply(exposure_dat$id.exposure,gsub,pattern="eqtl-a-",replacement=""))
exposure_dat<-merge(exposure_dat,genes,by="ensembl_gene_id")
exposure_dat<-exposure_dat[exposure_dat$chr.exposure==exposure_dat$chromosome_name & exposure_dat$pos.exposure>=exposure_dat$cisstart,]
exposure_dat<-exposure_dat[exposure_dat$chr.exposure==exposure_dat$chromosome_name & exposure_dat$pos.exposure<=exposure_dat$cisend,]

outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

dat<-unique(setDT(dat)[order(pval.exposure)], by = "id.exposure")


#Keep only the most significant SNP for each gene
dat <- dat[match(unique(dat$id.exposure), dat$id.exposure),]

dat$window_start<-dat$pos.exposure-1000000
dat$window_end<-dat$pos.exposure+1000000

for (i in 1:nrow(dat)){
  print(i)
  x<-dat[i,]
  print(x$id.exposure)
  print(x$SNP)
  variants<-paste(x$chr.exposure,x$window_start,sep=":")
  variants<-paste(variants,x$window_end,sep="-")
  snps<-associations(id=x$id.exposure,variants=variants)
  outcome_dat <- read_outcome_data(
    snps = snps$rsid,
    filename = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  fwrite(snps,paste(x$id.exposure,"_distal_","eQTLGen_","exposure.csv",sep=""),row.names = FALSE,quote=FALSE)
  fwrite(outcome_dat,paste(x$id.exposure,"_distal_","eQTLGen_","outcome.csv",sep=""),row.names = FALSE,quote=FALSE)
  }

#BARCeQTL
eQTL<-df[which(df$type=="BARCeQTL"),]
colon<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/BARCUVa-Seq/barcuvaseq.eqtls.allpairs.sumstats.txt")


#distal
CRC<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
distal_eQTL<-eQTL[which(eQTL$id.outcome=="distal"),]

for (i in 1:length(unique(distal_eQTL$id.exposure))){
  y<-unique(distal_eQTL$id.exposure)[i]
  print(i)
  print(y)
  
  exposure_dat <- read_exposure_data(
    filename = paste("/home/sw20203/MiniProject2/genexpression/MRResults/",y,"_BARC_eqtl_clumped.csv",sep=""),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure"
  )
  
  snps_for_proxy <- setdiff(exposure_dat$SNP,CRC$SNP)
  
  if (length(snps_for_proxy)!=0){
    proxySNPs=data.frame(matrix(ncol=4,nrow=0))
    for(j in snps_for_proxy) {
      print(j)
      proxy <- LDproxy(j, pop = "EUR", r2d = "r2", token = "a563a6f57fcc")
      proxy2 <- proxy %>% dplyr::filter(R2 > 0.8) %>% # minimum R2 0.8
        dplyr::filter(RS_Number != ".") %>% # remove ones without rsID
        dplyr::filter(RS_Number %in% CRC$SNP) # only keep SNPs that are in CRC data
      if (length(proxy2)!=0){
        proxy3 <- proxy2[1,] # keep rsID at the top with the highest LD
        proxy3$missing_snp <- j # label the SNP we were trying to proxy so we can merge back in when creating summary dat
        ea<-exposure_dat[which(exposure_dat$SNP==j),5]
        oa<-exposure_dat[which(exposure_dat$SNP==j),6]
        Correlated_Alleles<-data.frame(proxy3$Correlated_Alleles)
        Correlated_Alleles$proxy3.Correlated_Alleles<- as.character(Correlated_Alleles$proxy3.Correlated_Alleles)
        Correlated_Alleles<-Correlated_Alleles %>% separate(proxy3.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
        Correlated_Alleles<-na.omit(Correlated_Alleles)
        if (nrow(Correlated_Alleles)!=0){
          proxy4<-cbind(proxy3,Correlated_Alleles)
          proxy4<-cbind(proxy4,ea)
          proxy4<-cbind(proxy4,oa)
          
          if (proxy4$a1exp==proxy4$ea){
            proxy4$REF<-proxy4$a1out
            proxy4$ALT<-proxy4$a2out
          } else if(proxy4$a2exp==proxy4$ea){
            proxy4$REF<-proxy4$a2out
            proxy4$ALT<-proxy4$a1out
          } else {
            proxy4$REF<-"NA"
            proxy4$ALT<-"NA"
          }
          proxy5<-select(proxy4,missing_snp,RS_Number,REF,ALT)
          proxySNPs<-rbind(proxySNPs,proxy5)
        }
      } else{
        proxySNPs<-data.frame(matrix(ncol=17,nrow=0))
        colnames(proxySNPs)<-"missing_snp" 
      }
      
      
      
      exposure_datproxies<-exposure_dat[which(exposure_dat$SNP %in% snps_for_proxy),]
      
      if (nrow(proxySNPs)!=0){
        
        exposure_dat2<-merge(exposure_datproxies,proxySNPs,by.x="SNP",by.y="missing_snp")
        exposure_dat2$SNP<-exposure_dat2$RS_Number
        exposure_dat2$effect_allele.exposure<-exposure_dat2$REF
        exposure_dat2$other_allele.exposure<-exposure_dat2$ALT
        exposure_dat2<-exposure_dat2[,1:12]
        exposure_dat<-rbind(exposure_dat,exposure_dat2)
      }}}
  
  
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$rsid,
    filename = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  
  
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
  
  p<-unique(setDT(dat)[order(pval.exposure)])
  p$pos.exposure<-CRC$chrom_start[which(CRC$SNP==p$SNP)]
  p$window_start<-p$pos.exposure-1000000
  p$window_end<-p$pos.exposure+1000000
  p$chr.exposure<-CRC$chr_name[which(CRC$SNP==p$SNP)]
  snps<-colon[which(colon$chr==p$chr.exposure & colon$pos>=p$window_start & colon$pos<=p$window_end),]
  snps$id.exposure<-substr(snps$gene_id,1,15)
  snps<-snps[which(snps$id.exposure==y),]
  p$id.exposure<-y
  outcome_dat <- read_outcome_data(
    snps = snps$variant_id,
    filename = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  fwrite(snps,paste(y,"_distal_","BARCUVa_","exposure.csv",sep=""),row.names = FALSE,quote=FALSE)
  fwrite(outcome_dat,paste(y,"_distal_","BARCUVa_","outcome.csv",sep=""),row.names = FALSE,quote=FALSE)
}

#GTEx
eQTL<-df[which(df$type=="totaleQTL"),]
x<-dir("/home/sw20203/MiniProject2/genexpression/GTExMetaAnalysis",pattern=eQTL$id.exposure)
x<-x[grep("clumped",x)]
z<-substr(x,1,15)
colon<-fread("/home/sw20203/MiniProject2/genexpression/GTExMetaAnalysis/GTExMetaAnalysisResultsFixed.csv")
colon$id.exposure<-substr(colon$phenotype_id,1,15)

#distal
CRC<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
distal_eQTL<-eQTL[which(eQTL$id.outcome=="distal"),]
z1<-intersect(z,distal_eQTL$id.exposure)

for (i in 1:length(x)){
  y<-z1[i]
  print(i)
  print(z1[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("/home/sw20203/MiniProject2/genexpression/GTExMetaAnalysis/",x[i],sep=""),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure"
  )
  
  snps_for_proxy <- setdiff(exposure_dat$SNP,CRC$SNP)
  
  if (length(snps_for_proxy)!=0){
    proxySNPs=data.frame(matrix(ncol=4,nrow=0))
    for(j in snps_for_proxy) {
      print(j)
      proxy <- LDproxy(j, pop = "EUR", r2d = "r2", token = "a563a6f57fcc")
      proxy2 <- proxy %>% dplyr::filter(R2 > 0.8) %>% # minimum R2 0.8
        dplyr::filter(RS_Number != ".") %>% # remove ones without rsID
        dplyr::filter(RS_Number %in% CRC$SNP) # only keep SNPs that are in CRC data
      if (length(proxy2)!=0){
        proxy3 <- proxy2[1,] # keep rsID at the top with the highest LD
        proxy3$missing_snp <- j # label the SNP we were trying to proxy so we can merge back in when creating summary dat
        ea<-exposure_dat[which(exposure_dat$SNP==j),5]
        oa<-exposure_dat[which(exposure_dat$SNP==j),6]
        Correlated_Alleles<-data.frame(proxy3$Correlated_Alleles)
        Correlated_Alleles$proxy3.Correlated_Alleles<- as.character(Correlated_Alleles$proxy3.Correlated_Alleles)
        Correlated_Alleles<-Correlated_Alleles %>% separate(proxy3.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
        Correlated_Alleles<-na.omit(Correlated_Alleles)
        if (nrow(Correlated_Alleles)!=0){
          proxy4<-cbind(proxy3,Correlated_Alleles)
          proxy4<-cbind(proxy4,ea)
          proxy4<-cbind(proxy4,oa)
          
          if (proxy4$a1exp==proxy4$ea){
            proxy4$REF<-proxy4$a1out
            proxy4$ALT<-proxy4$a2out
          } else if(proxy4$a2exp==proxy4$ea){
            proxy4$REF<-proxy4$a2out
            proxy4$ALT<-proxy4$a1out
          } else {
            proxy4$REF<-"NA"
            proxy4$ALT<-"NA"
          }
          proxy5<-select(proxy4,missing_snp,RS_Number,REF,ALT)
          proxySNPs<-rbind(proxySNPs,proxy5)
        }
      } else{
        proxySNPs<-data.frame(matrix(ncol=17,nrow=0))
        colnames(proxySNPs)<-"missing_snp" 
      }
      
      
      
      exposure_datproxies<-exposure_dat[which(exposure_dat$SNP %in% snps_for_proxy),]
      
      if (nrow(proxySNPs)!=0){
        
        exposure_dat2<-merge(exposure_datproxies,proxySNPs,by.x="SNP",by.y="missing_snp")
        exposure_dat2$SNP<-exposure_dat2$RS_Number
        exposure_dat2$effect_allele.exposure<-exposure_dat2$REF
        exposure_dat2$other_allele.exposure<-exposure_dat2$ALT
        exposure_dat2<-exposure_dat2[,1:12]
        exposure_dat<-rbind(exposure_dat,exposure_dat2)
      }}}
  
  
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$rsid,
    filename = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  
  
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
  
  
  
  p<-unique(setDT(dat)[order(pval.exposure)], by = "id.exposure")
  p$pos.exposure<-CRC$chrom_start[which(CRC$SNP==p$SNP)]
  p$window_start<-p$pos.exposure-1000000
  p$window_end<-p$pos.exposure+1000000
  p$chr.exposure<-CRC$chr_name[which(CRC$SNP==p$SNP)]
  snps<-colon[which(colon$id.exposure==y),]
  colon_snp<-snps %>% separate(variant_id, c("CHR", "BP", "other_allele.exposure", "effect_allele.exposure"))
  colon_snp$chr<- as.numeric(gsub("chr", "", colon_snp$CHR)) #remove the "chr" letters
  snps<-colon_snp[which(colon_snp$id.exposure==y & colon_snp$chr==p$chr.exposure & colon_snp$BP>=p$window_start & colon_snp$BP<=p$window_end),]
  
  #Get rsIDs
  snp_info<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
  snps$BP<-as.numeric(snps$BP)
  snps<-merge(snps,snp_info,by.x=c("CHR","BP"),by.y=c("chr","variant_pos"))
  
  p$id.exposure<-y
  outcome_dat <- read_outcome_data(
    snps = snps$rs_id_dbSNP151_GRCh38p7,
    filename = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  fwrite(snps,paste(y,"_distal_","GTEx_","exposure.csv",sep=""),row.names = FALSE,quote=FALSE)
  fwrite(outcome_dat,paste(y,"_distal_","GTEx_","outcome.csv",sep=""),row.names = FALSE,quote=FALSE)
}

sink()