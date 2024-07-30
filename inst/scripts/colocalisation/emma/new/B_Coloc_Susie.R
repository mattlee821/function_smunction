library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(coloc)
library(tidyverse)
library(plinkbinr)
library(biomaRt)
library(beepr)
library(Rfast)

setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/Year 1/Mini Project 2/Coloc")

overall_results<-data.frame(matrix(ncol=12,nrow=0))

x<-dir(pattern="exposure.csv")

#eQTLGen
eQTLx<-x[grep("eQTLGen",x)]

for (a in eQTLx){

snps1 <- read_exposure_data(
  filename = a,
  snp_col = "rsid",
  sep=",",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  eaf_col = "eaf"
)

snps2 <- read_outcome_data(
  snps = snps1$rsid,
  filename = paste(gsub("exposure","outcome",a)),
  snp_col = "SNP",
  sep=",",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "eaf.outcome"
)

snps <- harmonise_data(exposure_dat = snps1, outcome_dat = snps2)
snps$id.outcome<-"overall"

#LD matrix
ld<-ld_matrix_local(snps$SNP, with_alleles = FALSE, bfile="/Volumes/MRC-IEU-research/projects/icep2/wp2/034/working/data/1000GenomesReferenceFiles/EUR",plink_bin=get_plink_exe())
ld<-ld[which(rownames(ld) %in% snps$SNP),which(colnames(ld) %in% snps$SNP)]
snps<-snps[which(snps$SNP %in% rownames(ld)),]
ld<-ld[match(snps$SNP,rownames(ld)),]
ld<-ld[,match(snps$SNP,colnames(ld))]
snps<-snps[match(rownames(ld),snps$SNP),]


dataset1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = 31300,snps=rownames(ld),LD=ld)
dataset2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "cc", N = 98715,snps=rownames(ld),LD=ld)

df1<-runsusie(dataset1,coverage=0.1)
df2<-runsusie(dataset2,coverage=0.1)
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(df1,df2)
  print(susie.res$summary)
}  
df3 <- data.frame(susie.res$summary)
df3$gene<-substr(gsub("eqtl-a-","",a),1,15)
df3$site<-"overall"
overall_results<-rbind(overall_results,df3)
}

#BARCUVa
eQTLx<-x[grep("BARCUVa",x)]

for (a in eQTLx){
  
  snps1 <- read_exposure_data(
    filename = a,
    snp_col = "variant_id",
    sep=",",
    beta_col = "slope",
    se_col = "slope_se",
    effect_allele_col = "REF",
    other_allele_col = "ALT",
    eaf_col = "maf"
  )
  
  snps2 <- read_outcome_data(
    snps = snps1$rsid,
    filename = paste(gsub("exposure","outcome",a)),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    eaf_col = "eaf.outcome"
  )
  
  snps <- harmonise_data(exposure_dat = snps1, outcome_dat = snps2)
  snps$id.outcome<-"overall"
  
  #LD matrix
  ld<-ld_matrix_local(snps$SNP, with_alleles = FALSE, bfile="/Volumes/MRC-IEU-research/projects/icep2/wp2/034/working/data/1000GenomesReferenceFiles/EUR",plink_bin=get_plink_exe())
  ld<-ld[which(rownames(ld) %in% snps$SNP),which(colnames(ld) %in% snps$SNP)]
  snps<-snps[which(snps$SNP %in% rownames(ld)),]
  ld<-ld[match(snps$SNP,rownames(ld)),]
  ld<-ld[,match(snps$SNP,colnames(ld))]
  snps<-snps[match(rownames(ld),snps$SNP),]
  
  
  dataset1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = 31300,snps=rownames(ld),LD=ld)
  dataset2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "quant", N = 98715,snps=rownames(ld),LD=ld)
  
  df1<-runsusie(dataset1,coverage=0.1)
  df2<-runsusie(dataset2,coverage=0.1)
  if(requireNamespace("susieR",quietly=TRUE)) {
    susie.res=coloc.susie(df1,df2)
    print(susie.res$summary)
  }  
  df3 <- data.frame(susie.res$summary)
  df3$gene<-substr(gsub("eqtl-a-","",a),1,15)
  df3$site<-"overall"
  overall_results<-rbind(overall_results,df3)
}


write.csv(overall_results,"Overall_eQTLGen_Coloc_Results.csv",row.names = FALSE,quote=FALSE)
