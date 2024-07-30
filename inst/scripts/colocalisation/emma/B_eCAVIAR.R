sink("caviardistal.txt")

#Create Z-scores
library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(coloc)
library(tidyverse)
#devtools::install_github("explodecomputer/plinkbinr")
library(plinkbinr)
setwd("~/MiniProject2/genexpression/Coloc/caviar")

#Get a list of analyses needed to perform
df<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/034/working/results/MiniProject2/MR/MR_Results.csv")
df<-df[which(df$bh_p<0.05),]

#eQTLGen
eQTL<-df[which(df$type=="eQTLGen"),]

#distal
distal_eQTL<-eQTL[which(eQTL$id.outcome=="distal"),]
exposure_dat<-extract_instruments(outcomes=paste("eqtl-a-",distal_eQTL$ensembl_gene_id,sep=""), p1=1e-06)

file <-fread("/projects/MRC_IEU/research/projects/icep2/wp2/034/working/results/MiniProject2/01_Gene Expression/Allgeneswithgencode_AN.csv")
exposure_dat$gene<-substr(exposure_dat$id.exposure,8,22)
file$gene<-file$ensembl_gene_id
file$Chromosome<-as.numeric(file$ensembl_chr)
file$cisstart<-as.numeric(file$ensembl_cisstart)
file$cisend<-as.numeric(file$ensembl_cisend)
file<-select(file,gene,Chromosome,cisstart,cisend)
exposure_dat<-merge(exposure_dat,file,by="gene")
exposure_dat<-exposure_dat[which((exposure_dat$chr.exposure=exposure_dat$Chromosome) & (exposure_dat$pos.exposure <=exposure_dat$cisend)& (exposure_dat$pos.exposure >=exposure_dat$cisstart)),]

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

dat <- dat[match(unique(dat$gene), dat$gene),]

for (i in 1:nrow(dat)){
  print(i)
  x<-dat[i,]
  x$id.exposure<-substr(x$id.exposure,8,22)
  print(x$id.exposure)
  print(x$SNP)
  x$window_start<-x$pos.exposure-1000000
  x$window_end<-x$pos.exposure+1000000
  variants<-paste(x$chr.exposure,x$window_start,sep=":")
  variants<-paste(variants,x$window_end,sep="-")
  snps<-associations(id=paste("eqtl-a-",x$id.exposure,sep=""),variants=variants)
  outcome_dat<- read_outcome_data(
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
  fwrite(snps,"snps.csv",row.names = FALSE,quote=FALSE)
  snps <- read_exposure_data(
    filename = "snps.csv",
    snp_col = "rsid",
    sep=",",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "ea",
    other_allele_col = "nea",
    eaf_col = "eaf"
  )
  
  snps$Z.eQTL<-snps$beta.exposure/snps$se.exposure
  snps<-select(snps,SNP,Z.eQTL)
  outcome_dat$Z.GWAS<-outcome_dat$beta.outcome/outcome_dat$se.outcome
  outcome_dat<-select(outcome_dat,SNP,Z.GWAS)
  snps<-na.omit(snps)
  outcome_dat<-na.omit(outcome_dat)
  snps<-merge(snps,outcome_dat,by="SNP")
  ld<-ld_matrix_local(snps$SNP, with_alleles = FALSE, bfile="/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/1000GenomesReferenceFiles/EUR",plink_bin=get_plink_exe())
  ld<-ld[which(rownames(ld) %in% snps$SNP),which(colnames(ld) %in% snps$SNP)]
  snps<-snps[which(snps$SNP %in% rownames(ld)),]
  ld<-ld[match(snps$SNP,rownames(ld)),]
  ld<-ld[,match(snps$SNP,colnames(ld))]
  snps<-snps[match(rownames(ld),snps$SNP),]
  snps1<-select(snps,SNP,Z.eQTL)
  snps2<-select(snps,SNP,Z.GWAS)
  write.table(snps1,paste("eQTLGen_",x$id.exposure,"_distal_SNPs_Z_eQTL",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(snps2,paste("eQTLGen_",x$id.exposure,"_distal_SNPs_Z_GWAS",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(ld,paste("eQTLGen_",x$id.exposure,"_distal_SNPs_LD",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  
}

#BARCeQTL
eQTL<-df[which(df$type=="BARCeQTL"),]
x<-dir("/home/sw20203/MiniProject2/genexpression/BARC-UVA",pattern="clumped.csv")
z<-substr(x,1,15)
colon<-fread("/home/sw20203/MiniProject2/genexpression/BARC-UVA/barcuvaseq.eqtls.allpairs.sumstats.txt")
colon$id.exposure<-substr(colon$gene_id,1,15)

#distal
CRC<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
distal_eQTL<-eQTL[which(eQTL$id.outcome=="distal"),]
z1<-intersect(z,distal_eQTL$ensembl_gene_id)
if (length(z1)>0){
distal_results<-data.frame(matrix(nrow=0,ncol=9))

for (i in 1:length(z1)){
  y<-z1[i]
  print(i)
  print(y)
  
  exposure_dat <- read_exposure_data(
    filename = paste("/home/sw20203/MiniProject2/genexpression/BARC-UVA/",y,"_BARC_eqtl_clumped.csv",sep=""),
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
  snps<-colon[which(colon$id.exposure==y & colon$chr==p$chr.exposure & colon$pos>=p$window_start & colon$pos<=p$window_end),]
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
  fwrite(snps,"snps.csv",row.names = FALSE,quote=FALSE)
  snps <- read_exposure_data(
    filename = "snps.csv",
    snp_col = "variant_id",
    sep=",",
    beta_col = "slope",
    se_col = "slope_se",
    effect_allele_col = "REF",
    other_allele_col = "ALT",
    eaf_col = "maf"
  )
  snps$Z.eQTL<-snps$beta.exposure/snps$se.exposure
  snps<-select(snps,SNP,Z.eQTL)
  outcome_dat$Z.GWAS<-outcome_dat$beta.outcome/outcome_dat$se.outcome
  outcome_dat<-select(outcome_dat,SNP,Z.GWAS)
  snps<-na.omit(snps)
  outcome_dat<-na.omit(outcome_dat)
  snps<-merge(snps,outcome_dat,by="SNP")
  ld<-ld_matrix_local(snps$SNP, with_alleles = FALSE, bfile="/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/1000GenomesReferenceFiles/EUR",plink_bin=get_plink_exe())
  ld<-ld[which(rownames(ld) %in% snps$SNP),which(colnames(ld) %in% snps$SNP)]
  snps<-snps[which(snps$SNP %in% rownames(ld)),]
  ld<-ld[match(snps$SNP,rownames(ld)),]
  ld<-ld[,match(snps$SNP,colnames(ld))]
  snps<-snps[match(rownames(ld),snps$SNP),]
  snps1<-select(snps,SNP,Z.eQTL)
  snps2<-select(snps,SNP,Z.GWAS)
  write.table(snps1,paste("BARC_",y,"_distal_SNPs_Z_eQTL",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(snps2,paste("BARC_",y,"_distal_SNPs_Z_GWAS",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(ld,paste("BARC_",y,"_distal_SNPs_LD",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  
}}

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
z1<-intersect(z,distal_eQTL$ensembl_gene_id)
if (length(z1)>0){
distal_results<-data.frame(matrix(nrow=0,ncol=9))

for (i in 1:length(z1)){
  y<-z1[i]
  print(i)
  print(z1[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("/home/sw20203/MiniProject2/genexpression/exposuredat/p0.000005/",z1[i],"_colontotal_eqtl_clumped.csv",sep=""),
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
  fwrite(snps,"snps.csv",row.names = FALSE,quote=FALSE)
  snps <- read_exposure_data(
    filename = "snps.csv",
    snp_col = "rs_id_dbSNP151_GRCh38p7",
    sep=",",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "ref",
    other_allele_col = "alt",
    eaf_col = "maf"
  )
  snps$Z.eQTL<-snps$beta.exposure/snps$se.exposure
  snps<-select(snps,SNP,Z.eQTL)
  outcome_dat$Z.GWAS<-outcome_dat$beta.outcome/outcome_dat$se.outcome
  outcome_dat<-select(outcome_dat,SNP,Z.GWAS)
  snps<-na.omit(snps)
  outcome_dat<-na.omit(outcome_dat)
  snps<-merge(snps,outcome_dat,by="SNP")
  ld<-ld_matrix_local(snps$SNP, with_alleles = FALSE, bfile="/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/1000GenomesReferenceFiles/EUR",plink_bin=get_plink_exe())
  ld<-ld[which(rownames(ld) %in% snps$SNP),which(colnames(ld) %in% snps$SNP)]
  snps<-snps[which(snps$SNP %in% rownames(ld)),]
  ld<-ld[match(snps$SNP,rownames(ld)),]
  ld<-ld[,match(snps$SNP,colnames(ld))]
  snps<-snps[match(rownames(ld),snps$SNP),]
  snps1<-select(snps,SNP,Z.eQTL)
  snps2<-select(snps,SNP,Z.GWAS)
  write.table(snps1,paste("GTExTotal_",y,"_distal_SNPs_Z_eQTL",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(snps2,paste("GTExTotal_",y,"_distal_SNPs_Z_GWAS",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(ld,paste("GTExTotal_",y,"_distal_SNPs_LD",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
}
}

#Sigmoid
eQTL<-df[which(df$type=="sigmoideQTL"),]
x<-dir("/home/sw20203/MiniProject2/genexpression/GTExMetaAnalysis",pattern=eQTL$id.exposure)
x<-x[grep("clumped",x)]
z<-substr(x,1,15)
colon<-fread("/home/sw20203/MiniProject2/genexpression/GTExMetaAnalysis/GTExMetaAnalysisResultsFixed.csv")
colon$id.exposure<-substr(colon$phenotype_id,1,15)

#distal
CRC<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
distal_eQTL<-eQTL[which(eQTL$id.outcome=="distal"),]
z1<-intersect(z,distal_eQTL$ensembl_gene_id)
if (length(z1)>0){
distal_results<-data.frame(matrix(nrow=0,ncol=9))

for (i in 1:length(z1)){
  y<-z1[i]
  print(i)
  print(z1[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("/home/sw20203/MiniProject2/genexpression/exposuredat/p0.000005/",z1[i],"_colonsigmoid_eqtl_clumped.csv",sep=""),
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
  fwrite(snps,"snps.csv",row.names = FALSE,quote=FALSE)
  snps <- read_exposure_data(
    filename = "snps.csv",
    snp_col = "rs_id_dbSNP151_GRCh38p7",
    sep=",",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "ref",
    other_allele_col = "alt",
    eaf_col = "maf"
  )
  snps$Z.eQTL<-snps$beta.exposure/snps$se.exposure
  snps<-select(snps,SNP,Z.eQTL)
  outcome_dat$Z.GWAS<-outcome_dat$beta.outcome/outcome_dat$se.outcome
  outcome_dat<-select(outcome_dat,SNP,Z.GWAS)
  snps<-na.omit(snps)
  outcome_dat<-na.omit(outcome_dat)
  snps<-merge(snps,outcome_dat,by="SNP")
  ld<-ld_matrix_local(snps$SNP, with_alleles = FALSE, bfile="/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/1000GenomesReferenceFiles/EUR",plink_bin=get_plink_exe())
  ld<-ld[which(rownames(ld) %in% snps$SNP),which(colnames(ld) %in% snps$SNP)]
  snps<-snps[which(snps$SNP %in% rownames(ld)),]
  ld<-ld[match(snps$SNP,rownames(ld)),]
  ld<-ld[,match(snps$SNP,colnames(ld))]
  snps<-snps[match(rownames(ld),snps$SNP),]
  snps1<-select(snps,SNP,Z.eQTL)
  snps2<-select(snps,SNP,Z.GWAS)
  write.table(snps1,paste("GTExSigmoid_",y,"_distal_SNPs_Z_eQTL",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(snps2,paste("GTExSimgoid_",y,"_distal_SNPs_Z_GWAS",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(ld,paste("GTExSigmoid_",y,"_distal_SNPs_LD",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
}
}

#Transverse
eQTL<-df[which(df$type=="transverseeQTL"),]
x<-dir("/home/sw20203/MiniProject2/genexpression/GTExMetaAnalysis",pattern=eQTL$id.exposure)
x<-x[grep("clumped",x)]
z<-substr(x,1,15)
colon<-fread("/home/sw20203/MiniProject2/genexpression/GTExMetaAnalysis/GTExMetaAnalysisResultsFixed.csv")
colon$id.exposure<-substr(colon$phenotype_id,1,15)

#distal
CRC<-fread("/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
distal_eQTL<-eQTL[which(eQTL$id.outcome=="distal"),]
z1<-intersect(z,distal_eQTL$ensembl_gene_id)
if (length(z1)>0){
distal_results<-data.frame(matrix(nrow=0,ncol=9))

for (i in 1:length(z1)){
  y<-z1[i]
  print(i)
  print(z1[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("/home/sw20203/MiniProject2/genexpression/exposuredat/p0.000005/",z1[i],"_colontransverse_eqtl_clumped.csv",sep=""),
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
  fwrite(snps,"snps.csv",row.names = FALSE,quote=FALSE)
  snps <- read_exposure_data(
    filename = "snps.csv",
    snp_col = "rs_id_dbSNP151_GRCh38p7",
    sep=",",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "ref",
    other_allele_col = "alt",
    eaf_col = "maf"
  )
  snps$Z.eQTL<-snps$beta.exposure/snps$se.exposure
  snps<-select(snps,SNP,Z.eQTL)
  outcome_dat$Z.GWAS<-outcome_dat$beta.outcome/outcome_dat$se.outcome
  outcome_dat<-select(outcome_dat,SNP,Z.GWAS)
  snps<-na.omit(snps)
  outcome_dat<-na.omit(outcome_dat)
  snps<-merge(snps,outcome_dat,by="SNP")
  ld<-ld_matrix_local(snps$SNP, with_alleles = FALSE, bfile="/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/1000GenomesReferenceFiles/EUR",plink_bin=get_plink_exe())
  ld<-ld[which(rownames(ld) %in% snps$SNP),which(colnames(ld) %in% snps$SNP)]
  snps<-snps[which(snps$SNP %in% rownames(ld)),]
  ld<-ld[match(snps$SNP,rownames(ld)),]
  ld<-ld[,match(snps$SNP,colnames(ld))]
  snps<-snps[match(rownames(ld),snps$SNP),]
  snps1<-select(snps,SNP,Z.eQTL)
  snps2<-select(snps,SNP,Z.GWAS)
  write.table(snps1,paste("GTExTransverse_",y,"_distal_SNPs_Z_eQTL",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(snps2,paste("GTExTransverse_",y,"_distal_SNPs_Z_GWAS",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
  write.table(ld,paste("GTExTransverse_",y,"_distal_SNPs_LD",".txt",sep=""),row.names = FALSE,quote=FALSE,sep="\t",col.names=FALSE)
}
}

sink()