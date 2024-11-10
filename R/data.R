#' mapping_BUILD_PHASE_somalogic
#'
#' this is a mapping file for somalgic IDs that gives IDs and genomic information
#' the GRC build in the file name is the build used to map somalogic IDs to other IDs
#' the code used to create this file is [here](https://github.com/mattlee821/functions/blob/main/inst/scripts/mapping-somalogic.R)
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{SeqId}{sequence ID from somalogic}
#'   \item{SomaId}{somalogic ID from somalogic}
#'   \item{UNIPROT}{UNIPROT ID from somalogic}
#'   \item{Target}{protein short name from somalogic}
#'   \item{TargetFullName}{protein long name from somalogic}
#'   \item{EntrezGeneID}{ENTREZ gene ID from somalogic}
#'   \item{EntrezGeneSymbol}{ENTREZ gene symbol from somalogic}
#'   \item{uniprot_gn_id}{uniprot ID mapped by me from biomart using ensembl}
#'   \item{uniprot_gn_symbol}{uniprot symbol mapped by me from biomart using ensembl}
#'   \item{entrezgene_id}{entrez gene ID mapped by me from biomart using ensembl}
#'   \item{entrezgene_accession}{entrez gene symbol mapped by me from biomart using ensembl}
#'   \item{hgnc_id}{hgnc gene ID mapped by me from biomart using ensembl}
#'   \item{hgnc_symbol}{hgnc gene symbol ID mapped by me from biomart using ensembl}
#'   \item{external_gene_name}{external gene name (unsure really as description from biomart isnt good) mapped by me from biomart using ensembl}
#'   \item{gene_biotype}{gene type (restricted to protein-coding) mapped by me from biomart using ensembl}
#'   \item{CHR}{chromosome the gene is on mapped by me from TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   \item{START_hg19}{gene start position on chromosome ID mapped by me from TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   \item{END_hg19}{gene end position on chromosome ID mapped by me from TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   \item{strand_hg19}{strand (positive or negative) for start and end position mapped by me from TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   \item{START_hg38}{gene start position on chromosome mapped by me from TxDb.Hsapiens.UCSC.hg38.knownGene}
#'   \item{END_hg38}{gene end position on chromosome mapped by me from TxDb.Hsapiens.UCSC.hg38.knownGene}
#'   \item{strand_hg38}{strand (positive or negative) for start and end position mapped by me from TxDb.Hsapiens.UCSC.hg38.knownGene}
#' }
#' @name mapping_GRCh38_p14_somalogic
"mapping_GRCh38_p14_somalogic"

#' mapping_BUILD_PHASE_olink
#'
#' this is a mapping file for somalgic IDs that gives IDs and genomic information
#' the GRC build in the file name is the build used to map somalogic IDs to other IDs
#' the code used to create this file is [here](https://github.com/mattlee821/functions/blob/main/inst/scripts/mapping-somalogic.R)
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{olinkID}{sequence ID from olink}
#'   \item{UNIPROT}{UNIPROT ID from olink}
#'   \item{Target}{protein short name from olink}
#'   \item{TargetFullName}{protein long name from olink}
#'   \item{uniprot_gn_id}{uniprot ID mapped by me from biomart using ensembl}
#'   \item{uniprot_gn_symbol}{uniprot symbol mapped by me from biomart using ensembl}
#'   \item{entrezgene_id}{entrez gene ID mapped by me from biomart using ensembl}
#'   \item{entrezgene_accession}{entrez gene symbol mapped by me from biomart using ensembl}
#'   \item{hgnc_id}{hgnc gene ID mapped by me from biomart using ensembl}
#'   \item{hgnc_symbol}{hgnc gene symbol ID mapped by me from biomart using ensembl}
#'   \item{external_gene_name}{external gene name (unsure really as description from biomart isnt good) mapped by me from biomart using ensembl}
#'   \item{gene_biotype}{gene type (restricted to protein-coding) mapped by me from biomart using ensembl}
#'   \item{CHR}{chromosome the gene is on mapped by me from TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   \item{START_hg19}{gene start position on chromosome ID mapped by me from TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   \item{END_hg19}{gene end position on chromosome ID mapped by me from TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   \item{strand_hg19}{strand (positive or negative) for start and end position mapped by me from TxDb.Hsapiens.UCSC.hg19.knownGene}
#'   \item{START_hg38}{gene start position on chromosome mapped by me from TxDb.Hsapiens.UCSC.hg38.knownGene}
#'   \item{END_hg38}{gene end position on chromosome mapped by me from TxDb.Hsapiens.UCSC.hg38.knownGene}
#'   \item{strand_hg38}{strand (positive or negative) for start and end position mapped by me from TxDb.Hsapiens.UCSC.hg38.knownGene}
#' }
#' @name mapping_GRCh38_p14_olink
"mapping_GRCh38_p14_olink"
