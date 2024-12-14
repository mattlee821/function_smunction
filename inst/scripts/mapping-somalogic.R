rm(list=ls())
set.seed(821)

# environment ====
# BiocManager::install(version = "3.20", force = T)
# BiocManager::version()
# BiocManager::valid()
# BiocManager::install(c("Biobase",
#                        "AnnotationDbi",
#                        "org.Hs.eg.db",
#                        "GenomicFeatures",
#                        "TxDb.Hsapiens.UCSC.hg19.knownGene",
#                        "TxDb.Hsapiens.UCSC.hg38.knownGene",
#                        "biomaRt",
#                        "UniProt.ws"),
#                      update = TRUE, ask = FALSE, force = TRUE)

library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)
library(UniProt.ws)
library(SomaDataIO)
library(dplyr)
library(tidyr)

# source ====
biomaRt_getBM_batch <- function(mart, attributes, filters, values, chunk_size) {
  # Remove `filters` from `attributes` if it's already present
  if (filters %in% attributes) {
    attributes <- setdiff(attributes, filters)
  }

  # Split attributes into chunks, each including the filters column
  attribute_chunks <- lapply(
    split(attributes, ceiling(seq_along(attributes) / chunk_size)),
    function(chunk) c(filters, chunk)
  )

  # Loop over each chunk and fetch data, storing results in a list
  result_list <- list()
  for (i in seq_along(attribute_chunks)) {
    result <- tryCatch(
      getBM(
        attributes = attribute_chunks[[i]],
        filters = filters,
        values = values,
        mart = mart
      ),
      error = function(e) NULL # If a query fails, return NULL for this chunk
    )

    # Only add result to the list if it's not NULL and contains the identifier column
    if (!is.null(result) && filters %in% colnames(result)) {
      result_list[[i]] <- result
    }
  }

  # Sequentially join each data frame in result_list
  sequential_join <- function(result_list, filters, values) {
    # Start with a data frame containing only the 'values' column
    final_result <- data.frame(setNames(list(values), filters))

    # Sequentially join each data frame in result_list
    for (i in seq_along(result_list)) {
      final_result <- full_join(final_result, result_list[[i]], by = filters)

      # Free up memory after each join
      gc()
    }

    return(final_result)
  }

  # Perform sequential join on result_list
  map_uniprot <- sequential_join(result_list, filters = filters, values = values)

  # Apply filtering on the resulting data frame
  map_uniprot <- map_uniprot %>%
    dplyr::filter(gene_biotype == "protein_coding") %>%
    dplyr::filter(chromosome_name %in% c(as.character(1:22), "X", "Y", "")) %>%
    dplyr::filter(uniprot_gn_symbol == entrezgene_accession |
                    uniprot_gn_symbol == hgnc_symbol |
                    uniprot_gn_symbol == external_gene_name) %>%
    dplyr::rename(CHR = chromosome_name) %>%
    distinct()

  return(map_uniprot)
}

# data ====
# wget https://github.com/SomaLogic/SomaLogic-Data/raw/main/example_data.adat # Retrieve just the 5k (v4.0) ADAT
# wget https://github.com/SomaLogic/SomaLogic-Data/raw/main/example_data_v4.1_plasma.adat # Retrieve just the 7k (v4.1) ADAT
# wget https://github.com/SomaLogic/SomaLogic-Data/raw/main/example_data_v5.0_plasma.adat # Retrieve just the 11k (v5.0) ADAT
data <- SomaDataIO::read_adat("inst/data/example_data_v5.0_plasma.adat")
data <- SomaDataIO::getAnalyteInfo(data)

## format data
data <- data %>%
  tidyr::separate_rows(UniProt, sep = "\\|") %>%
  tidyr::separate_rows(EntrezGeneID, sep = "\\|") %>%
  tidyr::separate_rows(EntrezGeneSymbol, sep = "\\|") %>%
  dplyr::filter(Organism == "Human",
                Type == "Protein", # remove controls
                !if_any(everything(), ~ grepl("internal use only", ., ignore.case = TRUE))) %>% # remove internal use only
  dplyr::select(SeqId, SomaId, UniProt, Target, TargetFullName, EntrezGeneID, EntrezGeneSymbol) %>%
  dplyr::rename(UNIPROT = UniProt) %>%
  dplyr::mutate(EntrezGeneID = as.integer(EntrezGeneID))

id_uniprot_id <- unique(data$UNIPROT)
id_entrez_id <- unique(data$EntrezGeneID)
id_hgnc_id <- unique(data$EntrezGeneSymbol)

# mart ====
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
VAR_build <- listDatasets(mart)[grep("hsapiens", listDatasets(mart)$dataset), "description"] # check build
VAR_build <- sub(".*\\((.*)\\).*", "\\1", VAR_build)
mart_attributes <- listAttributes(mart) %>%
  filter(page == "feature_page")

# map each ID ====
## map UNIPROT ====
map_uniprot <- biomaRt_getBM_batch(
  mart = mart,
  attributes = c("uniprot_gn_id", "uniprot_gn_symbol",
                 "entrezgene_id", "entrezgene_accession",
                 "hgnc_id", "hgnc_symbol",
                 "external_gene_name",
                 "gene_biotype",
                 "chromosome_name"),
  filters = "uniprot_gn_id",
  values = id_uniprot_id,
  chunk_size = 2
)

## map ENTREZ ====
map_entrez <- biomaRt_getBM_batch(
  mart = mart,
  attributes = c("uniprot_gn_id", "uniprot_gn_symbol",
                 "entrezgene_id", "entrezgene_accession",
                 "hgnc_id", "hgnc_symbol",
                 "external_gene_name",
                 "gene_biotype",
                 "chromosome_name"),
  filters = "entrezgene_id",
  values = as.integer(id_entrez_id),
  chunk_size = 2
)

## map HGNC ====
map_hgnc <- biomaRt_getBM_batch(
  mart = mart,
  attributes = c("uniprot_gn_id", "uniprot_gn_symbol",
                 "entrezgene_id", "entrezgene_accession",
                 "hgnc_id", "hgnc_symbol",
                 "external_gene_name",
                 "gene_biotype",
                 "chromosome_name"),
  filters = "hgnc_symbol",
  values = id_hgnc_id,
  chunk_size = 2
)

# ensemble primary id ====
ensembl_uniprot <- getBM(
  attributes = c("uniprot_gn_id", "ensembl_gene_id", "chromosome_name"),
  filters = "uniprot_gn_id",
  values = id_uniprot_id,
  mart = mart
)

ensembl_entrez <- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id", "chromosome_name"),
  filters = "entrezgene_id",
  values = as.integer(id_entrez_id),
  mart = mart
)

ensembl_hgnc <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name"),
  filters = "hgnc_symbol",
  values = id_hgnc_id,
  mart = mart
)

map_ensembl <- full_join(ensembl_uniprot, ensembl_entrez, by = c("ensembl_gene_id", "chromosome_name")) %>%
  full_join(ensembl_hgnc, by = c("ensembl_gene_id", "chromosome_name")) %>%
  filter(chromosome_name %in% c(as.character(1:22), "X", "Y"))

# combine maps ====
map <- rbind(map_uniprot,
             map_entrez,
             map_hgnc) %>%
  dplyr::distinct()

map <- full_join(map, map_ensembl)

counts <- map %>%
  summarise(across(everything(), ~ sum(is.na(.) | . == ""))) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "Column", values_to = "MissingCount")

# positions ====
## hg19 ====
id_hg19 <- subset(genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                      single.strand.genes.only = FALSE),
                gene_id %in% unique(map$entrezgene_id)) %>%
  as.data.frame() %>%
  dplyr::rename(entrezgene_id = group_name,
                CHR = seqnames,
                START_hg19 = start,
                END_hg19 = end,
                strand_hg19 = strand) %>%
  dplyr::group_by(entrezgene_id) %>%
  dplyr::filter(any(CHR %in% paste0("chr", c(1:22, "X", "Y")))) %>% # Keep only groups with at least one valid chromosome entry
  dplyr::filter(CHR %in% paste0("chr", c(1:22, "X", "Y")) | n() == 1) %>% # Keep valid chromosomes or keep if only one entry exists
  dplyr::ungroup() %>%
  dplyr::mutate(CHR = sub("^chr", "", CHR)) %>%
  dplyr::mutate(entrezgene_id = as.integer(entrezgene_id)) %>%
  droplevels() %>%
  dplyr::select(-group, -width)

## hg38 ====
id_hg38 <- subset(genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                        single.strand.genes.only = FALSE),
                  gene_id %in% unique(map$entrezgene_id)) %>%
  as.data.frame() %>%
  dplyr::rename(entrezgene_id = group_name,
                CHR = seqnames,
                START_hg38 = start,
                END_hg38 = end,
                strand_hg38 = strand) %>%
  dplyr::group_by(entrezgene_id) %>%
  dplyr::filter(any(CHR %in% paste0("chr", c(1:22, "X", "Y")))) %>% # Keep only groups with at least one valid chromosome entry
  dplyr::filter(CHR %in% paste0("chr", c(1:22, "X", "Y")) | n() == 1) %>% # Keep valid chromosomes or keep if only one entry exists
  dplyr::ungroup() %>%
  dplyr::mutate(CHR = sub("^chr", "", CHR)) %>%
  dplyr::mutate(entrezgene_id = as.integer(entrezgene_id)) %>%
  droplevels() %>%
  dplyr::select(-group, -width)

## combine builds ====
id_gene <- dplyr::full_join(id_hg19, id_hg38,
                     by = c("entrezgene_id" = "entrezgene_id",
                            "CHR" = "CHR"))

# join map and build info ====
map <- map %>%
  dplyr::full_join(id_gene,
                   by = c("entrezgene_id" = "entrezgene_id",
                            "CHR" = "CHR"))

# join map and data ====
data_map_uniprot <- data %>%
  dplyr::left_join(map,
                   by = c("UNIPROT" = "uniprot_gn_id"))

data_map_entrez <- data %>%
  dplyr::left_join(map,
                   by = c("EntrezGeneID" = "entrezgene_id"))

data_map_hgnc <- data %>%
  dplyr::left_join(map,
                   by = c("Target" = "hgnc_symbol"))

data_map <- bind_rows(
  data_map_uniprot,
  data_map_entrez,
  data_map_hgnc
) %>%
  dplyr::group_by(SeqId, SomaId, UNIPROT) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ dplyr::coalesce(.x[1], .x[2]))) %>%
  dplyr::ungroup() %>%
  distinct() %>%
  as.data.frame()

columns <- c(
  "SeqId", "SomaId", "UNIPROT", "Target", "TargetFullName", "EntrezGeneID", "EntrezGeneSymbol",
  "uniprot_gn_id", "uniprot_gn_symbol",
  "entrezgene_id","entrezgene_accession",
  "hgnc_id", "hgnc_symbol", "ensembl_gene_id", "external_gene_name","gene_biotype",
  "CHR", "START_hg19", "END_hg19", "strand_hg19","START_hg38", "END_hg38", "strand_hg38"
)
data_map <- data_map[, columns]

counts <- data_map %>%
  summarise(across(everything(), ~ sum(is.na(.) | . == ""))) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "Column", values_to = "MissingCount")

## write
write.table(x = data_map,
            file = paste0("inst/data/mapping_", VAR_build, "_somalogic.txt"),
            sep = "\t", col.names = T, row.names = F, quote = FALSE)
# save(data_map, file = paste0("data/mapping_", VAR_build, "_somalogic.RData"))
mapping_GRCh38_p14_somalogic <- data_map
usethis::use_data(mapping_GRCh38_p14_somalogic, overwrite = TRUE)
