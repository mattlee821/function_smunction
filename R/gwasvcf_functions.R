#' Find LD proxies for a set of SNPs
#'
#' @param rsid list of rs IDs
#' @param bfile ld reference panel
#' @param tag_kb =5000 Proxy parameter
#' @param tag_nsnp =5000 Proxy parameter
#' @param tag_r2 =0.6 Proxy parameter
#' @param searchspace Optional list of rs IDs to use as potential proxies
#' @param threads Number of threads to use (=1)
#' @param out temporary output file
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#' @return data frame
get_ld_proxies <- function(rsid, bfile, searchspace=NULL, tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1, out=tempfile())
{
  stopifnot(check_plink())
  searchspacename <- paste0(out, ".searchspace")
  targetsname <- paste0(out, ".targets")
  outname <- paste0(out, ".targets.ld.gz")
  utils::write.table(rsid, file=targetsname, row.names = FALSE, col.names = FALSE, quote = FALSE)
  if(!is.null(searchspace))
  {
    stopifnot(is.character(searchspace))

    utils::write.table(unique(c(rsid, searchspace)), file=searchspacename, row.names = FALSE, col.names = FALSE, quote = FALSE)
    extract_param <- paste0(" --extract ", searchspacename)
  } else {
    extract_param <- " "
  }
  cmd <- paste0(options()[["tools_plink"]],
                " --bfile ", bfile,
                extract_param,
                " --keep-allele-order ",
                " --r in-phase with-freqs gz",
                " --ld-snp-list ", targetsname,
                " --ld-window-kb ", tag_kb,
                " --ld-window-r2 ", tag_r2,
                " --ld-window ", tag_nsnp,
                " --out ", targetsname,
                " --threads ", threads,
                " 2>&1 > /dev/null"
  )
  message("Finding proxies...")
  system(cmd)

  if (Sys.info()["sysname"] == "Windows") {
    stop("Currently, this function only works on macOS and Linux")
  }
  if (!file.exists(outname)) {
    ld <- data.frame(CHR_A = integer(), BP_A = integer(), SNP_A = character(), MAF_A = double(), CHR_B = integer(), BP_B = integer(),
                     SNP_B = character(), PHASE = character(), MAF_B = double(), R = double())
    message("Index SNP not found in the reference panel")
    return(ld)
  }
  ld <- data.table::fread(cmd = paste0("gunzip -c ", shQuote(outname)), header = TRUE) %>%
    dplyr::as_tibble(.name_repair="minimal") %>%
    dplyr::filter(.data[["R"]]^2 > tag_r2) %>%
    dplyr::filter(.data[["SNP_A"]] != .data[["SNP_B"]]) %>%
    dplyr::mutate(PHASE=gsub("/", "", .data[["PHASE"]])) %>%
    dplyr::filter(nchar(.data[["PHASE"]]) == 4)
  unlink(searchspacename)
  unlink(targetsname)
  unlink(paste0(targetsname, c(".log", ".nosex")))
  unlink(outname)
  if(nrow(ld) == 0)
  {
    message("No proxies found")
    return(ld)
  }
  temp <- do.call(rbind, strsplit(ld[["PHASE"]], "")) %>% dplyr::as_tibble(.name_repair="minimal")
  names(temp) <- c("A1", "B1", "A2", "B2")
  ld <- cbind(ld, temp) %>% dplyr::as_tibble(.name_repair="minimal")
  # ld <- dplyr::arrange(ld, desc(abs(R))) %>%
  # 	dplyr::filter(!duplicated(SNP_A))
  ld <- dplyr::arrange(ld, dplyr::desc(abs(.data[["R"]])))
  message("Found ", nrow(ld), " proxies")
  return(ld)
}
