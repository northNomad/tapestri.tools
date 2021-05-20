#' Retrieves all the dna variants in an .h5 file as a
#'
#' @param h5f
#' @param format Returns a data.table by default
#' @examples
#' feature_data <- get_variants(h5f)

##Get all the variants
get_variants <- function(h5f,
                         format = c("data.table", "GRanges")){

  dt.var <- rhdf5::h5read(h5f, "/assays/dna_variants/ca") %>%
    do.call(cbind, .) %>%
    data.table::data.table()

  dt.var$SYMBOL <- dt.var$amplicon %>%
    stringr::str_split("_") %>%
    lapply(function(x){x[3]}) %>%
    unlist()

  if(isFALSE(grepl(pattern = "chr", ignore.case = TRUE, x = dt.var$CHROM))){
    dt.var$CHROM <- paste0("chr", dt.var$CHROM)
  }

  if(format == "data.table"){
    return(dt.var)
  }

  if(format == "GRanges"){
    dt.var <- GRanges(seqnames = dt.var$CHROM,
                      ranges = IRanges(start = as.numeric(dt.var$POS), width = 1),
                      mcols = dplyr::select(dt.var, -c("CHROM", "POS")))
    return(dt.var)
  }
}
