#' Getting variants called by haplotype caller in tapestri pipeline
#'
#' get_variant retrieves all the dna variants in an .h5 file as a data.table or GRanges
#'
#' @param h5f h5f
#' @param format \code{data.table} by default.
#' @examples
#' path <- "path/to/.h5"
#' h5f <- read_h5f(path)
#' feature_data <- get_variants(h5f)
get_variants <- function(h5f,
                         format = c("data.table", "GRanges")){

  dt.var <- rhdf5::h5read(h5f, "/assays/dna_variants/ca") %>%
    do.call(cbind, .) %>%
    data.table::data.table()

  dt.var$SYMBOL <- dt.var$amplicon %>%
    stringr::str_split("_") %>%
    lapply(function(x){x[3]}) %>%
    unlist()

  if(isFALSE(grepl(pattern = "chr", ignore.case = TRUE, x = dt.var$CHROM[1]))){
    dt.var$CHROM <- paste0("chr", dt.var$CHROM)
  }

  if(format == "GRanges"){
    dt.var <- GRanges(seqnames = dt.var$CHROM,
                      ranges = IRanges(start = as.numeric(dt.var$POS), width = 1),
                      mcols = dplyr::select(dt.var, -c("CHROM", "POS")))
    colnames(mcols(dt.var)) <- colnames(mcols(dt.var)) %>% gsub(pattern = "mcols\\.", replacement = "")
    return(dt.var)
  } else
    return(dt.var)
}




#' Getting FLT3 internal tandem duplication (ITD) variants
#'
#' get_flt3itd retrieves FLT3 ITD variants in an .h5 file
#'
#' @param h5f h5f
#' @param format \code{data.table} by default.
#' @param gr_coordinate_within A \code{GRanges} object.
#'  Specifies the genomic coordinates to search for ITD mutations in the FLT3 gene.
#'  Currently uses hg19 coordinates as reference.
#' @param insertion_size A numeric vector with exactly two values.
#'  Specifies the minimum and maximum insertion size (in basepairs) introduced by an ITD mutation.
#' @return A \code{data.table} or a \code{GRanges} object with FLT3-ITD variants.
get_flt3itd <- function(h5f,
                        format = c("data.table", "GRanges"),
                        gr_coordinate_within = GRanges(seqnames = "chr13",
                                                       ranges = IRanges(start = 28607990, end = 28608531)
                                                       ), #exon14-15
                        insertion_size = c(6, 180)){

  var_flt3 <- get_variants(h5f, format = "GRanges")
  var_flt3 <- var_flt3[var_flt3$SYMBOL == "FLT3"]

  var_flt3 <- var_flt3[var_flt3 %within% gr_coordinate_within] #Select the variants within `gr_coordinate_within`

  length_alt <- mcols(var_flt3)$id %>% str_split(pattern = "/") %>% lapply(function(x) x[2]) %>% unlist() %>% nchar()
  length_ref <- mcols(var_flt3)$REF %>% gsub(pattern = "\\*", replacement = "") %>% gsub(pattern = "\\.", replacement = "") %>% nchar()
  length_insert <- length_alt - length_ref
  indx <- length_insert %in% seq(from = insertion_size[1], to = insertion_size[2], by = 3)

  var_flt3 <- var_flt3[indx]

  if(format == "GRanges"){return(var_flt3)}

  if(format == "data.table"){
    var_flt3 <- var_flt3 %>% as.data.frame() %>% data.table()
    var_flt3 <- var_flt3 %>%
      dplyr::rename(POS = "start", CHROM = "seqnames") %>%
      dplyr::select(ALT, CHROM, POS, QUAL, REF, ado_gt_cells, ado_rate, amplicon, filtered, id, SYMBOL)
    return(var_flt3)
  }
}
