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





#' Getting the row index of interesting variants
#'
#' Given a vector of named variants, get_variants_index returns their row indexes in /assays/dna_variants/layers of an h5f.
#' It is useful to index which variants (rows) to return when retrieving NGT matrices.
#'
#' @param h5f h5f
#' @param variants A named vector.
#'  You should get the value from the \code{id} column after calling \code{get_variants}.
#'  If you don't provide the names of the variants, R will do a bad job for you.
#' @param sort TRUE or FALSE. Whether to sort indexes by genomic coordinates. Defaults to TRUE.
#' @return A named numeric vector
#' @examples
#' #Retrieving the NGT matrix for IDH1_R132H and NPM1c
#' \dontrun{
#' interesting_variants <- c(IDH1_R132H = "chr2:209113112:C/T", NPM1c = "chr5:170837543:C/CTCTG")
#' index <- get_variants_index(h5f, interesting_variants)
#' read_dna_variants(h5f, index_variants = index, included_assays = "NGT")
#' }
get_variants_index <- function(h5f,
                               variants,
                               sort = TRUE){

  dt.var <- get_variants(h5f, "data.table")
  indx.variants <- match(variants, dt.var$id)
  names(indx.variants) <- names(variants)

  if(sort == TRUE){
    indx.variants <- sort(indx.variants)
    return(indx.variants)
  }

  return(indx.variants)
}




#' Getting the column index of interesting cells
#'
#' get_cells_index finds which cells are called with the queried NGT for all given variants.
#'
#' @param h5f h5f
#' @param variants A named vector.
#'  You should get the value from the \code{id} column after calling \code{get_variants}.
#'  If you don't provide the names of the variants, R will do a bad job for you.
#' @param int_ngt A numeric vector. Specifies genotypes of interest.
#' @param cell_index A numeric vector. Use this if you want to restrict the search to a subset of cells.
#' @return Either a numeric vector or a named list of numeric vectors.
#' @examples
#' #Retrieving the NGT matrix for IDH1_R132H and NPM1c
#' \dontrun{
#'
#' #Want to find which cells are mutated for IDH1_R132H and NPM1c
#' interesting_variants <- c(IDH1_R132H = "chr2:209113112:C/T",
#'                           NPM1c = "chr5:170837543:C/CTCTG")
#' get_cells_index(h5f,
#'                 variants = interesting_variants,
#'                 int_ngt = c(1, 2))
#'
#' #Want to find which cells are unknown for IDH1_R132H
#' get_cells_index(h5f,
#'                 variants = c(IDH1_R132H = "chr2:209113112:C/T"),
#'                 int_ngt = 3)
#'
#' }
##Get index of cells based on variants and wanted NGT
get_cells_index <- function(h5f,
                            variants,
                            int_ngt, cell_index = NULL){

  dt_var <- get_variants(h5f, "data.table")
  indx_var <- get_variants_index(h5f, variants, sort = TRUE)

  ngt <- read_assays_variants(h5f, included_assays = "NGT", index_variants = indx_var, index_cells = cell_index, format = "list")
  ngt <- as.data.frame(ngt)

  apply(ngt, 1, function(x){
    x %in% int_ngt
  }) %>%
    apply(2, function(x){which(x == TRUE)}) -> index_list

  names(index_list) <- names(indx_var)

  if(length(variants) == 1){
    #return a numeric vector if only one variant is supplied
    attr(index_list, "names") <- NULL
    attr(index_list, "dim") <- NULL
  }

  return(index_list)
}


