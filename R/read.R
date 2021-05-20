#' Reading a tapestri h5 file
#'
#' Wrapper function for H5Fopen from rhdf5. Reads in a \code{.h5} file.
#'
#' @param path A character, path to .h5
#' @return h5f
#' @examples
#' path <- "path/to/.h5"
#' h5f <- read_h5f(path)
#' h5f <- read_h5(path)
read_h5 <- function(path){
  h5f <- rhdf5::H5Fopen(path, flags = "H5F_ACC_RDONLY")
  return(h5f)
}





#' Reading the assays of variants
#'
#' read_assays_variants retrieves the assays of variants in an h5f.
#'
#' @param h5f h5f
#' @param included_assays A character vector. Selects which assays to retrieve.
#'   Possible values are \describe{
#'     \item{AF}{variant allele frequencies (VAF) per cell and per variant}
#'     \item{DP}{read counts per cell and per variant}
#'     \item{FILTER_MASK}{binary value}
#'     \item{GQ}{genotype quality values per cell and per variant}
#'     \item{NGT}{genotypes per cell and per variant \itemize{
#'       \item 0 = Reference (Wildtype)
#'       \item 1 = Heterozygous Mutant
#'       \item 2 = Homozygous Mutant
#'       \item 3 = Missing Data
#'       }
#'     }
#'     \item{RGQ}{???}
#'   }
#' @param index_cells A named numeric vector. You should use \code{get_variants_index} to retrieve the indexes.
#' @param index_cells A numeric vector. Specifies which cells (columns) to retrieve.
#' @param format A character. Specifies the return format. Either a \code{list} or a \code{RangedSummarizedExperiment}.
#' @return Default is a named \code{list}. Can return a \code{RangedSummarizedExperiment} if specified.
#' @examples
#' #' #Retrieving the NGT matrix for IDH1_R132H and NPM1c
#' \dontrun{
#' interesting_variants <- c(IDH1_R132H = "chr2:209113112:C/T", NPM1c = "chr5:170837543:C/CTCTG")
#' index <- get_variants_index(h5f, interesting_variants)
#' read_dna_variants(h5f, index_variants = index, included_assays = "NGT")
#' }
read_assays_variants <- function(h5f,
                                 included_assays = c("AF", "DP", "FILTER_MASK", "GQ", "NGT", "RGQ"),
                                 index_variants = NULL,
                                 index_cells = NULL,
                                 format = c("list", "RangedSummarizedExperiment")){

  assays_list <- list()
    for(i in 1:length(included_assays)){
      fp <- paste0("/assays/dna_variants/layers/", included_assays[i])
      dt <- rhdf5::h5read(h5f, name = fp, index = list(index_variants, index_cells)) %>% as.matrix()
      rownames(dt) <-  names(index_variants)
      assays_list[[i]] <- dt
    }
  names(assays_list) <- included_assays

  if(format == "list"){
   return(assays_list)
  }

  if(format == "RangedSummarizedExperiment"){
    row_data <- get_variants(h5f, "data.table")[index_variants, ]
    rowRanges <- GRanges(seqnames = row_data$CHROM,
                         ranges = IRanges(start = as.numeric(row_data$POS),
                                          width = nchar(row_data$REF)),
                         mcols = row_data)
    colnames(mcols(rowRanges)) <- colnames(mcols(rowRanges)) %>% gsub(pattern = "mcols\\.", replacement = "")

    SummarizedExperiment(assays = assays_list,
                         rowRanges = rowRanges
                         ) -> se
    return(se)
  }
}

