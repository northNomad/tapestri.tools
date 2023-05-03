#' Reading a tapestri h5 file
#'
#' Wrapper function for H5Fopen from rhdf5. Reads in a \code{.h5} file.
#'
#' @param path A character, path to .h5
#' @return h5f
#' @examples
#' path <- "path/to/.h5"
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
#' @param index_variants A named numeric vector. You should use \code{get_variants_index} to retrieve the indexes.
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




#' Reading protein read counts
#'
#' read_protein_counts retrieves the protein read counts in an h5f.
#'
#' @param h5f h5f
#' @param transpose Boolean. If \code{TRUE}, columns represent proteins and rows represent cells.
#' @param index_cells A numeric vector. Specifies which cells to retrieve.
#' @param normalization Either \code{clr} or \code{raw}. Defaults to \code{clr}, as implemented by the \code{compositions} package.
#' @return matrix
#' @examples
#' path <- "path/to/.h5"
#' h5f <- read_h5(path)
#'
#' read_protein_counts(h5f)
read_protein_counts <- function(h5f,
                                transpose = c(TRUE, FALSE),
                                index_cells = NULL,
                                normalization = c("clr", "raw")){
  x <- rhdf5::h5read(h5f, "/assays/protein_read_counts/layers/read_counts", index = list(NULL, index_cells))
  id <- rhdf5::h5read(h5f, "/assays/protein_read_counts/ca/id")

  if(normalization == "clr"){
    x <- base::apply(x, 2, compositions::clr, simplify = FALSE) %>% do.call(cbind, .)
  }

  if(transpose){
    x <- t(x)
    colnames(x) <- id
    return(x)
  } else{
    rownames(x) <- id
    return(x)
  }
}




#' Read tapestri amplicons
#'
#' Read bed file describing coverage of tapestri amplicons
#'
#' @param bed path to a 4 column .bed file.
#'   The four columns are: chrom, chromStart, chromEnd, name.
#' @param format Either \code{GRanges} or \code{data.table}. Defaults to \code{GRanges}.
#' @examples
#' \dontrun{
#'   gr_amplicons <- read_amplicons_bed("path/to/bed")
#' }
read_amplicons_bed <- function(bed,
                               format = c("GRanges", "data.table")){

  name_col <- c("chrom", "chromStart", "chromEnd", "name")
  x <- read.delim(bed, col.names = name_col, header = FALSE) %>% data.table()

  if(format == "GRanges"){
    x <- GRanges(seqnames = x$chrom,
                 ranges = IRanges(start = x$chromStart, end = x$chromEnd),
                 mcols = data.frame(name = x$name)
                 )
    return(x)
  } else{
    return(x)
  }
}




#' Read ENSEMBL VEP output
#'
#' Reads standard VEP output.
#'
#' @param file path to vep .txt output file.
#' @examples
#' \dontrun{
#'   write_vep_input(dt.var, "vep_input.txt")
#'   #run vep locally
#'   vep <- read_vep_output("path_vep_output.txt")
#' }
read_vep_output <- function(file){

  linesToSkip <- grep("##", readLines(file))
  linesToSkip <- linesToSkip[length(linesToSkip)]

  dt.out <- read.delim(file, sep = "\t", skip = linesToSkip)
  dt.out <- data.table(dt.out)

  return(dt.out)
}





#' Read SCITE output
#'
#' Designate the clone of each cell based on SCITE. If ambiguous, a random clone from all possible clones is selected.
#'
#' @param file path to SCITE output file in .gv format.
#' @param cell_index_h5f A numeric vector. The index of cells used for scite in the .h5f
#' @param seed integer.
#' @return data.table
#' @examples
#' \dontrun{
#'  read_scite_out("scite_output.gv")
#' }
read_scite_output <- function(file, cell_index_h5f, seed = 123){

  #assumes -a is turned on in scite
  linesToSkip <- grep("node", readLines(file))
  linesToSkip <- linesToSkip[length(linesToSkip)]

  #read
  dt.out <- read.delim(file, sep = "\t", skip = linesToSkip, header = FALSE)
  dt.out <- dt.out[-nrow(dt.out), , drop = FALSE]
  dt.out <- data.table(dt.out)

  #regex manipulation
  dt.out$Clone <- dt.out$V1 %>% str_split(" -> ") %>% lapply(function(x){x[1]}) %>% unlist()
  #more regex manipulation
  dt.out$scite_index <- dt.out$V1 %>%
    str_split(" -> ") %>%
    lapply(function(x){
      x <- x[2]
      x <- gsub("s", "", x)
      x <- gsub(";", "", x)
      x <- as.numeric(x) + 1
      x <- as.character(x)
    }) %>%
    unlist()

  #map scite index back to original cell index in h5f
  dt.out$cell_index <- cell_index_h5f[as.numeric(dt.out$scite_index)]

  #remove scite index. No longer needed.
  dt.out <- dplyr::select(dt.out, Clone, cell_index)

  #randomly select clone for ambiguious cells
  dt.out %>%
    group_by(cell_index) %>%
    group_modify(.f = function(x, y){
      set.seed(seed = seed)
      n <- nrow(x)
      n <- sample(n, 1)
      x[n, ]
    }) -> dt.out
  dt.out <- data.table(dt.out)
  return(dt.out)
}
