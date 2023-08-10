#' Write input files for SCITE
#'
#' A helper function to write genotype matrix and .geneNames files for SCITE.
#'
#' @param h5f h5f
#' @param variants A named vector.
#'  You should get the value from the \code{id} column after calling \code{get_variants}.
#'  The names should be meaningful, they are the labels in your inferred clonal architecture .gv file.
#'  Avoid symbols like \code{., -, :} in the names, and use an underscore \code{_} instead.
#' @param cell_index A numeric vector. Specifies which cells (columns) to retrieve.
#' @param fp_ngt_csv Output file name for NGT matrix .csv
#' @param fp_geneNames Output file name for .geneNames
#' @return NULL. Writes two files.
#' @references https://github.com/cbg-ethz/SCITE
#' @references Kuipers J et al. Single-cell sequencing data reveals widespread recurrence and loss of mutational hits in the life histories of tumorus. Genome Research 2017; 27:1885-1894.
write_scite_inputs <- function(h5f,
                               variants,
                               cell_index = NULL,
                               fp_ngt_csv,
                               fp_geneNames){

  indx_var <- get_variants_index(h5f, variants, sort = TRUE)

  dt_ngt <- read_assays_variants(h5f, included_assays = "NGT", index_variants = indx_var, index_cells = cell_index, format = "list")
  dt_ngt <- as.data.frame(dt_ngt)

  write_delim(dt_ngt, file = fp_ngt_csv, delim = " ", col_names = FALSE)

  fileConn <- file(fp_geneNames)
  writeLines(names(indx_var), fileConn)
  close(fileConn)
}





#' Parse SCITE mutation tree in parent vector format into mermaid code
#'
#' A helper function to write mermaid code. Used for visualization with \code{DiagrammeR::grViz()}.
#'
#' @param ParentVector A (numeric) vector that uses zero-based indexing.
#' Ideally you'd get this vector from the posterior *.samples* file from SCITE output.
#' @param Variants A character vector containing the names of the variants.
#' It should be in the same order as the rows in the matrix you used as input for SCITE.
#' @return A string. Can be directly used for visualization with \code{DiagrammeR::grViz()}
#' #' @references https://github.com/cbg-ethz/SCITE
#' @references Kuipers J et al. Single-cell sequencing data reveals widespread recurrence and loss of mutational hits in the life histories of tumorus. Genome Research 2017; 27:1885-1894.
ParentVector_to_DiagrammeR <- function(ParentVector, Variants){
  ParentVector <- ParentVector[1:length(Variants)]
  ParentVector <- as.numeric(ParentVector)
  ParentVector <- ParentVector + 1 #SCITE uses zero indexing.
  out <- "digraph G {\nnode [color=deeppink4, style=filled, fontcolor=white];\n"

  for(i in 1:length(Variants)){
    Child <- Variants[i]
    if(ParentVector[i] > length(Variants)){
      Parent <- "Root"
    }else{
      Parent <- Variants[ParentVector[i]]
    }

    Code <- paste0(Parent, " -> ", Child, ";", "\n")
    out <- paste0(out, Code)
  }
  out <- paste0(out, "}")
  return(out)
}





#' Write input files for ENSEMBL VEP
#'
#' A helper function to write standard input files for ENSEMBL VEP.
#'
#' @param dt_variants A (filtered) \code{data.frame} or \code{data.table} of variants derived from \code{get_variants}
#' @param file file path to write .txt file.
#' @return NULL. Writes the standard VEP .txt input file.
#' #' @examples
#'
#' dt.var <- get_variants(h5f)
#' write_vep_input(dt.var, "vep_input.txt")
#' @references https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html
write_vep_input <- function(dt_variants,
                            file){

  dt_variants %>%
    mutate(chromosome = gsub("chr", "", CHROM),
           start = as.numeric(POS),
           end = start - (nchar(ALT) - nchar(REF)),
           allele = paste(REF, ALT, sep = "/"),
           strand = "+",
           identifier = paste(SYMBOL, id, sep = "_")
    ) %>%
    dplyr::select(chromosome, start, end, allele, strand, identifier) -> out

  write_delim(out, file, delim = " ", col_names = FALSE)
}


##vep --offline --dir_cache $cache_dir -i ~/example_input.txt -o ~/example_vep_output.txt -merged --assembly GRCh37 --use_given_ref
