
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
