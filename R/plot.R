#' Plot an upsetplot
#'
#' Upset plot to visualize which variants co-occur.
#'
#' @param h5f h5f
#' @param variants A named vector.
#'  You should get the value from the \code{id} column after calling \code{get_variants}.
#'  If you don't provide the names of the variants, R will do a bad job for you.
#' @param unknown A String. Specifies if NGT = 3 should be considered as reference or discarded for plotting.
#' @param return_data Boolean. Defaults to FALSE. If TRUE, returns the NGT data.table with cells as rows and variants as columns.
#' @param index_cell A numeric vector. Use this if you want to restrict the search to a subset of cells.
#' @param nsets An integer. How many sets (ie, variants) of the supplied variants to plot.
#' @examples
#' path <- "path/to/.h5"
#' h5f <- read_h5f(path)
#'
#' interesting_variants <- c(IDH1_R132H = "chr2:209113112:C/T", NPM1c = "chr5:170837543:C/CTCTG")
#' plot_upset_with_variants(h5f, variants = interesting_variants)
plot_upset_with_variants <- function(h5f,
                                     variants,
                                     unknown = c("unknown_to_ref", "unknown_discard"),
                                     return_data = FALSE,
                                     index_cell = NULL,
                                     nsets = 5){


  dt.var <- get_variants(h5f)

  index_var <- get_variants_index(h5f = h5f, variants = variants, sort = TRUE)

  dt.ngt <- read_assays_variants(h5f = h5f, included_assays = "NGT", index_cells = index_cell, index_variants = index_var)
  dt.ngt <- dt.ngt[[1]]
  dt.ngt <- data.table(t(dt.ngt))


  colnames(dt.ngt) <- names(variants)[order(match(variants, dt.var[index_var, ]$id))]

  if(unknown == "unknown_to_ref"){
    dt.ngt <- dt.ngt %>%
      apply(2, function(x){gsub(pattern = "2", replacement = "1", x = x)}) %>%
      apply(2, function(x){gsub(pattern = "3", replacement = "0", x = x)}) %>%
      data.table() %>%
      mutate_all(as.numeric)
    if(return_data == FALSE){
      return(UpSetR::upset(dt.ngt, nsets = nsets))
    } else {
      return(dt.ngt)
    }
  }
  if(unknown == "unknown_discard"){
    dt.ngt <- dt.ngt %>%
      apply(2, function(x){gsub(pattern = "2", replacement = "1", x = x)}) %>%
      data.table() %>%
      mutate_all(as.numeric)

    indx_unknown <- apply(dt.ngt, 1, function(x){any(x == 3)})
    dt.ngt <- dt.ngt[!indx_unknown, ]
    if(return_data == FALSE){
      return(UpSetR::upset(dt.ngt, nsets = nsets))
    } else {
      return(dt.ngt)
    }
  }
}
