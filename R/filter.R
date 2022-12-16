# filter_variants_by_gq <- function(h5f,
#                             gqc = 30,
#                             dpc = 10,
#                             afc = 20,
#                             mv = 50,
#                             mc = 50,
#                             mm = 1,
#                             gt_mask = FALSE){
#
#   assays <- read_assays_variants(h5f, included_assays = c("AF", "DP", "GQ", "NGT"), format = "list")
#
# }





#' Count the number of mutant cells given variants
#'
#' count_cells tally the number of called genotypes (NGT slot) of given variants
#'
#' @param h5f h5f
#' @param format \code{matrix} by default.
count_cells <- function(h5f, variants){

  #get called genotype of variants
  index <- get_variants_index(h5f, variants)
  dt.ngt <- tapestri.tools::read_assays_variants(h5.ena, included_assays = "NGT", index_variants = index, format = "list")[[1]]

  #count
  dt.ngt <- apply(dt.ngt, 1, table)

  #where to place the counts in output matrix
  col_index_list <- lapply(dt.ngt, function(x) as.integer(unlist(dimnames(x))) + 1)

  #initialize output matrix
  m <- matrix(ncol = 4, nrow = length(variants), dimnames = list(names(variants), paste0("NGT", 0:3)))

  #write output
  for(i in 1:length(variants)){
    m[i, col_index_list[[i]]] <- as.numeric(dt.ngt[[i]])
  }

  #replace NA with zeros
  m[is.na(m)] <- 0
  m
}
