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
count_cells <- function(h5f, variants, index_cells = NULL, percent_mutated = TRUE, percent_genotyped = TRUE){

  #Check if the variants are named
  if(is.null(names(variants))){
    names(variants) <- paste0("Variant", 1:length(variants))
  }

  #Check if variants are present
  detected <- variants %in% get_variants(h5f)$id

  #Report any undetected variants
  if(all(detected) == FALSE){
    variants_undetected <- paste(variants[!detected], sep = ";")
    message(paste0("These variants are not detected: ", variants_undetected))
    }

  #Keep detected variants
  variants <- variants[detected]

  #get called genotype of variants
  index <- get_variants_index(h5f, variants)
  dt.ngt <- tapestri.tools::read_assays_variants(h5f, included_assays = "NGT", index_cells = index_cells, index_variants = index, format = "list")[[1]]

  #Count
  NGT0 <- apply(dt.ngt, 1, function(x) sum(x==0))
  NGT1 <- apply(dt.ngt, 1, function(x) sum(x==1))
  NGT2 <- apply(dt.ngt, 1, function(x) sum(x==2))
  NGT3 <- apply(dt.ngt, 1, function(x) sum(x==3))

  m <- data.frame(Variant = names(NGT0), NGT0 = NGT0, NGT1 = NGT1, NGT2 = NGT2, NGT3 = NGT3)

  if(percent_mutated == TRUE){
    m$percent_mutated <- with(m, (NGT1+NGT2)*100 / (NGT0+NGT1+NGT2))
  }

  if(percent_genotyped == TRUE){
    m$percent_genotyped <- with(m, (NGT0+NGT1+NGT2)*100 / (NGT0+NGT1+NGT2+NGT3))
  }

  return(m)
  # #count
  # dt.ngt <- apply(dt.ngt, 1, table)
  #
  # #Check class (will be different if not all variants have values of 0, 1, 2, 3)
  # class <- class(dt.ngt)[1]
  #
  # #if matrix
  # if(class == "matrix"){
  #   m <- dt.ngt %>% as.matrix() %>% t() %>% as.data.frame()
  #   m[, 5] <- rownames(m)
  #   m <- m[, c(5, 1, 2, 3, 4)]
  #
  #   colnames(m) <- c("Variant", "NGT0", "NGT1", "NGT2", "NGT3")
  # }
  #
  # #if list
  # if(class == "list"){
  #   #where to place the counts in output matrix
  #   col_index_list <- lapply(dt.ngt, function(x) as.integer(unlist(dimnames(x))) + 1)
  #
  #   #initialize output matrix
  #   m <- matrix(ncol = 4, nrow = length(variants), dimnames = list(names(variants), paste0("NGT", 0:3)))
  #
  #   #write output
  #   for(i in 1:length(variants)){
  #     m[i, col_index_list[[i]]] <- as.numeric(dt.ngt[[i]])
  #   }
  #   m <- as.data.frame(m)
  #   m[, 5] <- rownames(m)
  #   m <- m[, c(5, 1, 2, 3, 4)]
  #
  #   colnames(m) <- c("Variant", "NGT0", "NGT1", "NGT2", "NGT3")
  # }
  #
  # #replace NA with zeros
  # m[is.na(m)] <- 0
  # m
}
