#' Count the number of mutant cells given variants
#'
#' count_cells tally the number of called genotypes (NGT slot) of given variants
#'
#' @param x input, either an .h5f or a SingleCellExperiment object.
#' @param slot Slot storing the genotype calls. Used if x is a SingleCellExperiment object.
#'             Default is "NGT".
#' @param format \code{matrix} by default.
count_cells <- function(x,
                        variants, 
                        slot = "NGT",
                        index_cells = NULL, 
                        percent_mutated = TRUE, 
                        percent_genotyped = TRUE){
  #Check if x is an hf5 file or an sce object
  if(class(x) != "H5IdComponent" && class(x) != "SingleCellExperiment"){
    stop("The input file needs to be an h5f or a SingleCellExperiment object.")
  }
  
  #Check if the variants are named
  if(is.null(names(variants))){
    names(variants) <- paste0("Variant", 1:length(variants))
    message("You didn't provide the names of the variants. The package is naming them for you.")
  }
  
  #Check if variants are present
  if(class(x) == "H5IdComponent"){
    detected <- variants %in% get_variants(x)$id
  }
  if(class(x) == "SingleCellExperiment"){
    detected <- variants %in% rowData(x)[, "id"]
  }
  
  #Report any undetected variants
  if(all(detected) == FALSE){
    variants_undetected <- paste(variants[!detected], sep = ";")
    message(paste0("These variants are not detected: ", variants_undetected))
  }
  
  #Keep detected variants
  variants <- variants[detected]
  
  #get called genotype of variants
  if(class(x) == "H5IdComponent"){
    index <- get_variants_index(x, variants)
    dt.ngt <- tapestri.tools::read_assays_variants(x, 
                                                   included_assays="NGT", 
                                                   index_cells=index_cells, 
                                                   index_variants=index, 
                                                   format="list")[[1]]
  }
  if(class(x) == "SingleCellExperiment"){
    index <- match(names(variants), rownames(x))
    index <- index[!is.na(index)]
    dt.ngt <- assays(x)[[slot]][index, ]
  }
  
  #make dt.ngt a matrix in case if user only supplies one variant.
  #Need to do this because 'apply' expects input with dimensions.
  if(is.null(dim(dt.ngt))){
    dt.ngt <- matrix(dt.ngt, nrow=1)
    rownames(dt.ngt) <- names(variants)
  }
  
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
}




