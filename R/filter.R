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
}





#' Filtering by genotype quality
#'
#' Removes genotype calls in cells with quality < cutoff
#'
#' @param sce A \code{SingleCellExperiment} object, with at least a \code{NGT} and a \code{GQ} slot.
#' @param cutoff Genotype quality cutoff (Phred scale). Default is 30.
#' @return A \code{SingleCellExperiment} object with modified \code{NGT} slot.
filter_GQ <- function(sce, cutoff = 30){

  low_GQ <- assays(sce)[["GQ"]] < cutoff
  ngt <- assays(sce)[["NGT"]]
  ngt[low_GQ] <- 3
  assays(sce)[["NGT"]] <- ngt

  return(sce)
}





#' Filtering by sequencing depth
#'
#' Removes genotype calls in cells with depth < cutoff
#'
#' @param sce A \code{SingleCellExperiment} object, with at least a \code{NGT} and a \code{DP} slot.
#' @param cutoff Sequencing depth cutoff. Default is 10.
#' @return A \code{SingleCellExperiment} object with modified \code{NGT} slot.
filter_DP <- function(sce, cutoff = 10){

  low_DP <- assays(sce)[["GQ"]] < cutoff
  ngt <- assays(sce)[["NGT"]]
  ngt[low_DP] <- 3
  assays(sce)[["NGT"]] <- ngt

  return(sce)
}





#' Filtering by variant allele frequency
#'
#' Removes genotype calls in cells with variant allele frequency < cutoff
#'
#' @param sce A \code{SingleCellExperiment} object, with at least a \code{NGT} and a \code{AF} slot.
#' @param cutoff Variant allele frequency cutoff. Default is 20. The unit is in percentage.
#' @return A \code{SingleCellExperiment} object with modified \code{NGT} slot.
filter_AF <- function(sce, cutoff = 20){

  low_AF <- assays(sce)[["AF"]] < cutoff
  ngt <- assays(sce)[["NGT"]]
  ngt[low_AF] <- 3
  assays(sce)[["NGT"]] <- ngt

  return(sce)
}





#' Filtering variants with poor genotyping performance
#'
#' Removes variants genotyped in < cutoff% of cells
#'
#' @param sce A \code{SingleCellExperiment} object, with at least a \code{NGT}.
#' @param cutoff The unit is in percentage. The default is 50.
#' @return A modified \code{SingleCellExperiment} (containing fewer variants).
filter_variants_percent_genotyped <- function(sce, cutoff = 50){
  n_Cells <- ncol(sce) #Number of cells

  ngt <- assays(sce)[["NGT"]]
  ngt <- data.table(ngt)

  ngt <- apply(ngt, 1, function(x) sum(x == 3)) #Count number of ungenotyped cells
  ngt <- ngt * 100 / n_Cells #Calculate percentage of ungenotyped cells
  ngt <- ngt > cutoff #Index of low quality variants

  sce <- sce[-ngt]

  return(sce)
}





#' Filtering cells with poor genotyping performance
#'
#' Removes cells with <cutoff% of variants genotyped
#'
#' @param sce A \code{SingleCellExperiment} object, with at least a \code{NGT}.
#' @param cutoff The unit is in percentage. The default is 50.
#' @return A modified \code{SingleCellExperiment} (containing fewer cells).
filter_cells_percent_genotyped <- function(sce, cutoff = 50){
  n_Variants <- nrow(sce) #Number of cells

  ngt <- assays(sce)[["NGT"]]
  ngt <- data.table(ngt)

  ngt <- apply(ngt, 2, function(x) sum(x == 3)) #Count number of ungenotyped variants
  ngt <- ngt * 100 / n_Variants #Calculate percentage of ungenotyped variants
  ngt <- ngt > cutoff #Index of low quality variants

  sce <- sce[, -ngt]

  return(sce)
}





#' Filtering variants present in very few cells
#'
#' Removes variants mutated in < cutoff% of cells
#'
#' @param sce A \code{SingleCellExperiment} object, with at least a \code{NGT}.
#' @param cutoff The unit is in percentage. The default is 1.
#' @return A modified \code{SingleCellExperiment} (containing fewer variants).
filter_variants_percent_mutated <- function(sce, cutoff = 1){
  n_Cells <- ncol(sce) #Number of cells

  ngt <- assays(sce)[["NGT"]]
  ngt <- data.table(ngt)

  ngt <- apply(ngt, 1, function(x) sum(x == 1 | x == 2)) #Count number of mutated cells
  ngt <- ngt * 100 / n_Cells #Calculate percentage of mutated cells
  ngt <- ngt < cutoff #Index of low frequency variants

  sce <- sce[-ngt]

  return(sce)
}
