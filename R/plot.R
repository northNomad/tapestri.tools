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
      return(UpSetR::upset(as.data.frame(dt.ngt), nsets = nsets))
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
      return(UpSetR::upset(as.data.frame(dt.ngt), nsets = nsets))
    } else {
      return(dt.ngt)
    }
  }
}





#' Plot protein abundance across clones
#'
#' Uses ggplot framework to plot protein abundance across groups of cells (aka, clones).
#'
#' @param h5f h5f
#' @param index_cells_list A named list of integers.
#'   The names are the names of the clones.
#'   The integers are the index of cells that belong to the corresponding clones.
#' @param select_proteins A character vector. Specifies the proteins to plot.
#' To see available proteins, see \code{\link{get_protein_ids}}
#' @param geom Either \code{ridges} or \code{boxplot}. Default is \code{ridges}.
#' @return A named list of ggplot objects
#' @examples
#' #Compare CD34 expression between two clones.
#' #clone1: IDH1_R132H single mutant
#' #clone2: IDH1_R132H & NPM1c double mutant
#'
#'
#' \dontrun{
#' idh1_het <- get_cells_index(h5f, variants = c(IDH1_R132H = "chr2:209113112:C/T"), int_ngt = "1")
#' npm1c <- get_cells_index(h5f, variants = c(NPM1c = "chr5:170837543:C/CTCTG"), int_ngt = "1")
#' npm1_wt <- get_cells_index(h5f, variants = c(NPM1c = "chr5:170837543:C/CTCTG"), int_ngt = "0")
#'
#' index_clone1 <- intersect(idh1_het, npm1_wt)
#' index_clone2 <- intersect(idh1_het, npm1c)
#' index_cells_list <- list(clone1 = index_clone1, clone2 = index_clone2)
#'
#' plot_protein_by_clones(h5f, index_cells_list, select_proteins = "CD34", geom = "ridges")
#' }
plot_protein_by_clones <- function(h5f,
                                   index_cells_list,
                                   select_proteins = NULL,
                                   geom = c("ridges", "boxplot")){

  x <- read_protein_counts(h5f, transpose = TRUE, index_cells = NULL, normalization = "clr")

  #Specifies protein
  if(!is.null(select_proteins)){
    #check if selected proteins are available
    protein_check <- !(select_proteins %in% colnames(x))
    if(any(protein_check)){stop("selected protein is not present in data")}

    x <- x %>% as.data.frame() %>% dplyr::select(select_proteins)
  }

  ###
  group_names <- names(index_cells_list)
  dt <- data.frame()
  for(i in group_names){
    dt_sub <- x[index_cells_list[[i]], ] %>% as.data.frame()
    dt_sub$group <- i

    dt <- dplyr::bind_rows(dt, dt_sub)
  }
  ###

  ###
  #Make data tidy for plotting
  dt <- dt %>%
    pivot_longer(cols = 1:(ncol(.)-1),
                 names_to = "protein",
                 values_to = "norm_count")
  ###

  ###
  #declare list to store plots
  plot_list <- list()

  #plot
  if(geom == "ridges"){
    for(i in select_proteins){
      dt_sub <- subset(dt, protein == i)
      p <- ggplot(dt_sub, aes(norm_count, group)) +
        geom_density_ridges() +
        ggtitle(i)

      plot_list[[i]] <- p
    }
  }

  if(geom == "boxplot"){
    for(i in select_proteins){
      dt_sub <- subset(dt, protein == i)
      p <- ggplot(dt_sub, aes(group, norm_count)) +
        geom_boxplot() +
        ggtitle(i)

      plot_list[[i]] <- p
    }
  }
  ###
  return(plot_list)
}





#' Plot VAF of variants across clones
#'
#' An in-house function used to identify multiplets in sequencing.
#'
#' @param h5f h5f
#' @param index_cells_list A named list of integers.
#'   The names are the names of the clones.
#'   The integers are the index of cells that belong to the corresponding clones.
#' @param indx_variants See \code{\link{get_variants_index}}
#' @param return_data Boolean. If \code{TRUE}, returns the tidy-data used for plotting
#' @return A named list of ggplot objects, or the underlying tidy-data
#' @examples
#' #Compare IDH1_R132H VAF between wild type and mutant cells
#'
#' \dontrun{
#' idh1mt <- get_cells_index(h5f, variants = c(IDH1_R132H = "chr2:209113112:C/T"), int_ngt = 1)
#' idh1wt <- get_cells_index(h5f, variants = c(IDH1_R132H = "chr2:209113112:C/T"), int_ngt = 0)
#' index_cells_list <- list(idh1wt = idh1wt, idh1mt = idh1mt)
#'
#' plot_vaf_by_clones(h5f,
#'   index_cells_list = indx_ls,
#'   index_variants = get_variants_index(h5f, c(IDH1_R132H = "chr2:209113112:C/T"))
#'   )
#' }
plot_vaf_by_clones <- function(h5f,
                               index_cells_list,
                               index_variants,
                               return_data = FALSE){

  x <- read_assays_variants(h5f,
                            c("AF"),
                            index_variants = index_variants,
                            format = "list")

  dt_vaf <- x$AF %>% t() %>% as.data.frame()
  #dt_ngt <- x$NGT %>% t() %>% as.data.frame()
  rm(x)

  ###
  group_names <- names(index_cells_list)
  dt <- data.frame()
  for(i in group_names){
    dt_sub <- dt_vaf[index_cells_list[[i]], , drop = FALSE]
    dt_sub$group <- i

    dt <- dplyr::bind_rows(dt, dt_sub)
  }
  ###

  ###
  #Make data tidy for plotting
  dt <- dt %>%
    pivot_longer(cols = 1:(ncol(.)-1),
                 names_to = "variants",
                 values_to = "vaf")
  ###

  if(return_data == TRUE){
    return(dt)
    }

  #declare list to store plots
  plot_list <- list()

  for(i in names(index_variants)){
    dt_sub <- subset(dt, variants == i)
    p <- ggplot(dt_sub, aes(group, vaf)) +
      geom_boxplot() +
      ggtitle(i)

    plot_list[[i]] <- p
  }
  return(plot_list)
}
