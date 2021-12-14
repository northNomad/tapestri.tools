
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tapestri.tools

<!-- badges: start -->
<!-- badges: end -->

A toolbox for the analysis of MissionBio’s tapestri single cell
multi-omics data.

## Installation

You can install the development version of tapestri.tools from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("northNomad/tapestri.tools")
```

## Read all the DNA variants

``` r
#Load tapestri multi-omics data stored in .h5 format into session
h5f <- tapestri.tools::read_h5(path = "data_private/mpdxD-18_MB_33_38.dna+protein.h5")
```

``` r
dt_variants <- get_variants(h5f = h5f, format = "data.table")
head(dt_variants)
#>    ALT CHROM       POS             QUAL REF ado_gt_cells ado_rate
#> 1:   A  chr1 115256513   44855.19921875   G           -1       -1
#> 2:   A  chr1 115256514 597.070007324219   G           -1       -1
#> 3:   T  chr1 115256515 570.809997558594   C           -1       -1
#> 4:   G  chr1 115256516             1436   A           -1       -1
#> 5:   T  chr1 115256516             1436   A           -1       -1
#> 6:   T  chr1 115256517  1098.0400390625   C           -1       -1
#>                 amplicon filtered                 id SYMBOL
#> 1: AML_v2_NRAS_115256512       00 chr1:115256513:G/A   NRAS
#> 2: AML_v2_NRAS_115256512       01 chr1:115256514:G/A   NRAS
#> 3: AML_v2_NRAS_115256512       01 chr1:115256515:C/T   NRAS
#> 4: AML_v2_NRAS_115256512       01 chr1:115256516:A/G   NRAS
#> 5: AML_v2_NRAS_115256512       01 chr1:115256516:A/T   NRAS
#> 6: AML_v2_NRAS_115256512       01 chr1:115256517:C/T   NRAS
```

## `tapestri.tools` makes data wrangling easy

In the example experiment, we transplanted a mix of two primary AML
samples into one immunodeficient mice.

12 weeks after transplantation, we harvested the bone marrow cells and
performed single cell DNA + protein sequencing using MissionBio’s
tapestri platform. One of the AML sample carries IDH1<sup>R132H</sup>
mutation, the other carries IDH2<sup>R140Q</sup>; we know these two
*IDH* variants are mutually exclusive.

``` r
#based on hg19
int_variants <- c(IDH1_R132H = "chr2:209113112:C/T", IDH2_R140Q = "chr15:90631934:C/T")

tapestri.tools::plot_upset_with_variants(h5f, variants = int_variants, unknown = "unknown_to_ref")
#> Warning in if (format == "GRanges") {: the condition has length > 1 and only the
#> first element will be used
#> Warning in if (format == "list") {: the condition has length > 1 and only the
#> first element will be used
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

### Filter or index the cells based on the genotypes of specific variants:

``` r
#which cells carry IDH1 mutation? 
index_cells_IDH1 <- get_cells_index(h5f, variants = int_variants[1], int_ngt = c(1, 2)) 

#which cells carry IDH2 mutation?
index_cells_IDH2 <- get_cells_index(h5f, variants = int_variants[2], int_ngt = c(1, 2))

#which cells are doublets 
x <- get_cells_index(h5f, variants = int_variants, int_ngt = c(1, 2)) 
index_doublets <- intersect(x$IDH1_R132H, x$IDH2_R140Q)
```

### Read in genotype matrix of FLT3ITD variants of IDH1 mutated cells

``` r
FLT3ITD <- get_flt3itd(h5f, format = "data.table")
```

``` r
id_FLT3ITD <- FLT3ITD$id
names(id_FLT3ITD) <- paste0("FLT3ITD", 1:length(id_FLT3ITD))

read_assays_variants(h5f, 
                     included_assays = c("AF", "NGT"),
                     index_variants = get_variants_index(h5f, id_FLT3ITD),
                     index_cells = index_cells_IDH1,  
                     format = "list") -> x

names(x) #AF, NGT
#> [1] "AF"  "NGT"
dim(x$AF)
#> [1]   18 3071
```

## Protein UMAPs

``` r
get_protein_ids(h5f)
#>  [1] "CD11b"  "CD14"   "CD19"   "CD3"    "CD33"   "CD34"   "CD38"   "CD45"  
#>  [9] "CD45RA" "CD90"
```

``` r
dt_protein <- read_protein_counts(h5f, normalization = "clr", transpose = TRUE) %>% as.data.frame() 


umap_protein <- umap::umap(dt_protein)
umap_protein <- cbind(umap_protein$data, as.data.frame(umap_protein$layout))

#
ls_umap <- list()
for(i in get_protein_ids(h5f)){
  ggplot(umap_protein, aes_string("V1", "V2", color = i)) + 
    geom_point() + 
    theme_minimal() + 
    labs(x = "UMAP1", y = "UMAP2") +
    scale_color_viridis_c() -> p
  
  ls_umap[[i]] <- p
}
#
patchwork::wrap_plots(ls_umap, ncol = 2)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />
….More plotting and filtering functions to come…
