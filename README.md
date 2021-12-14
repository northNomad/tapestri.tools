
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

## Load tapestri multi-omics data stored in .h5 format into session

``` r
h5f <- tapestri.tools::read_h5(path = "data_private/mpdxD-18_MB_33_38.dna+protein.h5")
```

## Read all the variants called by haplotype caller

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
tapestri platform.

One of the AML sample carries IDH1<sup>R132H</sup> mutation, the other
carries IDH2<sup>R140Q</sup>; we know these two *IDH* variants are
mutually exclusive.

``` r
#based on hg19
int_variants <- c(IDH1_R132H = "chr2:209113112:C/T", IDH2_R140Q = "chr15:90631934:C/T")

tapestri.tools::plot_upset_with_variants(h5f, variants = int_variants)
#> Warning in if (format == "GRanges") {: the condition has length > 1 and only the
#> first element will be used
#> Warning in if (format == "list") {: the condition has length > 1 and only the
#> first element will be used
#> Warning in if (unknown == "unknown_to_ref") {: the condition has length > 1 and
#> only the first element will be used
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
intersect(x$IDH1_R132H, x$IDH2_R140Q)
#>   [1]   16  122  155  224  246  252  256  261  291  313  319  324  372  375  405
#>  [16]  437  449  450  468  470  498  537  558  562  582  602  700  725  731  739
#>  [31]  804  815  958  972  979 1011 1031 1047 1117 1216 1223 1301 1366 1380 1383
#>  [46] 1423 1468 1516 1540 1558 1561 1567 1594 1680 1711 1748 1773 1791 1800 1832
#>  [61] 1844 1855 1882 1897 1912 1919 1923 1984 1992 1999 2025 2093 2138 2145 2174
#>  [76] 2220 2222 2266 2304 2363 2418 2421 2444 2446 2489 2496 2512 2557 2598 2600
#>  [91] 2606 2609 2631 2672 2689 2722 2807 2827 2843 2847 2855 2858 2891 2900 2919
#> [106] 2924 2932 2997 3010 3038 3044 3052 3148 3180 3189 3190 3206 3215 3233 3242
#> [121] 3278 3342 3344 3374 3400 3401 3449 3467 3543 3554 3567 3611 3624 3639 3684
#> [136] 3701 3702 3711 3734 3806 3816 3918 3952 3986 4007 4014 4027 4035 4053 4054
#> [151] 4093 4096
```

### Read in genotype matrix of FLT3ITD variants of IDH1 mutated cells

``` r
id_FLT3ITD <- get_flt3itd(h5f)$id
#> Warning in if (format == "GRanges") {: the condition has length > 1 and only the
#> first element will be used
#> Warning in if (format == "data.table") {: the condition has length > 1 and only
#> the first element will be used
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

….More plotting and filtering functions to come…
