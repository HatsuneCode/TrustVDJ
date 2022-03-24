
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TrustVDJ

The goal of TrustVDJ is to read immune repertoire data, especially from
TRUST4, 10x Genomics cellranger or AIRR format results.

## Installation

1.  TrustVDJ is available on CRAN:

``` r
install.packages('TrustVDJ')
```

2.  Install TrustVDJ using devtools:

``` r
devtools::install_github('HatsuneCode/TrustVDJ')
```

*** Maybe dependency 'Biostrings' is not available:

``` r
install.packages('BiocManager')
BiocManager::install('Biostrings')
```


## Version

-   v0.1.0: initial release.
-   v0.1.1: fix bug in build_IMGT_reference: when same gene name,
    combine their sequences.

## Example

There are some basic examples showing how to read 10x/TRUST4 data
commonly:

``` r
library(TrustVDJ)
## basic example code
# 10x cellranger:
airr10x   = system.file('extdata', '10x_airr_rearrangement.tsv.gz', package = 'TrustVDJ')
contig10x = system.file('extdata', '10x_filtered_contig_annotations.csv.gz', package = 'TrustVDJ')
vdj10x    = Read10x(airr_file = airr10x, contig_file = contig10x, verbose = FALSE)
summary(vdj10x[,1:3])
#>    cell_id            clone_id         sequence_id       
#>  Length:40          Length:40          Length:40         
#>  Class :character   Class :character   Class :character  
#>  Mode  :character   Mode  :character   Mode  :character

# TRUST4:
airrTrust = system.file('extdata', 'TRUST4_airr.tsv.gz', package = 'TrustVDJ')
bcTrust   = system.file('extdata', 'TRUST4_barcode_report.tsv.gz', package = 'TrustVDJ')
vdjTrust  = ReadTrust(airr_file = airrTrust, barcode_report_file = bcTrust, verbose = FALSE)
summary(vdjTrust[,1:3])
#>  sequence_id          sequence           rev_comp        
#>  Length:3245        Length:3245        Length:3245       
#>  Class :character   Class :character   Class :character  
#>  Mode  :character   Mode  :character   Mode  :character
```
