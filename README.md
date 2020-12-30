
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rentrezaddon

**Warning: This package is still under development \!**

A package comprising functions to execute [standalone
Blast+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) searches from
within R (though executive and search path have to be set up
separately), to aid workflows based on
[rentrez](https://docs.ropensci.org/rentrez/) and
[Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html).

## Installation

You can install the released version of rentrezaddon from GitHub.

``` r
# should not work, since the repo is private
# devtools::install_github("mirkko-hub/rentrezaddon")
```

``` r
library(tidyverse, warn.conflicts = FALSE)
#> ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
#> ✓ tibble  3.0.3     ✓ dplyr   1.0.2
#> ✓ tidyr   1.1.2     ✓ stringr 1.4.0
#> ✓ readr   1.4.0     ✓ forcats 0.5.0
#> ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
library(magrittr, warn.conflicts = FALSE)
library(rentrezaddon)
```

``` r
# system.file("data-raw", "example_xstringset.fna", package = "rentrezaddon", mustWork = TRUE)
```

## Examples

#### Blast searches - Prerequisites

  - local installation of current release of [standalone
    Blast+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
  - databases generated using ‘makeblastdb’ (add ‘-parse\_seqids’ to
    work with blast\_cmd)
  - some query file in fasta format
  - search path to executives are given in the body of the function,
    which have to be modified accordingly.

*Proper setup should enable something like this to work*

(files are local only and not provided here, since everyone else would
be lacking the BLAST+ executive anyways\!)

``` r
# databases to search:
nuccores <- readr::read_csv("~/Documents/NCBI_thesis_2/analysis/20201120_assembly_information_rbcl_form.csv") %$% caption
#> 
#> ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> cols(
#>   .default = col_character(),
#>   nuccore_id = col_double(),
#>   taxid = col_double(),
#>   slen = col_double(),
#>   parent_taxid = col_double(),
#>   assembly_id = col_double()
#> )
#> ℹ Use `spec()` for the full column specifications.
db_path <- paste0("~/Documents/NCBI_thesis_2/nuccore/db/", nuccores)[1:100]

# query
ccmm <- "~/Documents/NCBI_thesis_2/queries/pcc7120_ccmm.faa"

# search
blastn_ccmm <- purrr::map_dfr(db_path, query_blast, query = ccmm, blast = tblastn)
head(blastn_ccmm)
#> # A tibble: 6 x 13
#>   qseqid sseqid pident length mismatch gapopen qstart  qend sstart   send
#>   <chr>  <chr>   <dbl>  <int>    <int>   <int>  <int> <int>  <int>  <int>
#> 1 BAB72… NC_00…  100      555        0       0      1   555 9.96e5 9.94e5
#> 2 BAB72… NC_00…   55.7     88       39       0    356   443 1.80e6 1.80e6
#> 3 BAB72… NC_00…   59.8     82       33       0    474   555 1.80e6 1.80e6
#> 4 BAB72… NC_00…   51.2    129       48       2    191   317 1.80e6 1.80e6
#> 5 BAB72… NC_00…   45.3    552      288       4      6   555 2.18e5 2.17e5
#> 6 BAB72… NC_00…   44.4    338      170       2    228   555 2.17e5 2.16e5
#> # … with 3 more variables: evalue <dbl>, bitscore <int>, database <chr>

# works with conserved domain searches, too
blast_cdd <- query_blast("/home/mirkko/Documents/NCBI/ncbi_databases/cdd/Cdd_NCBI", query = ccmm, blast = blastp)
head(blast_cdd)
#> # A tibble: 6 x 13
#>   qseqid sseqid pident length mismatch gapopen qstart  qend sstart  send
#>   <chr>  <chr>   <dbl>  <int>    <int>   <int>  <int> <int>  <int> <int>
#> 1 BAB72… CDD:1…   56.9    174       65       1     19   192      2   165
#> 2 BAB72… CDD:2…   74.4     82       21       0    234   315      3    84
#> 3 BAB72… CDD:2…   70.2     84       25       0    357   440      1    84
#> 4 BAB72… CDD:2…   65.5     84       29       0    472   555      1    84
#> 5 BAB72… CDD:1…   33.1    136       78       4     22   156      3   126
#> 6 BAB72… CDD:1…   30.9    139       78       5     21   158      2   123
#> # … with 3 more variables: evalue <dbl>, bitscore <int>, database <chr>

### check ccds with query_blastdbcmd_cdd -- see below!
```

``` r
# works only with indexed databases
# to generate an indexed daabase, you have to include '-parse_seqids' in 'makeblastdb' command!
query_blastdbcmd(db = "/home/mirkko/Documents/r_packages/test_data/example_db",
                 entry = "sequence_1",
                 dbtype = "nucl",
                 filepath = "data-raw/example_blastdbcmd", range = "1-10", strand = "minus")
Biostrings::readDNAStringSet("data-raw/example_blastdbcmd")
#> DNAStringSet object of length 7:
#>     width seq                                               names               
#> [1]    10 ATGATGGAAG                                        sequence_1:c10-1 
#> [2]    10 ATGATGGAAG                                        sequence_1:c10-1 
#> [3]    10 ATGATGGAAG                                        sequence_1:c10-1 
#> [4]    10 ATGATGGAAG                                        sequence_1:c10-1 
#> [5]    10 ATGATGGAAG                                        sequence_1:c10-1 
#> [6]    10 ATGATGGAAG                                        sequence_1:c10-1 
#> [7]    10 ATGATGGAAG                                        sequence_1:c10-1

# if no indexed database available / you don't want to generate it / above function does not work fo you
# -> try 'blast_to_xstringset()'
```

``` r
# get sequences (xstringset) from identified ranges (blast)
# generate example_xstringset and example_blast:
# example_xstringset <- system.file("data-raw", "example_xstringset.fna" , package = "rentrezaddon", mustWork = TRUE) %>% Biostrings::readDNAStringSet(.)
# example_xstringset
# example_blast
# 
# # if DNA, output is ALWAYS in 5'-3' direction no matter which DNA strand the region is on!
# example_blast %>%
# blast_to_xstringset(., xstringset = example_xstringset)
```

``` r
# get further info to cdd motifs
blast_cdd$sseqid %>% purrr::map_dfr(., query_blastdbcmd_cdd, db = "/home/mirkko/Documents/NCBI/ncbi_databases/cdd/Cdd_NCBI")
#> # A tibble: 6 x 3
#>   cdd     cdd_title                          description                        
#>   <chr>   <chr>                              <chr>                              
#> 1 CDD:10… cd00710, LbH_gamma_CA, Gamma carb… >CDD:100039 cd00710, LbH_gamma_CA,…
#> 2 CDD:23… cd00307, RuBisCO_small_like, Ribu… >CDD:238187 cd00307, RuBisCO_small…
#> 3 CDD:23… cd00307, RuBisCO_small_like, Ribu… >CDD:238187 cd00307, RuBisCO_small…
#> 4 CDD:23… cd00307, RuBisCO_small_like, Ribu… >CDD:238187 cd00307, RuBisCO_small…
#> 5 CDD:10… cd04745, LbH_paaY_like, paaY-like… >CDD:100058 cd04745, LbH_paaY_like…
#> 6 CDD:10… cd04650, LbH_FBP, Ferripyochelin … >CDD:100055 cd04650, LbH_FBP, Ferr…
```

#### Some more Biostrings helpers

  - *have a fasta/uniprot/genbank/etc. header?*
  - **just wanna rename the files with id and remove that other crap?**

why don’t you try this?

``` r
# example_header <- system.file("data-raw", "example_headers.faa" , package = "rentrezaddon", mustWork = TRUE) %>% Biostrings::readAAStringSet(.)
# example_header %>% names()
# example_header %>% names() %>% map_dfr(., extract_header_2)
# example_header %>% `names<-`(names(.) %>% map_dfr(., extract_header_2) %$% faa_id)
```

#### Some more rentrez help?

``` r
esummary_2 <- rentrez::entrez_summary(db = "protein", id = c("119370761", "73919686"))
esummary_2 %>% purrr::map_dfr(., extract_esummary_protein)
#> # A tibble: 2 x 12
#>   uid   caption accessionversion extra     gi title  taxid  slen moltype
#>   <chr> <chr>   <chr>            <chr>  <int> <chr>  <int> <int> <chr>  
#> 1 1193… Q05164  Q05164.2         gi|1… 1.19e8 RecN… 559292   967 aa     
#> 2 7391… Q9NS18  Q9NS18.1         gi|7… 7.39e7 RecN…   9606   164 aa     
#> # … with 3 more variables: sourcedb <chr>, completeness <chr>, organism <chr>
```
