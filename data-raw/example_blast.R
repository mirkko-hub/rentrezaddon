## code to prepare `example_blast` dataset goes here

library(tidyverse)
library(rentrez)
library(magrittr)
library(Biostrings)
library(rentrezaddon)

# query the first 100 nuccores (cyano set) using 7120ccmm as query
nuccores <- read_csv("~/Documents/NCBI_thesis_2/analysis/20201120_assembly_information_rbcl_form.csv") %$% caption
db_path <- paste0("~/Documents/NCBI_thesis_2/nuccore/db/", nuccores)[1:100]
blastn_ccmm <- map_dfr(db_path, query_blast, query = "~/Documents/NCBI_thesis_2/queries/pcc7120_ccmm.faa", blast = tblastn)

# select both on + and - strand, each 5
df <- blastn_ccmm %>% select(sseqid, sstart, send) %>%
  mutate(df_strand = pmap(., define_strand)) %>%
  unnest(df_strand) %>%
  group_by(strand) %>%
  sample_n(5)

# create small subset that can be included in rentrezaddons
# if this is changed, REMEMBER to recreate the databases using makeblastdb!!!
dna <- paste0("~/Documents/NCBI_thesis_2/nuccore/fasta/", df$sseqid %>% unique()) %>% readDNAStringSet(.)
example_xstringset <- df %>%
  ungroup %>%
  select(sseqid, sstart, send) %>%
  blast_to_xstringset(.,  ext_3 = 500, ext_5 = 500, xstringset = dna) %>%
  `names<-`(paste0("sequence_", 1:10))

writeXStringSet(example_xstringset, filepath = "data-raw/example_xstringset.fna")
writeXStringSet(example_xstringset, filepath = "~/Documents/r_packages/test_data/example_xstringset.fna")

example_blast <- query_blast(db = "~/Documents/r_packages/test_data/example_db", query = "~/Documents/r_packages/test_data/pcc7120_rcassu.faa", blast = tblastn)
usethis::use_data(example_blast, overwrite = TRUE)
