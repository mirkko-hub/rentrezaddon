
# Extract id from XString object ------------------------------------------

#' extract_id_from_xstringset
#'
#' Extract id from fasta header
#'
#' @param string fasta header as character vector provided by names(xstringset)
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @examples
#' AAStringSet <- readAAStringSet(filepath = "/home/mirkko/Documents/NCBI/ncbi_sequences/groel_putative.fasta", format = "fasta")
#' names(AAStringSet) <- names(AAStringSet) %>% map_dfr(., extract_id_from_xstringset) %>% pull(id)


extract_id_from_xstringset <- function(string, ...){
  # grab everything before the first blank - blank not included !
  tibble::tibble(id = string %>% stringr::str_extract(., ".*?(?= )"))
}


# Extract ncbi header -----------------------------------------------------

#' extract_ncbi_header
#'
#' Extract id, title and organism from ncbi header
#'
#' @param string fasta header as character vector provided by e.g. 'names(xstringset)' derived from ncbi
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom stringr str_extract
#' @examples
#' AAStringSet <- readAAStringSet(filepath = "/home/mirkko/Documents/NCBI/ncbi_sequences/groel_putative.fasta", format = "fasta")
#' names(AAStringSet) %>% map_dfr(., extract_ncbi_header)

extract_ncbi_header <- function(string, ...){
  # grab everything before the first blank - blank not included !
  tibble::tibble(faa_id = string %>% str_extract(., ".*?(?= )"),
                 faa_title = string %>% str_extract(., "(?<= ).+?(?= \\[)"),
                 faa_organism = string %>% str_extract(., "(?<=\\[).+?(?=\\])"))
}


# Extract uniprot header --------------------------------------------------

#' extract_uniprot_header
#'
#' Extract id, title and organism from uniprot header. Use as in 'extract_ncbi_header'
#'
#' @param string fasta header as character vector provided by 'names(xstringset)' derived from uniprot
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom stringr str_extract

extract_uniprot_header <- function(string, ...){
  # grab everything before the first blank - blank not included !
  tibble::tibble(faa_id = string %>% str_extract(., "(?<=\\|).+?(?=\\|)"),
                 faa_title = string %>% str_extract(., "(?<= ).+?(?= OS=)"),
                 faa_organism = string %>% str_extract(., "(?<=OS=).+?(?= OX=)"))
}


# Extract mixed header (gi / gb / emb) ------------------------------------

#' extract_mixed_header
#'
#' Extract id, title and organism from (gi / gb / emb) header. Use as in 'extract_ncbi_header'
#'
#' @param string fasta header as character vector provided by 'names(xstringset)' occasionally derived from ncbi
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom stringr str_detect str_extract

extract_mixed_header <- function(string, ...) {
  if(string %>% str_detect(., 'gb\\|')) {
    faa_id <-  string %>% str_extract(., "(?<=gb\\|).+?(?=\\|)")
  } else if(string %>% str_detect(., 'emb\\|')) {
    faa_id <-  string %>% str_extract(., "(?<=emb\\|).+?(?=\\|)")
  } else if (string %>% str_detect(., 'gi\\|')) {
    faa_id <-  string %>% str_extract(., "(?<=gi\\|).+?(?=\\|)")
  }
  tibble::tibble(
    faa_id = faa_id,
    # faa_pos = string %>% str_extract(., "(?<=\\:).+?(?= )"),
    faa_title = string %>% str_extract(., "(?<= ).+?(?= \\[)"),
    faa_organism = string %>% str_extract(., "(?<=\\[).+?(?=\\])"))
}


# Extract pdb or prf headers ----------------------------------------------

#' extract_pdbprf_header
#'
#' Extract id, title but not the organism from (pdb|'id'| and prf||) headers. Use as in 'extract_ncbi_header'
#'
#' @param string fasta header as character vector provided by names(xstringset) occasionally derived from ncbi
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom stringr str_detect str_extract

extract_pdbprf_header <- function(string, ...) {
  if(string %>% str_detect(., 'pdb\\|')) {
    faa_id <-  string %>% str_extract(., "(?<=pdb\\|).+?(?=\\|)")
  } else if(string %>% str_detect(., 'prf\\|\\|')) {
    faa_id <-  string %>% str_extract(., "(?<=prf\\|\\|).+?(?= )")
  }
  tibble::tibble(
    faa_id = faa_id,
    # faa_pos = string %>% str_extract(., "(?<=\\:).+?(?= )"),
    faa_title = string %>% str_extract(., "(?<= ).*"),
    faa_organism = NA)
}

# Extract uniprot/ncbi headers --------------------------------------------

#' extract_header
#'
#' Extract id, title and organism from uniprot/ncbi headers. Use as in 'extract_ncbi_header'
#'
#' @param string fasta header as character provided by names(xstringset) derived from uniprot/ncbi
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom stringr str_extract str_detect
#' @examples
#'
#' # example fasta / uniprot
#' example_header <- system.file("data-raw", "example_headers.faa" , package = "rentrezaddon", mustWork = TRUE) %>% Biostrings::readAAStringSet(.)
#'
#' example_header %>% names() %>% map_dfr(., extract_header_2)
#' example_header %>% `names<-`(names(.) %>% map_dfr(., extract_header_2) %$% faa_id)

extract_header <- function(string, ...){

  if(string %>% str_extract(., ".*?(?=\\|)") %in% c("sp", "tr")) {
    string %>% extract_uniprot_header(.)
  }
  else if (string %>% str_detect(., 'gi\\|') & string %>% str_detect(., 'gb\\|') | string %>% str_detect(., 'emb\\|')) {
    string %>% extract_mixed_header(.)
  }
  else if (string %>% str_detect(., '\\[.*?\\]') == TRUE) {
    string %>% extract_ncbi_header(.)
  }
  else {"header format not of type fasta / mixed or uniprot!"}
}


# extract_header_v2 -------------------------------------------------------

#' extract_header_2
#'
#' Extract id, title and maybe organism from uniprot/ncbi headers. Use as in 'extract_ncbi_header'
#' Should work for extract_pdbprf_headers, too - in beta.
#'
#' @param string fasta header as character provided by names(xstringset) derived from uniprot/ncbi
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom stringr str_extract str_detect
#' @examples
#'
#' # example fasta / uniprot
#' example_header <- system.file("data-raw", "example_headers.faa" , package = "rentrezaddon", mustWork = TRUE) %>% Biostrings::readAAStringSet(.)
#'
#' example_header %>% names() %>% map_dfr(., extract_header_2)
#' example_header %>% `names<-`(names(.) %>% map_dfr(., extract_header_2) %$% faa_id)

extract_header_2 <- function(string, ...){

  if (string %>% str_detect(., ".*?(?=pdb\\|)") | string %>% str_detect(., ".*?(?=prf\\|\\|)")) {
    string %>% extract_pdbprf_header(.)
  }
  else if(string %>% str_extract(., ".*?(?=\\|)") %in% c("sp", "tr")) {
    string %>% extract_uniprot_header(.)
  }
  else if (string %>% str_detect(., 'gi\\|') & string %>% str_detect(., 'gb\\|') | string %>% str_detect(., 'emb\\|')) {
    string %>% extract_mixed_header(.)
  }
  else if (string %>% str_detect(., '\\[.*?\\]') == TRUE) {
    string %>% extract_ncbi_header(.)
  }
  else {"header format not of type fasta / mixed / pdb / prf or uniprot!"}
}

# define_strand -----------------------------------------------------------

#' define_strand
#'
#' Helper function for \code{\link{blast_to_xstring}}
#'
#' @param sstart sequence start
#' @param send sequence end
#' @param ... for use with pmap()
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom checkmate assert check_count
#' @examples
#'
#' define_strand(sstart = 10, send = 20)
#' define_strand(sstart = 20, send = 10)
#' define_strand(sstart = 20, send = 20)
#'
#' # for rowwise application in a tibble:
#' example_blast %>% dplyr::transmute(sseqid, df_strand = purrr::pmap(., define_strand)) %>% tidyr::unnest(df_strand)
#'
#'
#'
#' @return Output is a tibble with columns 'start', 'end' and 'strand' -
#' strand being either '+' (5' to 3') or '-' (3' to 5')

define_strand <- function(sstart, send, ...) {

  assert(
    check_count(sstart, null.ok = FALSE, positive = TRUE),
    check_count(send, null.ok = FALSE, positive = TRUE),
    combine = "and"
  )
  if (send == sstart) {
    stop(" 'sstart == send' is not supported ")
  }
  else if (sstart > send) {
    dplyr::tibble(start = send,
                  end = sstart,
                  strand = '-')
  }
  else if (send > sstart) {
    dplyr::tibble(start = sstart,
                  end = send,
                  strand = '+')
  }
}


# define_extraction_borders -----------------------------------------------

#' define_extraction_borders
#'
#' Helper function for \code{\link{blast_to_xstring}}
#'
#' Checks whether specified extraction margins are in accordance with the width of the XString
#' and returns either the specified range or the maximus available range according to the width of the XString.
#' The specified margins 'ext_3' and 'ext_5' are interpreted in a 'strand-sensitive' manner.
#'
#' @param sseqid the name of a single object of the xstringset
#' @param start a positive integer
#' @param end a positive integer > start
#' @param strand either '+' (5' to 3') or '-' (3' to 5')
#' @param xstringset object of type XStringSet (BioStrings package)
#' @param ext_3 defines extraction limits; positive number of type count that extends the intervall in 3' direction.
#' @param ext_5 defines extraction limits; positive number of type count that extends the intervall in 3' direction
#' @param ... for use with pmap()
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom checkmate assert check_numeric check_choice check_count check_subset
#' @importFrom Biostrings width
#' @importFrom dplyr case_when
#' @importFrom tibble tibble
#' @return Output is a tibble with additional columns 'start_extr', 'end_extr', possible extraction margins as specified in the input
#' @examples
#'
#' # generate example_xstringset and example_blast:
#' example_xstringset <- system.file("data-raw", "example_xstringset.fna" , package = "rentrezaddon", mustWork = TRUE) %>% Biostrings::readDNAStringSet(.)
#' example_blast
#'
#' # define extraction borders
#' define_extraction_borders(sseqid = "sequence_1",  start = 19, end = 43, strand = "-",
#' xstringset = example_xstringset, ext_5 = 20, ext_3 = 15)
#'
#' # for rowwise application in a tibble with pmap - fixed 5' and 3' extensions for all:
#' example_blast %>%
#' dplyr::transmute(sseqid, df_strand = purrr::pmap(., define_strand)) %>%
#' tidyr::unnest(df_strand) %>%
#' dplyr::mutate(extraction = purrr::pmap(., define_extraction_borders, ext_5 = 10, ext_3 = 20, xstringset = example_xstringset)) %>%
#' unnest(extraction)
#'
#' # for rowwise application in a tibble with pmap - every hit/row with individual 5' or 3' extensions:
#' example_blast %>%
#' dplyr::transmute(sseqid, df_strand = purrr::pmap(., define_strand)) %>%
#' tidyr::unnest(df_strand) %>%
#' dplyr::mutate(ext_5 = as.double(1:17), ext_3 = as.double(17:1)) %>%
#' dplyr::mutate(extraction = purrr::pmap(., define_extraction_borders, xstringset = example_xstringset)) %>%
#' unnest(extraction)

define_extraction_borders <- function(sseqid, start, end, strand, xstringset, ext_5 = 0, ext_3 = 0, ...) {

  assert(
    check_choice(class(xstringset)[1], c("DNAStringSet", "AAStringSet")),
    # 'start' and 'end' values must be within the width of the specific contig
    check_numeric(start, upper = width(xstringset[xstringset@ranges@NAMES == sseqid])),
    check_count(start, null.ok = F, positive = T),
    check_numeric(end, upper = width(xstringset[xstringset@ranges@NAMES == sseqid])),
    check_count(end, null.ok = F, positive = T),
    check_count(ext_5, null.ok = T),
    check_count(ext_3, null.ok = T),
    check_choice(strand, c("+", "-")),
    # sseqid must be a name of one of the elements of the xstringset
    check_subset(sseqid, choices = names(xstringset), empty.ok = F),
    combine = "and"
  )

  if (end <= start) {
    stop(" 'start >= end' is not supported")
  }
  if (strand == "+") {
    left = ext_5
    right = ext_3}
  else if (strand == "-") {
    right = ext_5
    left = ext_3}

  tibble(start_extr = case_when(start >= left ~ start - left,
                                start < left ~ 1),
         end_extr = case_when(width(xstringset[xstringset@ranges@NAMES == sseqid]) >= end + right ~ end + right,
                              width(xstringset[xstringset@ranges@NAMES == sseqid]) < end + right ~ width(xstringset[xstringset@ranges@NAMES == sseqid])
                              %>% as.double()) )
}


# extract_subseq ----------------------------------------------------------

#' extract_subseq
#'
#' Helper function for \code{\link{blast_to_xstring}}
#'
#' Extracts a subsequence of a XString object with margins defined by start_extr
#' and end_extr, eg. as determined by \code{\link{define_extraction_borders}}.
#' If 'strand == "-" ', the subsequnece is returned as reverseComplement().
#'
#'
#' @param start_extr a positive integer
#' @param end_extr a positive integer > start_extr
#' @param strand either '+' (5' to 3') or '-' (3' to 5')
#' @param xstringset object of type XStringSet (BioStrings package)
#' @param ... for use with pmap()
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @examples
#'
#' # generate example_xstringset and example_blast:
#' example_xstringset <- system.file("data-raw", "example_xstringset.fna" , package = "rentrezaddon", mustWork = TRUE) %>% Biostrings::readDNAStringSet(.)
#' example_blast
#'
#' example_blast %>%
#' dplyr::transmute(sseqid, df_strand = purrr::pmap(., define_strand)) %>%
#' tidyr::unnest(df_strand) %>%
#' dplyr::mutate(extraction = purrr::pmap(., define_extraction_borders, ext_5 = 10, ext_3 = 20, xstringset = example_xstringset)) %>%
#' unnest(extraction) %>%
#' dplyr::mutate(seq = purrr::pmap(., extract_subseq, xstringset = example_xstringset)) %$% seq %>%
#' Biostrings::DNAStringSetList(.) %>% unlist(.)

extract_subseq <- function(sseqid, start_extr, end_extr, strand, xstringset, ...) {

  assert(
    check_subset(sseqid, choices = names(xstringset), empty.ok = F),
    check_choice(strand, choices = c("+", "-")),
    check_choice(class(xstringset)[1], c("DNAStringSet", "AAStringSet")),
    combine = "and"
  )

  # check for valid limits?
  # no, subseq does this anyways!
  subsequence <- xstringset[xstringset@ranges@NAMES == sseqid] %>% Biostrings::subseq(., start = start_extr, end = end_extr)

  if (class(xstringset)[1] == "AAStringSet") {
    subsequence
  }
  else if(class(xstringset)[1] == "DNAStringSet") {
    if (strand == "-") {
      subsequence <- subsequence %>% Biostrings::reverseComplement(.)
      subsequence
    } else {
      subsequence
    }
  }
}

# blast_to_xstring -----------------------------------------------------

#' blast_to_xstring
#'
#' Extract DNA or AA sequences from a specified XStringSet on the basis of a
#' blast search and specify extraction margins.
#'
#' @param df_blast a tibble presenting results from a blast search, eg. as
#' provided by \code{\link{query_blast}}. Minimal columns that need to be
#' present are: 'sseqid', 'sstart' and 'send'.
#' @param ext_3 defines extraction limits; positive number of type count that
#' extends the intervall in 3' direction.
#' @param ext_5 defines extraction limits; positive number of type count that
#' extends the intervall in 5' direction.
#' @param xstringset an object of type XStringSet (BioStrings package). It
#' contains the sequences that correspond to the ncbi database, that has been
#' searched.
#'
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom purrr pmap map_dfr
#' @importFrom dplyr mutate select transmute
#' @importFrom tidyr unnest
#' @examples
#'
#' # generate example_xstringset and example_blast:
#' example_xstringset <- system.file("data-raw", "example_xstringset.fna" , package = "rentrezaddon", mustWork = TRUE) %>% Biostrings::readDNAStringSet(.)
#' example_blast
#'
#' example_blast %>%
#' blast_to_xstringset(., xstringset = example_xstringset)
#'
#' tibble(sseqid = "sequence_1", sstart = 10, send = 5) %>% blast_to_xstringset(., ext_3 = 0, ext_5 = 0, example_xstringset)
#' tibble(sseqid = "sequence_1", sstart = 5, send = 10) %>% blast_to_xstringset(., ext_3 = 0, ext_5 = 0, example_xstringset)

blast_to_xstringset <- function(df_blast, ext_3 = 0, ext_5 = 0, xstringset, ...) {

  if (intersect(df_blast$sseqid, names(xstringset)) %>% length() != 0) {
  } else {
    names(xstringset) <- names(xstringset) %>% map_dfr(., extract_id_from_xstringset) %$% id
  }

  seq_list <- df_blast %>%
    # define strand
    transmute(sseqid, df_strand = pmap(., define_strand)) %>%
    unnest(df_strand) %>%
    # define extraction borders with respect to hit position (strand) and xstring limits
    mutate(extraction = pmap(., define_extraction_borders, ext_3 = ext_3, ext_5 = ext_5, xstringset = xstringset)) %>%
    unnest(extraction) %>%
    # extract sequences within extraction limits and return subsequences in 5' to 3' direction
    # hits on ' strand == "-" ' are given as reverseComplement()
    mutate(seq = pmap(., extract_subseq, xstringset = xstringset)) %$%
    seq

  if (class(xstringset)[1] == "DNAStringSet") {
    seq_list %>% Biostrings::DNAStringSetList(.) %>% unlist(.)
  }
  else {
    seq_list %>% Biostrings::AAStringSetList(.) %>% unlist(.)
  }
}



# tibble(sseqid = names(AAStringSet)[1:5], sstart = c(23, 243, 32, 76, 54), send = c(62, 332, 87, 59, 100)) %>%
# mutate(strand = pmap(., define_strand)) %>%
# select(sseqid, strand) %>%
# unnest(strand) %>%
# mutate(extraction = pmap(., define_extraction_borders, ext_3 = 150, ext_5 = 150, xstringset = AAStringSet)) %>%
# unnest(extraction) %>%
# mutate(seq = pmap(., extract_subseq, xstringset = AAStringSet)) %>%
# pull(seq) %>%
# AAStringSetList(.) %>% unlist(.)
#
#
# cyanophora_aa <- readAAStringSet(filepath = "/home/mirkko/Documents/NCBI/cyanophora/Cyanophora_paradoxa_MAKER_gene_predictions-022111-aa.fasta")
# cyanophora_db <- query_blast(db = "/home/mirkko/Documents/NCBI/cyanophora/cyanophora",
#                  query = "/home/mirkko/Documents/NCBI/cyanophora/query_all/query",
#                  blast = blastp)
# blast_to_xstringset(cyanophora_db, ext_3 = 150, ext_5 = 150, xstringset = cyanophora_aa)
#
# nostoc_fna <- readDNAStringSet(filepath = "/home/mirkko/Documents/NCBI/ncbi_nuccore/fasta/17227497.fna")
# nostoc_db <- query_blast(db = "/home/mirkko/Documents/NCBI/ncbi_nuccore/db_nuccore/17227497",
#                  query = "/home/mirkko/Documents/NCBI/cyanophora/query_all/query",
#                  blast = tblastn)
#
# blast_to_xstringset(nostoc_db, ext_3 = 150, ext_5 = 150, xstringset = nostoc_fna) %>% writeXStringSet(filepath = "/home/mirkko/Desktop/tmp", format = "fasta")
#


# WriteXStringSetById -----------------------------------------------------

#' writeXStringSetById
#'
#' This function writes separate .fasta files from a XStringSet object to a specified path.
#'
#' @param identifier identifier that matches the identifier in the XStringSet object
#' @param xstringset XStringSet object from package 'Biostrings' (should contain the respective ids)
#' @param path path to preexisting directory; location of the generated files
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @examples
#' writeXStringSetById()

writeXStringSetById <- function(identifier, xstringset, path) {
  writeXStringSet(xstringset[xstringset@ranges@NAMES == identifier],
                  filepath = paste0(path, identifier))
}
