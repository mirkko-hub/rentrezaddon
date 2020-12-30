# tibble entrez summary entries -------------------------------------------


#' extract_esummary_nuccore
#'
#' This function extracts rentrez esummary records from nuccore database
#'
#' @param esummary_record single record of an esummary object as returned by rentrez::entrez_summary(db = "nuccore", id) with a single id.
#'
#' @keywords esummary nuccore
#' @export
#' @importFrom magrittr %>% %$%
#' @family functions for the extraction of esummaries
#' @examples
#' esummary_1 <- rentrez::entrez_summary(db = "nuccore", id = "17227497")
#' extract_esummary_nuccore(esummary_record = esummary_1)
#' esummary_2 <- rentrez::entrez_summary(db = "nuccore", id = c("17227497", "33239452"))
#' esummary_2 %>% map_dfr(., extract_esummary_nuccore)
#' # Shorter way to say a similar thing; in this case set output of rentret::entrez_summary always to list, even for single records
#' esummary_2 %>% map_dfr(., magrittr::extract, c("uid", "caption", "title", "taxid", "slen", "biomol", "moltype", "topology", "sourcedb", "genome", "strand", "assemblyacc", "completeness"))


extract_esummary_nuccore <- function(esummary_record){
    tibble::tibble(
      uid = esummary_record$uid,
      caption = esummary_record$caption,
      title = esummary_record$title,
      taxid = esummary_record$taxid,
      slen = esummary_record$slen,
      biomol = esummary_record$biomol,
      moltype = esummary_record$moltype,
      topology = esummary_record$topology,
      sourcedb = esummary_record$sourcedb,
      genome = esummary_record$genome,
      strand = esummary_record$strand,
      assemblyacc = esummary_record$assemblyacc,
      completeness = esummary_record$completeness)
}

#' extract_esummary_protein
#'
#' This function extracts rentrez esummary records from protein database
#'
#' @param esummary_record single record of an esummary object as returned by rentrez::entrez_summary(db = "protein", id) with a single id.
#'
#' @keywords esummary protein
#' @export
#' @importFrom magrittr %>% %$%
#' @family functions for the extraction of esummaries
#' @examples
#' esummary_1 <- rentrez::entrez_summary(db = "protein", id = "119370761")
#' extract_esummary_protein(esummary_record = esummary_1)
#' esummary_2 <- rentrez::entrez_summary(db = "protein", id = c("119370761", "73919686"))
#' esummary_2 %>% map_dfr(., extract_esummary_protein)
#' # Alternatively, this is shorter:
#' esummary_2 %>% map_dfr(., magrittr::extract, c("uid", "caption", "accessionversion", "extra", "gi", "title", "taxid", "slen", "moltype", "sourcedb", "completeness", "organism"))


extract_esummary_protein <- function(esummary_record){
  # column strain made a problem - was NULL
  # maybe need to add a test, whether any column is NULL
  tibble::tibble(
    uid = esummary_record$uid,
    caption = esummary_record$caption,
    accessionversion = esummary_record$accessionversion,
    extra = esummary_record$extra,
    gi = esummary_record$gi,
    title = esummary_record$title,
    taxid = esummary_record$taxid,
    slen = esummary_record$slen,
    moltype = esummary_record$moltype,
    sourcedb = esummary_record$sourcedb,
    completeness = esummary_record$completeness,
    organism = esummary_record$organism
    )
}


#' extract_esummary_assembly
#'
#' This function extracts rentrez esummary records from assembly database
#'
#' @param esummary_record single record of an esummary object as returned by rentrez::entrez_summary(db = "assembly", id) with a single id.
#'
#' @keywords esummary assembly
#' @importFrom magrittr %>% %$%
#' @export
#' @family functions for the extraction of esummaries
#' @examples
#' esummary_1 <- rentrez::entrez_summary(db = "assembly", id = "29508")
#' extract_esummary_assembly(esummary_record = esummary_1)
#' esummary_2 <- rentrez::entrez_summary(db = "assembly", id = c("29508", "31208"))
#' esummary_2 %>% map_dfr(., extract_esummary_assembly)

extract_esummary_assembly <- function(esummary_record){
  # to be extended!
  tibble::tibble(
    uid = esummary_record$uid,
    assemblyaccession = esummary_record$assemblyaccession,
    organism = esummary_record$organism,
    speciesname = esummary_record$speciesname,
    taxid = esummary_record$taxid,
    assemblystatus = esummary_record$assemblystatus,
    coverage = esummary_record$coverage,
    refseq_category = esummary_record$refseq_category
  )
}



