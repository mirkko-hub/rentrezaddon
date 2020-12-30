# extract cyanobacterial group from NCBI lineage --------------------------

#' extract_cyano_group
#'
#' This function extracts the cyanobacterial group from a NCBI lineage.
#'
#' @param string single single character string with lineage information of a cyanobacterial species
#'
#' @importFrom magrittr %>%
#' @export
#' @examples
#' tax_entry <- entrez_fetch(db = "taxonomy", id = 167539, rettype = "xml", parsed = TRUE)
#'  map_dfr(., function(x) tibble(taxid = x$TaxId,
#' tax_entry_tibble <- XML::xmlToList(tax_entry) %>%
#'  parent_taxid = x$ParentTaxId,
#'  lineage = x$Lineage,
#'  name = x$ScientificName))
#'
#'  # if you don't like XML just use following string
#'  # "cellular organisms; Bacteria; Terrabacteria group; Cyanobacteria/Melainabacteria group; Cyanobacteria; Synechococcales; Prochloraceae; Prochlorococcus; Prochlorococcus marinus; Prochlorococcus marinus subsp. marinus"
#'  extract_cyano_group(tax_entry_tibble$lineage)

extract_cyano_group <- function(string) {
  # extract everything between 'Cyanobactera; ' and next ';'
  tibble::tibble(tax_group = string %>% stringr::str_extract(., "(?<=Cyanobacteria\\; ).*?(?=\\; )"))
}


#' extract_tree_label
#'
#' This function extracts accessions from fasta headers
#'
#' @param string fasta header
#'
#' @importFrom magrittr %>%
#' @export
#' @examples
#'extract_tree_label()

extract_tree_label <- function(.) {
  tibble(
    # (?<=\[) - positive lookbehind for [
    # .*?  - non greedy match for the content
    # (?=\]) - positive lookahead for ]
    faa_label = .,
    faa_id = str_extract(., ".*?(?= )"),
    faa_title = str_extract(., "(?<= ).+?(?= \\[)"),
    faa_organism = str_extract(., "(?<=\\[).+?(?=\\])")
  )
}
