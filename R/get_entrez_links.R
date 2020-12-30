# get links to entrez database entries ------------------------------------

#' get_entrez_links
#'
#' This function extracts links from an entrez elink_element for assembly - nuccore links.
#'
#' Links are extracted according to availability in a prioritized manner - in this order: refseq, insdc, other.
#' If no links are available for a given assembly, "NULL" is returned in the ouput tibble!
#'
#' @param elink_element single record of an elink_element as returned by rentrez::entrez_link(db = "nuccore", dbfrom = "assembly", id) with a single id.
#' @return tibble with nuccore_links (nuccore_id) and the source (refseq, insdc, other) where these links come from
#' @keywords elink
#' @export
#' @importFrom magrittr %>% %$%
#' @examples
#' elink <- rentrez::entrez_link(db = "nuccore", dbfrom = "assembly", id = "29508")
#' elink %>% get_entrez_links(.)
#' elink2 <- rentrez::entrez_link(db = "nuccore", dbfrom = "assembly", id = c("29508", "31208"), by_id = T )
#' elink2 %>% map_dfr(., get_entrez_links) %>% mutate(assembly_id = c("29508", "31208")) %>% unnest(.)

get_entrez_links <- function(elink_element) {
  tibble::tibble(tibble_links = base::list(
    if (elink_element %$% links %>% base::length() != 0) {
      if(checkmate::check_subset("assembly_nuccore_refseq", elink_element %$% links %>% names()) == TRUE) {
        tibble::tibble(link_type = "assembly_nuccore_refseq",
               link = elink_element %$% links %$% assembly_nuccore_refseq)
      } else if(checkmate::check_subset("assembly_nuccore_insdc", elink_element %$% links %>% names()) == TRUE) {
        tibble::tibble(link_type = "assembly_nuccore_insdc",
               link = elink_element %$% links %$% assembly_nuccore_insdc)
      } else if(checkmate::check_subset("assembly_nuccore", elink_element %$% links %>% names()) == TRUE) {
        tibble::tibble(link_type = "assembly_nuccore",
               link = elink_element %$% links %$% assembly_nuccore)
      }
    } else {
      tibble::tibble(link_type = "NULL",
             link = "NULL")
      # "No links found for specified id !"
    } ))
}


# # further idea: extend this principle to all possible link combinations - involves enquo and stuff!
# elink_1 <- rentrez::entrez_link(db = "nuccore", dbfrom = "assembly", id = "29508")
# get_entrez_links <- function(elink_element) {
#   tibble(tibble_links = list(
#     if (elink_element %$% links %>% length() != 0) {
#       if(check_subset("assembly_nuccore_refseq", elink_element %$% links %>% names()) == TRUE) {
#         tibble(link_type = "assembly_nuccore_refseq",
#                link = elink_element %$% links %$% assembly_nuccore_refseq)
#       } else {
#         tibble(link_type = "NULL",
#                link = "NULL")
#         # "No links found for specified id !"
#       }} ))
# }
# elink_1 %>% get_entrez_links(.) %>% unnest()
#
# get_entrez_links_2 <- function(elink_element, dbfrom_db) {
#   dbfrom_db = rlang::enexpr(dbfrom_db)
#
#   tibble(tibble_links = list(
#     if (elink_element %$% links %>% length() != 0) {
#       if(check_subset(x = paste0(!! dbfrom_db, "_refseq"), elink_element %$% links %>% names()) == TRUE) {
#         tibble(link_type = paste0(!! dbfrom_db, "_refseq"),
#                link = elink_element %$% links %$% assembly_nuccore_refseq)
#       } else {
#         tibble(link_type = "NULL",
#                link = "NULL")
#         # "No links found for specified id !"
#       }} ))
# }
# elink_1 %>% get_entrez_links_2(., dbfrom_db = "assembly_nuccore") %>% unnest()
#
# get_entrez_links_3 <- function(elink_element, dbfrom_db) {
#   dbfrom_db = rlang::enexpr(dbfrom_db)
#   dbfrom_db_refseq = rlang::enexpr(paste0(dbfrom_db, "_refseq"))
#
#   tibble(tibble_links = list(
#     if (elink_element %$% links %>% length() != 0) {
#       if(check_subset(x = paste0(!! dbfrom_db, "_refseq"), elink_element %$% links %>% names()) == TRUE) {
#         tibble(link_type = paste0(!! dbfrom_db, "_refseq"),
#                link = elink_element %$% links %$%  !! dbfrom_db_refseq)
#       } else {
#         tibble(link_type = "NULL",
#                link = "NULL")
#         # "No links found for specified id !"
#       }} ))
# }
#
# elink_1 %>% get_entrez_links_3(., dbfrom_db = "assembly_nuccore") %>% unnest()
#
# a <- "elephant"
# checkmate::check_subset(x = paste0(a, "3"), choices = c("elephant1", "elephant2", "elephant3") )
# b <- "assembly_nuccore"
#
# elink_1$links %$% quo(paste0(b, "_refseq")) # here is the problem !


