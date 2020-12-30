# extract genomic distances -----------------------------------------------

#' get_genomic_distance
#'
#' Extracts genomic distances between genomic elements whose positions are known!
#'
#' @param data tibble with minimum columns: assembly_accession, sstart, send, protein -
#' where sstart and send are of type doubles and refer to start and end of eg. a tblastn querry - genomic position.
#' protein and assembly_accession are of type character, assembly_accession ... to be continued
#' @param prot1 name of a protein / id
#' @param prot2 name of a protein / id
#'
#' @importFrom magrittr %>%
#' @export
#' @examples
#' get_genomic_distance()
#' protein_presence <- read_csv("/home/mirkko/Documents/NCBI/ncbi_analysis/protein_presence2.csv")
#' get_genomic_distance(protein_presence2, prot1 = "rbcl", prot2 = "rbcs") %>%
#' ggplot() +
#' geom_histogram(mapping = aes(x = distance))

get_genomic_distance <- function(data, prot1, prot2) {
  prot1 = rlang::enquo(prot1)
  prot2 = rlang::enquo(prot2)
  df <- data %>% dplyr::select(assembly_accession, sstart, send, protein)
  df_intersect <- base::intersect(df %>% dplyr::filter(., protein == !! prot1) %$% assembly_accession,
                                  df %>% dplyr::filter(., protein == !! prot2) %$% assembly_accession)
  df <- df %>% dplyr::filter(., assembly_accession %in% df_intersect)

  gr_prot1 <- df %>% dplyr::filter(., protein == !! prot1) %>%
    dplyr::mutate(strand = purrr::pmap(., rentrezaddon::define_strand)) %>%
    tidyr::unnest(strand) %>%
    dplyr::select(-c(sstart, send)) %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "assembly_accession",
                                            start.field = "start",
                                            end.field = "end",
                                            strand.field = "strand",
                                            keep.extra.columns = T)

  gr_prot2 <- df %>% dplyr::filter(., protein == !! prot2) %>%
    dplyr::mutate(strand = purrr::pmap(., rentrezaddon::define_strand)) %>%
    tidyr::unnest(strand) %>%
    dplyr::select(-c(sstart, send)) %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "assembly_accession",
                                            start.field = "start",
                                            end.field = "end",
                                            strand.field = "strand",
                                            keep.extra.columns = T)
  tibble::tibble(assembly_accession = df$assembly_accession %>% base::unique(),
                 distance = GenomicRanges::distance(gr_prot1, gr_prot2, ignore.strand = T))

}

# for now on assembly_accession; ideally on nuccore to make sure that plasmid and genome positions are not mixed up in case that both proteins are on different DNA molecule
# input: tibble with min cols - assembly_accession, sstart, ssend, protein
#         protein1, protein2
