# extract features from feature table -------------------------------------

#' featuretable_to_tibble
#'
#' This function produces a tidy tibble of output from 'rentrez::entrez_fetch(db = "protein", rettype = "ft", retmode = "text")'
#'
#' @param ft output from 'rentrez::entrez_fetch(db = "protein", rettype = "ft", retmode = "text")'
#' @export
#' @importFrom magrittr %>% %$%
#' @examples
#' # NOT WORKING (PROPERLY)!
#' rentrez::entrez_fetch(db = "protein", id = "BAB77899.1", rettype = "ft", retmode = "text") %>% featuretable_to_tibble(.)


featuretable_to_tibble <- function(ft){
  orig_wd <- base::getwd()
  base::setwd("/home/mirkko")
  base::cat(ft,  file = "../../tmp/ft.txt")
  df <- readr::read_delim(file = "../../tmp/ft.txt",
                   delim = "\t", skip_empty_rows = F, col_names = c("start", "end", "feature", "modifier", "value")) %>%
    dplyr::mutate(feature = case_when(str_detect(start, '>Feature') ~ start)) %>%
    dplyr::select(feature, start, end, modifier, value) %>%
    tidyr::fill(fill_cols = "feature")
  df <- df[- which(df$start %>% str_detect(., '>Feature')),] %>%
    tidyr::fill(fill_cols = c("start", "end")) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(accession_version = feature %>% str_extract(., "(?<=\\|).+?(?=\\|)")) %>%
    dplyr::select(accession_version, everything())
  base::system2("rm", args = "../../tmp/ft.txt")
  base::on.exit(base::setwd(orig_wd))
  df
}
