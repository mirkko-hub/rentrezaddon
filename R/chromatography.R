# downsample crowded data for eg. plotting --------------------------------

#' down_sample()
#'
#' This function can be used for convenient downsampling of data in long format.
#'
#' Takes tidy, long data, spreads it, downsamples it by factor and gathers it again to tidy data!
#' Important: function works on redundancy reduced data, eg. if additional groupwise
#' information are given, introduce a unique ID for out- (before) and rejoining (after) calling down_sample.
#' See example below!s
#'
#' @param data tidy long data with only one
#' @param factor numeric deterining the degree of downsampling
#' @param key see tidyr::spread and tidyr::gather
#' @param value see tidyr::spread and tidyr::gather
#' @param var_exclude see tidyr::spread and tidyr::gather
#'
#' @export
#' @importFrom magrittr %>% %$%
#' @examples
#' # test data:
#' tab <- tibble(Group = rep(LETTERS[1:3], each = 10), Info = rep(letters[1:3], each = 10), Time = rep(1:10, 3), Observation = rnorm(30) )
#' # reduce redundancy:
#' tab_info <- select(tab, Group, Info) %>% distinct()
#' tab_to_ds <- select(tab, -Info)
#' # downsample
#' tab_ds <- down_sample(tab_to_ds, key = "Group", value = "Observation", var_exclude = Time, factor = 2)
#' # rejoin
#' left_join(tab_ds, tab_info, by = "Group")

down_sample <- function(data, key, value, var_exclude, factor) {
  var_exclude <- enquo(var_exclude)
  spread_data <- data %>%
    tidyr::spread(., key = key, value = value)
  spread_data[seq(1, dim(spread_data)[1], factor),]  %>%
    tidyr::gather(., key = !!key, value = !!value, -!!var_exclude)
}
