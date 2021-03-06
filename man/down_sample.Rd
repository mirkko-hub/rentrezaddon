% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chromatography.R
\name{down_sample}
\alias{down_sample}
\title{down_sample()}
\usage{
down_sample(data, key, value, var_exclude, factor)
}
\arguments{
\item{data}{tidy long data with only one}

\item{key}{see tidyr::spread and tidyr::gather}

\item{value}{see tidyr::spread and tidyr::gather}

\item{var_exclude}{see tidyr::spread and tidyr::gather}

\item{factor}{numeric deterining the degree of downsampling}
}
\description{
This function can be used for convenient downsampling of data in long format.
}
\details{
Takes tidy, long data, spreads it, downsamples it by factor and gathers it again to tidy data!
Important: function works on redundancy reduced data, eg. if additional groupwise
information are given, introduce a unique ID for out- (before) and rejoining (after) calling down_sample.
See example below!s
}
\examples{
# test data:
tab <- tibble(Group = rep(LETTERS[1:3], each = 10), Info = rep(letters[1:3], each = 10), Time = rep(1:10, 3), Observation = rnorm(30) )
# reduce redundancy:
tab_info <- select(tab, Group, Info) \%>\% distinct()
tab_to_ds <- select(tab, -Info)
# downsample
tab_ds <- down_sample(tab_to_ds, key = "Group", value = "Observation", var_exclude = Time, factor = 2)
# rejoin
left_join(tab_ds, tab_info, by = "Group")
}
