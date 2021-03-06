% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tracks.R
\name{map2bed}
\alias{map2bed}
\title{Converts a MAP data.frame to a BED data.frame}
\usage{
map2bed(map)
}
\arguments{
\item{map}{A MAP data.frame.}
}
\value{
A BED data.frame.
}
\description{
Takes a map file and:
 \itemize{
\item{column 1: Used as the chromosome column in the BED file..}
\item{column 4: Used as start and end in the BED data.frame (as we work with SNPs).}
}
}
