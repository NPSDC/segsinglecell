% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readSegCounts.R
\name{getReadEnd}
\alias{getReadEnd}
\title{Gets the type of read end}
\usage{
getReadEnd(arg)
}
\arguments{
\item{arg}{A character/number which is the 1st entry of the second Column of the TSV file}
}
\value{
character containing either 'S' or 'P' representing single or paired end read
}
\description{
Outputs whether a file is single or paired end read by checking whether arg starts with "SEG"
}
