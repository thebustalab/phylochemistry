% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public_functions.R
\name{normalize}
\alias{normalize}
\title{Normalizes a vector of numbers to a range of zero to one.}
\usage{
normalize(x, old_min = NULL, old_max = NULL, new_min = 0, new_max = 1)
}
\arguments{
\item{x}{The vector to normalize}

\item{old_min}{The minimum of the old range}

\item{old_max}{The maximum of the old range}

\item{new_min}{The minimum of the new range, defaults to 0}

\item{new_max}{The maximum of the new range, defaults to 1}
}
\description{
Normalizes a vector of numbers to a range of zero to one.
}
