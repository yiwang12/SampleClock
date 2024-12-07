% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{myQuantileNorm}
\alias{myQuantileNorm}
\title{Performs quantile normalization on a matrix of expression values}
\usage{
myQuantileNorm(x, preserve_zero = T)
}
\arguments{
\item{x}{Matrix of expression values (genes in rows, samples in columns)}

\item{preserve_zero}{Logical indicating whether to preserve zero values (default: TRUE)}
}
\value{
Matrix with same dimensions as input, containing normalized values
}
\description{
Performs quantile normalization on a matrix of expression values
}
\details{
Implementation based on Dave Tang's Blog method:
\enumerate{
\item Ranks values within each sample
\item Takes mean across samples for each rank
\item Assigns mean values back based on original ranks
\item Optionally preserves zero values from original data
}
}
\examples{
counts <- matrix(rpois(1000, lambda = 10), ncol = 10)
normalized_counts <- myQuantileNorm(counts)

}
\references{
Dave Tang's Blog: https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
}
