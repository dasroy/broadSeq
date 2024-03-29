% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{round_df}
\alias{round_df}
\title{Applies round function only on numeric columns of a data.frame.}
\usage{
round_df(df, digits)
}
\arguments{
\item{df}{data.frame object}

\item{digits}{passed to \code{\link{round}}}
}
\value{
data.frame object
}
\description{
Applies round function only on numeric columns of a data.frame.
}
\examples{
data("iris")
iris \%>\% round_df(digits = 0) \%>\% head()
}
