% Generated by roxygen2 (4.0.1): do not edit by hand
\name{mewGetMean}
\alias{mewGetMean}
\title{Extract MEW average value}
\usage{
mewGetMean(av)
}
\arguments{
\item{av}{The current state of the MEW average}
}
\value{
(vector double length n_xx) the current value of the MEW
average if it is up-to-date
}
\description{
Return the current value of the moving expanding
window (MEW) average if it is up-to-date; otherwise, raise an error
}
\examples{
## see the examples for the function \\code{mewMean}
}

