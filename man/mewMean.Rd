% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mew_mean.R
\name{mewMean}
\alias{mewMean}
\title{Update the moving expanding window average}
\usage{
mewMean(av)
}
\arguments{
\item{av}{(class mewTyp) the current state of the MEW average}
}
\value{
the updated instance of the argument \code{av}
}
\description{
When desired, the \code{x_mean} slot in an S4 object
of class \code{mewTyp} may be updated to contain the correct moving
expanding window (MEW) average (it is not updated by the function
\code{mewAccum} to save computation).  If the slot \code{know_mean}
is unity, the slot \code{x_mean} is up-to-date; otherwise; it is
not.
}
\examples{
n_iter <- 100

i_to_print <- 10

results <- matrix(data = double(2*n_iter/i_to_print),
                  nrow = n_iter/i_to_print,
                  ncol = 2)

av <- mewInit(n_bin = 4, n_xx = 2, ff = 0.5)

for (i in 1:n_iter) {

  value <- runif(n=2)
  value[1] <- ((cos(value[1]*2*pi))^2)*(1 - exp(-0.01*i))
  value[2] <- (-((sin(value[2]*2*pi))^2))*(1 - exp(-0.01*i))

  av <- mewAccum(xx = value, av = av)

  if (i\%\%i_to_print == 0) {

    av <- mewMean(av)
    show(av)
    results[i/i_to_print, ] <- mewGetMean(av)
  }
}

## plot the results

plot(c(1, (n_iter/i_to_print)),
     c(min(results), max(results)),
     type = "n")
points(1:(n_iter/i_to_print), results[, 1])
points(1:(n_iter/i_to_print), results[, 2])

## Now, a larger example, and we pause part way through to assess
## convergence

n_iter <- 1000
av <- mewInit(n_bin = 4, n_xx = 5000, ff = 0.5)
for (i in 1:n_iter) {

  new_samp <- runif(n = 5000)
  av <- mewAccum(xx = new_samp, av = av)
}

av <- mewMean(av = av)

## of course each element of the mean sould converge to 0.5.  After
## 1000 iterations, the first six elements of the mean vector are
show(av)

## run another 1000 iterations
for (i in 1:1000) {

  new_samp <- runif(n = 5000)
  av <- mewAccum(xx = new_samp, av = av)
}

av <- mewMean(av)

## check the mean of the first six elements again
show(av)

}
