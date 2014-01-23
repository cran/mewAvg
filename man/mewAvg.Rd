\name{mewAvg}
\alias{mewAvg}
\title{A Fixed Memeory Moving Expanding Window Average}
\usage{
  mewAvg(f, n.bin, n.xx, ff, n.save, n.iter, i.to.save,
    ...)
}
\arguments{
  \item{f}{The user defined R function.  See the 'Details'
  section for more on defining this function.}

  \item{n.bin}{The fixed number of bins to use to define
  the moving expanding window.}

  \item{n.xx}{The length of the numeric vector returned by
  \code{f}.}

  \item{ff}{The fraction of the samples to included in each
  window.}

  \item{n.save}{The number of estimates to save and
  return.}

  \item{n.iter}{The number of times to call \code{f}.}

  \item{i.to.save}{A vector of zeros and ones of length
  \code{n.iter} where position \code{i} is 1 if an average
  should be calculated and saved at iteration i, and zero
  otherwise.}

  \item{...}{The initial arguments to \code{f}.}
}
\value{
  A matrix of dimension \code{n.save} by \code{n.xx}
  containing the saved averages
}
\description{
  Computes the average of a sequence of random vectors in a
  moving expanding window using a fixed amount of storage
}
\details{
  The function \code{f} should generate the sequence of
  radom vectors one at a time.  The returned value from a
  single call should be a list with at least one element.
  The first element should be a numeric vector of length
  \code{n.xx} (the next vector in the sequence), and the
  remaining elements should be the updated arguments for
  the next call to \code{f}.  The 'Examples' section
  provides further guidance.
}
\examples{
MyFun <- function (k) {

 value <- runif(n=2)
 value[1] <- ((cos(value[1]*2*pi))^2)*(1 - exp(-0.01*k))
 value[2] <- (-((sin(value[2]*2*pi))^2))*(1 - exp(-0.01*k))

 k <- k + 1

 return(list(value=value, k=k))
}

i.to.save <- seq(from=1, to=1025, by=32)
tmp <- rep(x=0, times=1025)
tmp[i.to.save] <- 1
i.to.save <- tmp

mean.vals <- mewAvg(f=MyFun,
                    n.bin=4,
                    n.xx=2,
                    ff=0.5,
                    n.save=sum(i.to.save),
                    n.iter=length(i.to.save),
                    i.to.save=i.to.save,
                    k=1)

plot(c(1:sum(i.to.save),
       1:sum(i.to.save)),
     c(mean.vals[, 1],
       mean.vals[, 2]),
     type="n",
     xlab="Saved Iter",
     ylab="Mean")
points(1:sum(i.to.save),
       mean.vals[, 1])
points(1:sum(i.to.save),
       mean.vals[, 2])

## an AR(1) process

ArOne <- function (x.old, phi, sig.eps) {

  value <- phi*x.old + rnorm(n=1, mean=0, sd=sig.eps)

  return(list(value=value, x.old=value))
}

mean.vals.ar1 <- mewAvg(f=ArOne,
                        n.bin=4,
                        n.xx=1,
                        ff=0.5,
                        n.save=sum(i.to.save),
                        n.iter=length(i.to.save),
                        i.to.save=i.to.save,
                        x.old=0,
                        phi=0.5,
                        sig.eps=1)

plot(x=c(1, sum(i.to.save)),
     y=c(-0.5, 0.5),
     xlab="Saved Iter",
     ylab="Mean",
     type="n")
points(x=1:sum(i.to.save),
       y=mean.vals.ar1)
abline(h=0, col="red")
}
\author{
  Adam L. Pintar <adam.pintar@nist.gov>

  Zachary H. Levine <zachary.levine@nist.gov>
}
\references{
  Add JSS reference when published.
}

