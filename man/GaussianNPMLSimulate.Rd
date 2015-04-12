\name{GaussianEBSimulate}
\alias{GaussianEBSimulate}

\title{ Get posterior draws for theta }
\description{
  Simulates draws from f(theta|y) using a given empirical Bayes estimate for g(theta).
}
\usage{
GaussianEBSimulate(Y, sigma2=rep(1,length(Y)),mu,alpha,size=100,verbose=FALSE)
}

\arguments{
  \item{Y}{ vector of data points }
  \item{sigma2}{ vector of known variances }
  \item{mu}{ vector of mass locations for g(theta) }
  \item{alpha}{ vector of probability mass for each mu }
  \item{size}{ number of Monte Carlo draws}
  \item{verbose}{ indicator of whether to print progress information }
}
\details{
Draws samples from f(theta|y) via Monte Carlo. Assumes that Y is normal with mean theta and variance sigma2.
}
\value{
  \item{theta}{empirical Bayes posterior mean estimates of theta}
  \item{sigma2}{estimates of sigma2 (currently not estimated, assumed known)}
  \item{mu}{mass points from the NPML}
  \item{alpha}{probability mass for each mu}
}

\author{Greg Ridgeway \email{gregr@rand.org}}

\seealso{ \code{\link{GaussianNPML}},\code{\link{GaussianSBR}}}
\examples{
k <- 100
theta <- rnorm(k,0,1)
sigma <- rep(1,k)
Y <- rnorm(k,theta,sigma)

out.npml  <- GaussianNPML(Y,sigma^2)
GaussianEBSimulate(Y,sigma^2,mu=out.npml$mu,alpha=out.npml$alpha,size=100)
}
\keyword{models}
