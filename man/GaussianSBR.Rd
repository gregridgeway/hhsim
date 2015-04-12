\name{GaussianSBR}
\alias{GaussianSBR}

\title{ Smoothing by Roughening, Gaussian }
\description{
  Implements the smoothing by roughening (SBR) estimate of g(theta).
}
\usage{
GaussianSBR(Y, sigma2 = rep(1, length(Y)), n.iters = round(2 * log(length(Y))), n.masspoints = 200, verbose = FALSE)
}

\arguments{
  \item{Y}{ vector of data points }
  \item{sigma2}{ vector of known variances }
  \item{n.iters}{ number of SBR iterations }
  \item{n.masspoints}{ number of mass points for discrete SBR approximation }
  \item{verbose}{ indicator of whether to print progress information }
}
\details{
See Shen and Louis (1999) for complete details.

Assumes a model of the form \eqn{Y_k~N(\theta_k,\sigma_k)}{Y[k]~N(theta[k],sigma[k])} where \eqn{\theta_k~g(\theta)}{theta[k]~g(theta)}. This function gets posterior means for \eqn{\theta}{theta} using an empirical Bayes estimate for \eqn{g}{g}. It sets up \code{n.masspoints} equally spaced points on the range of \code{Y}. Using the EM algorithm, it iteratively refines the estimate of \eqn{g}{g} stopping after \code{n.iters} iterations. After many 1000s of iterations this will converge to the NPML for \eqn{g}{g}.
}
\value{
  \item{theta}{empirical Bayes posterior mean estimates of theta}
  \item{sigma2}{estimates of sigma2 (currently not estimated, assumed known)}
  \item{mu}{mass points for the discrete approximation to SBR}
  \item{alpha}{probability mass for each mu}
}
\references{ Laird NM, Louis TA (1991). Smoothing the non-parametric estimate of a prior distribution by roughening: An empirical study. Comput. Statist. and Data Analysis, 12:27-38.

Shen W, Louis TA (1999). Empirical Bayes Estimation via the Smoothing by Roughening
Approach. J. Computational and Graphical Statistics, 8: 800-823.}
\author{Greg Ridgeway \email{gregr@rand.org}}

\seealso{ \code{\link{GaussianNPML}}}
\examples{
k <- 100
theta <- rnorm(k,0,1)
sigma <- rep(1,k)
Y <- rnorm(k,theta,sigma)

out.sbr  <- GaussianSBR(Y,sigma^2,
                        n.iters=round(3*log(k)))
x <- with(out.sbr, sample(mu,size=100000,replace=TRUE,prob=alpha))
hist(x,main="SBR",prob=TRUE,xlab="theta",ylab="g(theta)")

plot(out.sbr$theta,theta)
}
\keyword{models}
