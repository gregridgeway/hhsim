\name{GaussianNPML}
\alias{GaussianNPML}

\title{ Nonparametric Maximum Likelihood estimation, Gaussian }
\description{
  Implements the nonparametric maximum likelihood estimation (NPML) estimate of g(theta).
}
\usage{
GaussianNPML(Y,
             sigma2=rep(1,length(Y)),
             mu=Y,
             alpha=rep(1/length(Y),length(Y)),
             tolerance=0.0001,
             verbose=FALSE)
}

\arguments{
  \item{Y}{ vector of data points }
  \item{sigma2}{ vector of known variances }
  \item{mu}{ vector of initial mass points for the NPML }
  \item{alpha}{ vector of initial weights on each mass point, mu }
  \item{tolerance}{ absolute change in thetas before considered converged }
  \item{verbose}{ indicator of whether to print progress information }
}
\details{
Assumes a model of the form \eqn{Y_k~N(\theta_k,\sigma_k)}{Y[k]~N(theta[k],sigma[k])} where \eqn{\theta_k~g(\theta)}{theta[k]~g(theta)}. This function gets posterior means for \eqn{\theta}{theta} using the NPML estimate for \eqn{g}{g}. It uses the EM algorithm to solve for \eqn{g}{g}. No attempt is made to merge very close mass points since our main interest is in estimates of \eqn{\theta_k}{theta[k]}.
}
\value{
  \item{theta}{empirical Bayes posterior mean estimates of theta}
  \item{sigma2}{estimates of sigma2 (currently not estimated, assumed known)}
  \item{mu}{mass points from the NPML}
  \item{alpha}{probability mass for each mu}
}
\references{ Laird N.M. (1982). Empirical Bayes estimate using the non-parametric maximum likelihood estimate of the prior. Journal of Statistical Computation and Simulation, 15:211-220.}
\author{Greg Ridgeway \email{gregr@rand.org}}

\seealso{ \code{\link{GaussianSBR}}}
\examples{
k <- 100
theta <- rnorm(k,0,1)
sigma <- rep(1,k)
Y <- rnorm(k,theta,sigma)

out.npml  <- GaussianNPML(Y,sigma^2)
x <- with(out.npml, sample(mu,size=100000,replace=TRUE,prob=alpha))
hist(x,main="NPML",prob=TRUE,xlab="theta",ylab="g(theta)")

plot(out.npml$theta,theta)
}
\keyword{models}
