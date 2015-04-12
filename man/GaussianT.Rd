\name{GaussianT}
\alias{GaussianT}
\title{Gaussian-t hierarchical model}
\description{Simulates from the posterior of a Gaussian-t model.}
\usage{
# Susan: fix the parameters here to match the R function
GaussianT(Y,
	  sigma2=rep(1,length(Y)),
          nDraws=100,
          s1=0.001,s2=0.001,
          a=0,A=1000,
          t1=0.001,t2=0.001,nu=7)}
\arguments{
\item{Y}{a vector of numbers}
\item{sigma2}{vector of \eqn{\sigma2_i}}
\item{nDraws}{number of MCMC draws}
\item{s1,s2}{prior parameters for \eqn{\sigma^-2}}
\item{a,A}{mean and variance of the prior for m, the mean of the thetas}
\item{t1,t2}{prior parameters for \eqn{\tau^-2}}
\item{nu}{degrees of freedom of the t-distribution}
}

\details{No details now.}

\value{
\code{GaussianT} is still being tested.
}


\author{Susan Paddock \email{paddock@rand.org}}
\seealso{
\code{\link{GaussianGaussian}}
}
\examples{
N <- 20 
Y <- rnorm(N,0,1)
niter=10
sigma2 <- rep(1,N)
tmp=GaussianT(Y,sigma2,niter,nu=6)

}
\keyword{models}
