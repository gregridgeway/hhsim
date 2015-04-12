\name{GaussianGaussian}
\alias{GaussianGaussian}
\title{Gaussian-Gaussian hierarchical model}
\description{Simulates from the posterior of a Gaussian-Gaussian model.}
\usage{
GaussianGaussian(Y,
		sigma2=rep(1,length(Y)),
                 nDraws=100,
                 s1=0.001,s2=0.001,
                 a=0,A=1000,
                 t1=0.001,t2=0.001)}
\arguments{
\item{Y}{a vector of numbers}
\item{sigma2}{vector of known \eqn{\sigma2}}
\item{nDraws}{number of MCMC draws}
\item{s1,s2}{prior parameters for \eqn{\sigma^-2}}
\item{a,A}{mean and variance of the prior for m, the mean of the thetas}
\item{t1,t2}{prior parameters for \eqn{\tau^-2}}
}

\details{No details now.}

\value{
\code{GaussianGaussian} is still being tested.
}


\author{Greg Ridgeway \email{gridge@sas.upenn.edu}}
\seealso{
}
\examples{
N <- 50
theta <- rnorm(N,0,1)
Y <- rnorm(N,theta,1)
sigma2 <- rgamma(N,10,rate=10)
niter=100 
GaussianGaussian(Y,sigma2,niter)

}
\keyword{models}
