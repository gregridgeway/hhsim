
test <- function(Y)
{
   theta <- .Call("test",
                  Y=as.double(Y),
                  PACKAGE = "hhsim")

   return(theta)
}


GaussianGaussian <- function(Y,
                             sigma2=rep(1,length(Y)),
                             nDraws=100,
                             s1=0.001,s2=0.001,
                             a=0,A=1000,
                             t1=0.001,t2=0.001)
{
   results <- .Call("GaussianGaussian",
                     adY=as.double(Y),
                     adSigma2=as.double(sigma2),
                     aPriorParams=as.double(c(s1,s2,a,A,t1,t2)),
                     cDraws=as.integer(nDraws),
                     PACKAGE = "hhsim")

   names(results) <- c("theta","sigma","tau2","m")
   return(results)
}

GaussianT <- function(Y,
                      sigma2=rep(1,length(Y)),
                      nDraws=100,
                      s1=0.001,s2=0.001,  
                      a=0,A=1000,         
                      t1=0.001,t2=0.001,nu=7)
{
   results <- .Call("GaussianT",
                     adY=as.double(Y),
                     adSigma2=as.double(sigma2),
                     aPriorParams=as.double(c(s1,s2,a,A,t1,t2,nu)),
                     cDraws=as.integer(nDraws),
                     PACKAGE = "hhsim")

   names(results) <- c("theta","sigma","tau2","m")
   return(results)
}


GaussianMDP <- function(Y,
                        sigma2=rep(1,length(Y)),
                        nDraws=100,
                        numconfig=3,
                        m=0,
                        w=1,
                        W=1,
                        s2=1/rgamma(1,w,rate=W),
                        par1=8,
                        par2=2,
                        alpha=rgamma(1,par1,rate=par2))
{
   config = sample(c(1:numconfig, sample(1:numconfig,length(Y)-numconfig,replace=TRUE)))
   phi = rnorm(length(Y),m,sqrt(s2))

   results <- .Call("GaussianMDP",
                    adY           = as.double(Y),
                    adSigma2      = as.double(sigma2),
                    iK            = as.integer(numconfig),
                    aiConfig      = as.integer(config), # DEBUG get the number of configs in return
                    adPhi         = as.double(phi),
                    adPriorParams = as.double(c(m,s2,par1,par2,alpha,w,W)),
                    cDraws        = as.integer(nDraws),
                    PACKAGE = "hhsim")

   names(results) <- c("theta","m","s2","alpha","numconfig")
   return(results)
}


GaussianNPML <- function(Y,
                         sigma2=rep(1,length(Y)),
                         mu=Y,
                         alpha=rep(1/length(Y),length(Y)),
                         tolerance=0.0001,
                         verbose=FALSE)
{
   results <- .Call("GaussianNPML",
                     adY=as.double(Y),
                     adSigma2=as.double(sigma2),
                     radMu=as.double(mu),
                     radAlpha=as.double(alpha),
                     rdTolerance=as.double(tolerance),                     
                     fVerbose=as.integer(verbose),
                     PACKAGE = "hhsim")

   names(results) <- c("theta","sigma2","mu","alpha")
   return(results)
}




GaussianSBR <- function(Y,
                        sigma2=rep(1,length(Y)),
                        n.iters=round(2*log(length(Y))),
                        n.masspoints=200,
                        verbose=FALSE)
{
   mu <- seq(min(Y),max(Y),length=n.masspoints)
   results <- .Call("GaussianSBR",
                     adY         = as.double(Y),
                     adSigma2    = as.double(sigma2),
                     cIterations = as.integer(n.iters),
                     adMu        = as.double(mu),
                     fVerbose    = as.integer(verbose),
                     PACKAGE = "hhsim")

   names(results) <- c("theta","sigma2","alpha")
   results$mu <- mu

   return(results)
}


GaussianEBSimulate <- function(Y,
                         sigma2=rep(1,length(Y)),
                         mu,
                         alpha,
                         size=100,
                         verbose=FALSE)
{
   theta <- apply(cbind(Y,sigma2),1,
                  function(x,mu,alpha,size)
                  {
                     p <- dnorm(x[1],mu,sqrt(x[2]))*alpha
                     sample(mu,size=size,replace=TRUE,prob=p)
                  },
                  mu=mu,alpha=alpha,size=size)

   return(list(theta=theta))
}




gr.quantile <- function(kmax,
                        alow,
                        ahigh,
                        p,
                        postmean,
                        ensemble,
                        probens,
                        nmcmc,
                        niter)
{
   result <-
      .C("grquantl",
         kmax = as.integer(kmax),
         alow = as.double(alow),
         ahigh = as.double(ahigh),
         p = as.double(p),
         postmean = as.double(postmean),
         ensemble = as.double(ensemble),
         probens = as.double(probens),
         nmcmc = as.integer(nmcmc),
         niter = as.integer(niter),
         PACKAGE = "hhsim")

   return(result)
}


theta.ssel <- function(kmax,
                       jmax,
                       thetatrue,
                       est,
                       weight,
                       sel,
                       ssel,
                       wsel,
                       wssel)
{
   result <-
      .C("thetssel",
         kmax = as.integer(kmax),
         jmax = as.integer(jmax),
         thetatrue = as.double(thetatrue),
         est = as.double(est),
         weight = as.double(weight),
         sel = as.double(sel),
         ssel = as.double(ssel),
         wsel = as.double(wsel),
         wssel = as.double(wssel),
         PACKAGE = "hhsim")

   return(result)
}



theta.simu <- function(jmax,
                       kmax,
                       sel,
                       ssel,
                       simusel,
                       simussel,
                       wsel,
                       wssel,
                       simuwsel,
                       simuwsse)
{
   result <-
      .C("thetsimu",
         jmax = as.integer(jmax),
         kmax = as.integer(kmax),
         sel = as.double(sel),
         ssel = as.double(ssel),
         simusel = as.double(simusel),
         simussel = as.double(simussel),
         wsel = as.double(wsel),
         wssel = as.double(wssel),
         simuwsel = as.double(simuwsel),
         simuwsse = as.double(simuwsse),
         PACKAGE = "hhsim")

   return(result)
}


ensmom <- function(kmax,
                   jmax,
                   est,
                   estmean,
                   estvar)
{
   temp <- matrix(est,ncol=kmax,byrow=TRUE)

   result <-
      list(estmean = rowMeans(temp),
           estvar = apply(temp,1,var))

   if(FALSE) # DEBUG: replaced this call with the R code
   {
      result <-
         .C("ensmom",
            kmax = as.integer(kmax),
            jmax = as.integer(jmax),
            est = as.double(est),
            estmean = as.double(estmean),
            estvar = as.double(estvar),
            PACKAGE = "hhsim")
   }

   return(result)
}

ensrisk <- function(kmax,
                    thetatrue,
                    est,
                    isel)
{
   result <-
      .C("ensrisk",
         kmax = as.integer(kmax),
         thetatrue = as.double(thetatrue),
         est = as.double(est),
         isel = as.double(isel),
         PACKAGE = "hhsim")

   return(result)
}

enshis <- function(nhis,
                   kmax,
                   alowhis,
                   ahighhis,
                   pmhis,
                   pmmean)
{
   breaks <- seq(alowhis,ahighhis,length=nhis+1)
   result <- list(his = table(cut(pmmean,breaks=breaks,right=FALSE)))

   return(result)
}
