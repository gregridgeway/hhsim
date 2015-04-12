
# set.seed(11172)

library(hhsim)

simulparam <- expand.grid(truth=c("Gaussian","Mdp","Tdist"),
                          assumed=c("Gaussian","Mdp","Tdist","NPML","SBR"),
                          rgt=c(0.1,0.33,1,10),
                          rsr=c(1,5,10),
                          n=100,
                          df=5)
i <- with(simulparam, order(n,rgt,rsr,truth,assumed))
simulparam <- simulparam[i,]
simulparam$fname <- rep("",nrow(simulparam))
for(j in 1:nrow(simulparam))
{
   simulparam$fname[j] <-
      with(simulparam[j,],
         paste(paste(truth,assumed,"_",rgt,"_",rsr,"_",n,sep="")))
}


GenerateGaussian <- function(n,rgt,rsr)
{
   sigma2max = rgt*rsr
   sigma2min = sigma2max/(rsr^2)
   sigma2 = exp(seq(from=log(sigma2min),to=log(sigma2max),length=n))
   sd = sqrt(sigma2)
   tau = sqrt(1)
   mn = 0
   data.obj = list()
   data.obj$theta = rnorm(n,mn,tau)
   data.obj$Y = rnorm(n,data.obj$theta,sd)
   data.obj$sigma2 = sigma2
   return(data.obj)
}

GenerateT <- function(n,rgt,rsr,nu)
{
   sigma2max = rgt*rsr
   sigma2min = sigma2max/(rsr^2)
   sigma2 = exp(seq(from=log(sigma2min),to=log(sigma2max),length=n))
   sd = sqrt(sigma2)
   mn = 0
   data.obj = list()
   data.obj$theta = rt(n,nu)*sqrt((nu-2)/nu)
   data.obj$Y = rnorm(n,data.obj$theta,sd)
   data.obj$sigma2 = sigma2
   return(data.obj)
}

GenerateSmix <- function(n,rgt,rsr)
{
   sigma2max = rgt*rsr
   sigma2min = sigma2max/(rsr^2)
   sigma2 = exp(seq(from=log(sigma2min),to=log(sigma2max),length=n))
   sd = sqrt(sigma2)
   # generate 75 pct from a t(8) and 25 pct from a t(50)
   data.obj = list()
   nnorm = floor(0.25*n)
   data.obj$theta = rt(n,8)
   data.obj$theta[1:nnorm] = rt(nnorm,50)
   data.obj$Y = rnorm(n,data.obj$theta,sd)
   data.obj$sigma2 =  sigma2
   return(data.obj)
}

GenerateMdp=function(n,rgt,rsr)
{
   data.obj = list()
   data.obj=NULL
   delta=4
   ups = 1 # want components to have equal variance
   eps=0.2
# simulate a mixture of 2 normals with mean 0 and var 1 with this param-ization
   a = sqrt((1-eps) + eps*ups^2 + eps*(1-eps)*delta^2)
   ind = runif(n) < (1-eps)
   data.obj$theta=ind*rnorm(n,-eps*delta/a,sqrt(1/a^2)) + (1-ind)*rnorm(n,(1-eps)*delta/a,sqrt(ups^2/a^2))
# simulate true distribution for quantiles etc
   ind = runif(n) < (1-eps)
   data.obj$theta=ind*rnorm(n,-eps*delta/a,sqrt(1/a^2)) + (1-ind)*rnorm(n,(1-eps)*delta/a,sqrt(ups^2/a^2))
   tailp = c(0.05,0.1,0.25,0.75,0.9,0.95)
   sigma2max=rgt*rsr
   sigma2min=sigma2max/(rsr^2)
   sigma2=exp(seq(from=log(sigma2min),to=log(sigma2max),length=n))
   sd=sqrt(sigma2)
   data.obj$Y = rnorm(n,data.obj$theta,sd)
   data.obj$sigma2 = sigma2
#   ind = runif(1000) < (1-eps)
#   quant=ind*rnorm(1000,-eps*delta/a,sqrt(1/a^2)) + (1-ind)*rnorm(1000,(1-eps)*delta/a,sqrt(ups^2/a^2))
   data.obj$priquan = mdppriquan
   return(data.obj)
}


# get true cdf of thetas from mixture, via simulation
# these pars will change dep on what changes in GenerateMdp function

tailp = c(0.05,0.1,0.25,0.75,0.9,0.95)
delta=4
ups = 1 # want components to have equal variance
eps=0.2
a = sqrt((1-eps) + eps*ups^2 + eps*(1-eps)*delta^2)
#ind = runif(2000000) < (1-eps)
#quant=ind*rnorm(2000000,-eps*delta/a,sqrt(1/a^2)) + (1-ind)*rnorm(2000000,(1-eps)*delta/a,sqrt(ups^2/a^2))
ind = runif(20000) < (1-eps)
quant=ind*rnorm(20000,-eps*delta/a,sqrt(1/a^2)) + (1-ind)*rnorm(20000,(1-eps)*delta/a,sqrt(ups^2/a^2))

mdppriquan = quantile(quant,tailp)

######################################################

LL=length(simulparam$truth)

for(LLind in 1:LL)
{
  ##########################################################
  ### initial values and required parameters for simulation
  evartheta=NULL
  nsimu = 500    # number of Monte Carlo iterations to perform
  nmcmc = 500  # number of MCMC iterations per each model fit
  nburn = 100  # burn-in iterations to discard at each round
  nmcmcall = nmcmc + nburn
  nmcm=as.vector(as.integer(nmcmc)) # need for R to C conversion
  kmax = simulparam$n[LLind] # sample size to generate
  alow = -6
  ahigh = 6
  alowhis = -6
  ahighhis= 6
  nquan = 6
  jmax = 3 #  number of estimators to consider (bayes, gr, mle)
  weight = rep(1,kmax) # weight (not used yet)
  tailp = c(0.05,0.1,0.25,0.75,0.9,0.95)
  nhis= 50 # number of bars in average histogram
  simurankest=rep(0,jmax)

  # initialize key vectors before simulation study
  kindex=c(0:(kmax-1))
  p=(2*(kindex+1) - 1)/(2*kmax)
  irstar=rep(0,kmax)
  simumean=simuvar=simuisel=simussel=simuwsse=rep(0,jmax)
  simusel=rep(0,jmax*kmax)
  simuwsel=rep(0,jmax*kmax)
  simuselrank=0
  simutail=rep(0,nquan*jmax)
  simuhis = rep(0,nhis*jmax)
  avgwsel=avgsel = rep(0,jmax*kmax)
  avgmean=avgvar=avgisel=weff=eff=avgwssel=avgssel = rep(0,jmax)
  avghis = rep(0,nhis*jmax)
  avgtail = rep(0,nquan*jmax)
  niter = 50

  genname = paste("Generate",simulparam$truth[LLind],sep="")
  anlname = paste("Gaussian",simulparam$assumed[LLind],sep="")

  for(nmccount in 1:nsimu)
  {
      cat(genname,anlname,nmccount,"\n")
      # generate data
      if(genname=="GenerateGaussian") {
         priquan = qnorm(tailp)
         data.obj = GenerateGaussian(kmax,
                                     simulparam$rgt[LLind],
                                     simulparam$rsr[LLind])
      } else
      if(genname=="GenerateTdist") {
         priquan =  qt(tailp,simulparam$df[LLind])*
                    sqrt((simulparam$df[LLind]-2)/
                    simulparam$df[LLind])
         data.obj = GenerateT(kmax,
                              simulparam$rgt[LLind],
                              simulparam$rsr[LLind],
                              simulparam$df[LLind])
      } else
      if(genname=="GenerateMdp") {
         data.obj = GenerateMdp(kmax,
                                simulparam$rgt[LLind],
                               simulparam$rsr[LLind])
         priquan = data.obj$priquan
      }

      evartheta=c(evartheta,c(var(data.obj$theta)))

      # posterior sample
      if(anlname=="GaussianGaussian")
         outp=GaussianGaussian(data.obj$Y,
                              data.obj$sigma2,
                              nmcmcall)
      else if(anlname=="GaussianTdist")
         outp=GaussianT(data.obj$Y,
                        data.obj$sigma2,
                        nDraws=nmcmcall,
                        nu=simulparam$df[LLind])
      else if(anlname=="GaussianMdp")
         outp = GaussianMDP(data.obj$Y,
                           data.obj$sigma2,
                           nDraws=nmcmcall,
                           numconfig=3,
                           m=0,      # initial value of m
#                           w=1,
#                           W=10,
#                           par1=10,
#                           par2=0.1)
                           w=1,
                           W=1,
                           par1=4,
                           par2=4)
      else if(anlname=="GaussianNPML")
      {
         EBestimate = GaussianNPML(data.obj$Y,data.obj$sigma2)
         outp = GaussianEBSimulate(data.obj$Y,
                                   data.obj$sigma2,
                                   mu=EBestimate$mu,
                                   alpha=EBestimate$alpha,
                                   size=nmcmcall)
      }
      else if(anlname=="GaussianSBR")
      {
         EBestimate = GaussianSBR(data.obj$Y,
                                  data.obj$sigma2,
                                  n.iters=round(3*log(kmax)))
         outp = GaussianEBSimulate(data.obj$Y,
                                   data.obj$sigma2,
                                   mu=EBestimate$mu,
                                   alpha=EBestimate$alpha,
                                   size=nmcmcall)
      }
# for preferring the Gaussian (Kmax=100, sigma2=1): w=1,W=10,par1=10,par2=0.1
# to prefer the MDP / mxiture: w=1,W=1,par1=4,par2=4

      postmean = data.matrix(data.frame(outp$theta))

      if(is.element(anlname,
                    c("GaussianGaussian","GaussianTdist",
                      "GaussianMdp")))
      {
         postmean = postmean[,c((nburn+1):nmcmcall)]
         pmmean=rowMeans(postmean)
      } else
      {
         pmmean = EBestimate$theta
         postmean = t(postmean[c((nburn+1):nmcmcall),])
      }

      # compute GR estimates
      grest = ensemble = probens = rep(0,kmax)
      tmp = gr.quantile(kmax,alow,ahigh,p,postmean,ensemble,probens,nmcmc,niter)
      ensemble = tmp$ensemble

      rhat = rep(0,kmax)

      for(k in 1:kmax)
      {
         for(j in 1:kmax)
         {
            tmp=ifelse(j != k,
                     sum(postmean[k,] >=postmean[j,])/length(postmean[k,]), 1)
            rhat[k] = rhat[k] + tmp
         }
      }
      irstar <- rank(rhat)
      grest <- ensemble[irstar]

      thetatrue = data.obj$theta
      ranktrue  = rank(thetatrue)

      selrank    = sum((irstar/(kmax+1) - ranktrue/(kmax+1))^2)/kmax
      selpmrank  = sum((rank(pmmean)/(kmax+1) - ranktrue/(kmax+1))^2)/kmax
      selgrrank  = sum((rank(grest)/(kmax+1) - ranktrue/(kmax+1))^2)/kmax
      selmlrank  = sum((rank(data.obj$Y)/(kmax+1) - ranktrue/(kmax+1))^2)/kmax

      est=cbind(pmmean,grest,data.obj$Y)

      rankest=c(selpmrank, selgrrank, selmlrank)
      simurankest= simurankest + rankest

      # evaluate estimates with loss functions at this iteration of sim study
      sel = wsel= rep(0,kmax*jmax)
      ssel = wssel = rep(0,jmax)

      # SSEL for individual estimation
      tmp = theta.ssel(kmax,jmax,thetatrue,est,weight,sel,ssel,wsel,wssel)
      sel = tmp$sel
      ssel = tmp$ssel
      wsel = tmp$wsel
      wssel = tmp$wssel

      tmp =  theta.simu(jmax,kmax,sel,ssel,simusel,simussel,
                        wsel,wssel,simuwsel,simuwsse)
      simusel = tmp$simusel
      simuwsel = tmp$simuwsel
      simussel = tmp$simussel
      simuwsse = tmp$simuwsse

      estvar = estmean = rep(0,jmax)
      tmp = ensmom(kmax,jmax,est,estmean,estvar)
      estmean = tmp$estmean
      estvar  = tmp$estvar

      grisel = mlisel = pmisel = 0

      pmisel  = ensrisk(kmax,thetatrue,pmmean,pmisel)$isel
      grisel  = ensrisk(kmax,thetatrue,grest,grisel)$isel
      mlisel  = ensrisk(kmax,thetatrue,data.obj$Y,mlisel)$isel

      grhis = mlhis = pmhis = rep(0,nhis)

      pmhis  = enshis(nhis,kmax,alowhis,ahighhis,pmhis,pmmean)$his
      mlhis  = enshis(nhis,kmax,alowhis,ahighhis,mlhis,data.obj$Y)$his
      grhis  = enshis(nhis,kmax,alowhis,ahighhis,grhis,grest)$his

      tail = rep(0,nquan*jmax)

      for (j in 1:jmax) {
         for (i in 1: nquan) {
            for (k in 1:kmax) {
                  if (est[(j-1)*kmax+k] <= priquan[i])
                      tail[(i-1)*jmax+j] = tail[(i-1)*jmax+j] + 1
            }
         }
      }

      for (j in 1:jmax)
         for (i in 1:nquan)
            tail[(i-1)*jmax+j] =  tail[(i-1)*jmax+j] / kmax

     allisel = c(pmisel,grisel,mlisel)


      simuisel = simuisel + allisel
      simumean = simumean + estmean
      simuvar = simuvar + estvar
      simuhis = simuhis + cbind(pmhis,grhis,mlhis)

      simutail = simutail + tail
      simuselrank = simuselrank + selrank
   }

   # average over all simulations

   avgsel = simusel/nsimu
   avgwsel = simuwsel/nsimu

   avgssel = simussel / nsimu
   avgwssel = simuwsse / nsimu

   avgisel = simuisel / nsimu
   avgmean = simumean / nsimu
   avgvar = simuvar / nsimu
   avghis = simuhis/(nsimu*kmax)
   avgtail = simutail/nsimu

   avgselrank = simuselrank / nsimu
   avgrankest = simurankest / nsimu
   outpt = NULL
   outpt$nsimulations = nsimu
   outpt$grboundary = c(alow,ahigh)
   outpt$nmcmc = nmcmc
   outpt$nburn = nburn
   outpt$numintervalhis = nhis

   avgsel = avgssel  ## important: want this in the output
   tt=c(1,4,7,10,13,16)
   tailprob = cbind(c(0.05,0.10,0.25,0.75,0.90,0.95),
                     avgtail[tt],avgtail[c(tt+1)], avgtail[c(tt+2)])

   output = rbind(avgsel,avgisel,avgmean,avgvar,avgrankest)*10000 # to better read output

   fname <- simulparam$fname[LLind]
   fname1 = paste(fname,".txt",sep="")
   fname2 = paste(fname,"_hist.txt",sep="")
   fname3 = paste(fname,"selrnk.txt",sep="")
   fname4 = paste(fname,"trueout.txt",sep="")

   avgselrank = 10000*avgselrank # just for output to be clearer

   dput(outpt,fname1)
   write.table(round(output,2),fname1,append=TRUE,quote=FALSE,
               row.names=FALSE,col.names=FALSE)
   write.table(round(tailprob,2),fname1,append=TRUE,quote=FALSE,
               row.names=FALSE,col.names=FALSE)
   write.table(avghis,fname2,quote=FALSE,
               row.names=FALSE,col.names=FALSE)
   write.table(round(avgselrank,2),fname3,quote=FALSE,
               row.names=FALSE,col.names=FALSE)

   tmp = c(LLind, mean(evartheta),c(range(evartheta)), median(evartheta) )
   tmp = matrix(tmp,nrow=1)
   write.table(round(tmp,4),fname4,quote=FALSE,
               row.names=FALSE,col.names=FALSE,append=TRUE)
}
warnings()




# create tables from saved results

corefc=function(nn)
{
   g1a=matrix(scan(nn,skip=8),ncol=4,byrow=T,quiet = TRUE)
   g1a=g1a[c(2:5),-1]
   g1b=matrix(scan(nn,skip=3),ncol=3,byrow=T,quiet = TRUE)[c(1,2,5),]
   g1=rbind(g1b,g1a)
   h1=g1
   h1[1:3,c(1:2)]=g1[1:3,c(1:2)]/g1[1:3,3]*100
   ind=c(3,1,2)
   return(h1[,ind])
}

i <- with(simulparam, order(n,rgt,rsr,truth,assumed))
simulparam <- simulparam[i,]

temp <- apply(simulparam[,c("rgt","rsr","n")],1,paste,collapse="+")
temp <- as.numeric(factor(temp,levels=unique(temp),ordered=TRUE))
simulparam$group <- temp

nind = 1
for(i in unique(simulparam$group))
{
   A <- with(subset(simulparam,group==i),
             split(fname,truth))
   A <- lapply(A,function(x){paste(x,".txt",sep="")})

   all <- NULL
   for(j in 1:length(A))
   {
      all1 <- NULL
      for(k in 1:length(A[[j]]))
      {
         temp <- corefc(A[[j]][k])
         if(k>1) temp <- temp[,-1]
         all1 <- cbind(all1,temp[-3,]) # remove selr
      }
      all <- rbind(all,all1)
   }

   print(round(all))
}

