
#include "gaussianmdp.h"

CGaussianMDP::CGaussianMDP()
{
}

CGaussianMDP::~CGaussianMDP()
{
   if(nconfig != NULL) delete [] nconfig;
   if(prob != NULL) delete [] prob;
   if(ysamp != NULL) delete [] ysamp;
   if(sigmasamp != NULL) delete [] sigmasamp;
}


HHRESULT CGaussianMDP::getTheta
(
   double *adTheta
)
{
   HHRESULT hr = HH_OK;

   unsigned long i = 0;
   for(i=0; i<cY; i++)
   {
      adTheta[i] = theta[i];
   }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


HHRESULT CGaussianMDP::initialize
(
    double *Y,           // data
    double *sigma2,      // fixed sigma2 params
    unsigned long cY,
    int K,
    int *config,
    double *phi,
    double *adPriorParams  // prior parameters m, s2, alpha, w, W
)
{
   HHRESULT hr = HH_OK;
   fInitialized = true;

   m     = adPriorParams[0];
   s2    = adPriorParams[1];
   par1  = adPriorParams[2];
   par2  = adPriorParams[3];
   alpha = adPriorParams[4];
   w     = adPriorParams[5];
   W     = adPriorParams[6];
   
   this->Y      = Y;
   this->sigma2 = sigma2;
   this->cY     = cY;
   this->K      = K;
   this->config = config;
   this->phi    = phi;
   
   theta = NULL;
   theta = new double[cY];  
   if(theta==NULL)
   {
      Rprintf("memory alloc failed theta\n");
      goto Error;
   }
   nconfig = NULL;
   nconfig = new int[cY];
   if(nconfig==NULL)
   {
      Rprintf("memory alloc failed nconfig\n");
      goto Error;
   }
   prob = NULL;
   prob = new double[cY+1];
   if(prob==NULL)
   {
      Rprintf("memory alloc failed prob\n");
      goto Error;
   }
   ysamp = NULL;
   ysamp = new double[cY];
   if (ysamp==NULL) 
   { 
      Rprintf("memory alloc failed ysamp\n");
      goto Error;
   }
   sigmasamp = NULL;
   sigmasamp = new double[cY];
   if (sigmasamp==NULL)
   { 
      Rprintf("memory alloc failed sigmasamp\n");
      goto Error;
   }

Cleanup:
   return hr;
Error:
   goto Cleanup;
} // end initialize




HHRESULT CGaussianMDP::sample()
{
   int j = 0;
   HHRESULT hr = HH_OK;

   numconfig = 0;

   for(j=0; j<cY; j++)
   {
      hr = sample_config(config, j, sigma2, cY, Y, phi, alpha);
      numconfig = imax2(numconfig, config[j]);
   }
   
   hr = sample_phi(config, Y, sigma2, s2, m, cY, phi, numconfig);
   hr = sample_m(s2, phi, numconfig, m);
   hr = sample_s2(s2, phi, m, numconfig, w, W);
   hr = sample_theta(config, phi, cY, theta);
   hr = sample_alpha(par1, par2, cY,numconfig, alpha);
      
Cleanup:
   return hr;
Error:
   goto Cleanup;
}



/* # code for computing a hierarchical model, with normally distributed
   # level 1 errors (variance known) and level 2 follows a DP
   
   # y[i]:      observed datum for obs i
   # theta[i]:  level 1 mean for obs i
   # phi:       vector of unique values of theta (i.e., clusters)
   # config[i]: cluster label / configuration indicator

   ####################################################################
*/
HHRESULT CGaussianMDP::sample_config
(
   int *&config,
   int obs,
   double *sigma2,
   int n,
   double *y,
   double *phi,
   double alpha
)
{
  /*
    # config: vector of configuration indicators
    # obs:    index of observation under study
    # sigma2: (known) level 1 variances(
    # n:      sample size
  */
  
   int i,j,nclus,oldconfig,ind;
   int sumconfig = 0;
   double sumprob;
   double tempphi = 0.0;
   
   HHRESULT hr = HH_OK;
  
   /* get number of configurations/clusters 
      also set up other things to check */
   sumconfig = 0;
   nclus = 0; /* number of configurations */
   for(i=0; i<n; i++)
   {
      if(config[i]==config[obs]) sumconfig++;
      nclus = imax2(config[i], nclus);
   }

   /* ## STEP 1: nothing changes if obs under study (obs) has its own 
         cluster w/prob */
   if( (sumconfig == 1) && (runif(0.0,1.0) < (nclus-1.0)/nclus))
   { 
      goto Cleanup;
   }
  
   
   // nconfig counts obs in clusters, current obs not included

   for(i=0; i<nclus; i++) 
   {      
      nconfig[i] = 0; 
   }
   for(j=0; j<n; j++) 
   {
      nconfig[config[j]-1]++;
   }
   nconfig[config[obs]-1]--; /* #nclus-star */
   
   /* STEP 2: if there are more than 1  obs in case i's cluster, then: */

   if(sumconfig > 1)
   {        
      sumprob = 0;
   
      for(j=0; j<nclus; j++)
      { 
         prob[j] = nconfig[j] *
                   dnorm(y[obs], phi[j], sqrt(sigma2[obs]), 0);
         sumprob += prob[j];
      }
      prob[nclus] = (alpha/(nclus+1)) *
                    dnorm(y[obs], phi[nclus], sqrt(sigma2[obs]), 0); 
      sumprob+=prob[nclus];
      if(sumprob==0)
      { 
         for(j=0; j<=nclus; j++) prob[j]=1.0;
      }

      /* need to add in a sample-type function */
      config[obs] = multinomial(nclus+1,prob);
   
      goto Cleanup;
   }

/* STEP 3: if there is just one obs in cluster but need to sample new clustr:*/
   /*         else  s(i)=1 and need to sample new cluster */
   if(sumconfig==1)  /* # s/b unnec line */
   {
      oldconfig=config[obs];
      for(i=0; i<n; i++)
      {
         if(config[i] > oldconfig) config[i]--;
      }
      config[i]=nclus;

      for(i=1; i<nclus; i++)
      {
         if(i>=oldconfig) 
         {
            nconfig[i-1]=nconfig[i];/* last elt of nconfig now useless */
         }
      }
      
      // shifting the phis down by one, move phi[oldconfig-1] to the end
      if((oldconfig < nclus) && (nclus>1))
      { 
         tempphi = phi[oldconfig-1];
         
         for(i=oldconfig; i<nclus; i++)
         {
            phi[i-1] = phi[i];
         }
         phi[nclus-1] = tempphi;
      }
         
      nclus--;
      sumprob = 0.0;
      for(j=0; j<nclus; j++)
      {
         prob[j] = nconfig[j] * 
                   dnorm(y[obs], phi[j], sqrt(sigma2[obs]), 0);
         sumprob += prob[j];
      }
         
      prob[nclus] = (alpha/(nclus+1)) *
                    dnorm(y[obs], phi[nclus], sqrt(sigma2[obs]), 0);
      sumprob += prob[nclus];
      if(sumprob == 0) 
      {
         for(i=0; i<=nclus; i++) prob[i] = 1.0;
      }
         
      config[obs] = multinomial(nclus+1,prob);
   }
   
Cleanup:
   return hr;
   
Error:
   goto Cleanup;
}

double CGaussianMDP::sample_phi0
(
   int n,
   double *ysamp,
   double *sigmasamp,
   double s2,
   double m
)
{
   int i;
   double s=0.0;
   double ys=0.0;
   double var;
   double mn;
   
   for(i=0; i<n; i++)
   {
      s += 1.0/sigmasamp[i];
      ys += ysamp[i]/sigmasamp[i];
   }
   
   var = 1.0/(1.0/s2 + s);
   mn = (m/s2 + ys)*var;
   return rnorm(mn,sqrt(var));
}


/************************************************************/

HHRESULT CGaussianMDP::sample_phi
(
   int *config,
   double *y,
   double *sigma2,
   double s2,
   double mn,
   int n,
   double *&phi,
   int K
)
{
   int i,j,ind;
   
   HHRESULT hr = HH_OK;

   /* count how many cases have each configuration */
   for(j=0; j<K; j++) 
   {      
      nconfig[j] = 0; 
   }
   for(i=0; i<n; i++) 
   {
      nconfig[config[i]-1]++;
   }

   for(j=0; j<K; j++)
   {
      ind = 0;
      for(i=0; i<n; i++)
      {
         if(config[i]==(j+1))
         {
            ysamp[ind] = y[i];
            sigmasamp[ind] = sigma2[i];
            ind++;
         }
      }
      phi[j] = sample_phi0(nconfig[j],ysamp,sigmasamp,s2,mn);
      /* sample empty cluster phis from prior */
   }
   for(j=K;j<n;j++) 
   {
      phi[j] = rnorm(mn,sqrt(s2));
   }
   
Cleanup:
   return hr;

Error:
   goto Cleanup;   
}

/*********************************************************/

HHRESULT CGaussianMDP::sample_theta
(
   int *config,
   double *phi,
   int n, 
   double *&theta
)
{
   int i=0;
   HHRESULT hr = HH_OK;
   
   for(i=0; i<n; i++)
   { 
      theta[i] = phi[config[i]-1]; 
   }
   
   return hr;
}

/*********************************************************/

HHRESULT CGaussianMDP::sample_m
(
   double s2,
   double *phi,
   int sizephi,
   double &m
)
{
   int i;
   double mn = 0.0;
   HHRESULT hr = HH_OK;
   
   for(i=0; i<sizephi; i++)
   {
      mn += phi[i];
   }
   mn /= sizephi;
   m = rnorm(mn,sqrt(s2/sizephi));
   
   return hr;
}

/*********************************************************/
HHRESULT CGaussianMDP::sample_s2
(
   double &s2, 
   double *phi,
   double m,
   int sizephi,
   double w,
   double W
)
{
   double p1 = (double)sizephi/2.0 + w;
   double p2,s;
   int i;
   HHRESULT hr = HH_OK;
   
   s=0.0;
   for(i=0; i<sizephi; i++)
   {
      s += (phi[i]-m)*(phi[i]-m);
   }
   
   p2 = (s/2.0 + W); 
   
   s2 = 1.0/ rgamma2(p1,p2);
   
   return hr;
}


/*********************************************************/
HHRESULT CGaussianMDP::sample_alpha
(
   double par1,
   double par2,
   int n,
   int k,
   double &alpha
)
{
   double b,odds,prob;
   int ind;
   HHRESULT hr = HH_OK;
   
   b = rbeta(alpha+1,n);
   odds = (par1+k-1)/(n*(par2-log(b)));
   prob = odds/(odds+1);
   ind = (int)rbinom(1,prob);
   alpha = ind     * rgamma2(par1+k,   (par2-log(b))) +
           (1-ind) * rgamma2(par1+k-1, (par2-log(b)));
   
   return hr;
}  

/************/


int CGaussianMDP::multinomial(int ncell, double * nvec)
{
   /* draws just one from a multinomial distribution */
   int i, bindraw;
   double denom,tmp;
   
   /* draw multinomial via binomials */
   denom=0.0;
   
   for(i=0; i<ncell; i++)
      denom+=nvec[i];
   
   for(i=0; i<(ncell-1); i++) 
   {
      tmp = nvec[i]/denom;
      denom -= nvec[i];
      bindraw = runif(0.0,1.0)<=tmp;
      if(bindraw==1)
      {
         bindraw *= (i+1);
         return(bindraw);
      }
   }

   /* if 1,..,k-1 cells don't contain draw, then the last cell contains the draw*/
   bindraw = ncell;
   return(bindraw);
}





