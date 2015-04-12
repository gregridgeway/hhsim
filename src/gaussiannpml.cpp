// hhsim by Greg Ridgeway and Susan Paddock

#include <cmath>
#include "buildinfo.h"

extern "C" {

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP GaussianNPML
(
    SEXP radY,            // data
    SEXP radSigma2,       // known sigma2 vector
    SEXP radMu,           // initial point mass locations
    SEXP radAlpha,        // initial mass weights
    SEXP rdTolerance,     // convergence abs(theta-theta.old)
    SEXP rfVerbose        // flag for printing progress
)
{
   unsigned long hr = HH_OK;

   SEXP rAns = NULL;
   SEXP rTheta = NULL;
   SEXP rSigma2 = NULL;
   SEXP rMu = NULL;
   SEXP rAlpha  = NULL;

   unsigned long i  = 0;
   unsigned long k  = 0;
   double dNum      = 0.0;
   double dDen      = 0.0;

   bool fConverged     = false;
   double dThetaChange = 0.0;
   unsigned long cY    = length(radY);
   double *adY         = REAL(radY);
   double *adTheta     = NULL;
   double *adSigma2    = NULL;

   double *adP     = NULL;
   double *adAlpha = NULL;
   double *adMu    = NULL;
   double *adMuNum = NULL;
   double *adMuDen = NULL;

   if(radY==NULL)
   {
      hr = HH_INVALIDARG;
      goto Error;
   }
   if(radSigma2==NULL)
   {
      hr = HH_INVALIDARG;
      goto Error;
   }

   adP = new double[length(radMu)];
   if(adP==NULL)
   {
      hr = HH_OUTOFMEMORY;
      goto Error;
   }
   adMuNum = new double[length(radMu)];
   if(adMuNum==NULL)
   {
      hr = HH_OUTOFMEMORY;
      goto Error;
   }
   adMuDen = new double[length(radMu)];
   if(adMuDen==NULL)
   {
      hr = HH_OUTOFMEMORY;
      goto Error;
   }

   // setup the return value
   PROTECT(rAns = allocVector(VECSXP, 4));
   // allocate theta, a vector
   PROTECT(rTheta = allocVector(REALSXP, cY));
   SET_VECTOR_ELT(rAns,0,rTheta);
   UNPROTECT(1); // rTheta
   // allocate Sigma, a vector
   PROTECT(rSigma2 = allocVector(REALSXP, cY));
   SET_VECTOR_ELT(rAns,1,rSigma2);
   UNPROTECT(1); // rSigma2
   // allocate Mu, a vector, point masses for g
   PROTECT(rMu = allocVector(REALSXP, length(radMu)));
   SET_VECTOR_ELT(rAns,2,rMu);
   UNPROTECT(1); // rMu
   // allocate Alpha, a vector, mass on mu[k]
   PROTECT(rAlpha = allocVector(REALSXP, length(radAlpha)));
   SET_VECTOR_ELT(rAns,3,rAlpha);
   UNPROTECT(1); // rAlpha

   adTheta  = REAL(rTheta);
   adSigma2 = REAL(rSigma2);
   adMu     = REAL(rMu);
   adAlpha  = REAL(rAlpha);
   adY      = REAL(radY);

   // initialize mass points at y[k]
   for(i=0; i<cY; i++)
   {
      adTheta[i]  = adY[i];
      adSigma2[i] = REAL(radSigma2)[i];
   }
   for(k=0; k<length(radMu); k++)
   {
      adMu[k]     = REAL(radMu)[k];
      adAlpha[k]  = REAL(radAlpha)[k];
   }
   
   // NPML begins here
   fConverged = false;
   while(!fConverged)
   {
      dThetaChange = 0.0;
      for(k=0; k<length(radMu); k++)
      {
         adMuNum[k] = 0.0;
         adMuDen[k] = 0.0;
      }
      for(i=0; i<cY; i++)
      {
         dNum = 0.0;
         dDen = 0.0;
         for(k=0; k<length(radMu); k++)
         {
            adP[k] = adAlpha[k]*dnorm(adY[i],adMu[k],sqrt(adSigma2[i]),0);
            dNum  += adP[k]*adMu[k];
            dDen  += adP[k];
         }
         for(k=0; k<length(radMu); k++)
         {
            adP[k] /= dDen;
            adMuNum[k] += adP[k]*adY[i];
            adMuDen[k] += adP[k];
         }
         dThetaChange = fmax2(dThetaChange, fabs(adTheta[i] - dNum/dDen));
         adTheta[i] = dNum/dDen;
      }
      dDen = 0.0;
      for(k=0; k<length(radMu); k++)
      {
         dDen += adMuDen[k];
      }
      for(k=0; k<length(radMu); k++)
      {
         adMu[k]    = adMuNum[k]/adMuDen[k];
         adAlpha[k] = adMuDen[k]/dDen;
      }

      fConverged = (dThetaChange < *REAL(rdTolerance));
      if(*INTEGER(rfVerbose)==1) Rprintf("%g\n",dThetaChange);
   }

   UNPROTECT(1); // rAns

Cleanup:
   if(adP!=NULL)
   {
      delete [] adP;
      adP == NULL;
   }
   if(adMuNum!=NULL)
   {
      delete [] adMuNum;
      adMuNum == NULL;
   }
   if(adMuDen!=NULL)
   {
      delete [] adMuDen;
      adMuDen == NULL;
   }

   return rAns;
Error:
   goto Cleanup;
} // end GaussianNPML



} // end extern C



