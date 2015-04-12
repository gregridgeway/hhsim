// hhsim by Greg Ridgeway and Susan Paddock

#include <cmath>
#include "buildinfo.h"

extern "C" {

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP GaussianSBR
(
    SEXP radY,            // data
    SEXP radSigma2,       // known sigma2 vector
    SEXP rcIterations,    // number of SBR iterations
    SEXP radMu,           // mass points for discrete SBR approximation
    SEXP rfVerbose        // flag for printing progress
)
{
   unsigned long hr = HH_OK;

   SEXP rAns = NULL;
   SEXP rTheta = NULL;
   SEXP rSigma2 = NULL;
   SEXP rAlpha  = NULL;

   unsigned long i    = 0;
   unsigned long k    = 0;
   unsigned long iSBR = 0;
   double dNum        = 0.0;
   double dDen        = 0.0;

   bool fConverged     = false;
   unsigned long cY    = length(radY);
   unsigned long cMu   = length(radMu);
   double *adY         = REAL(radY);
   double *adTheta     = NULL;
   double *adSigma2    = NULL;

   double *adP         = NULL;
   double *adAlpha     = NULL;
   double *adAlphaTemp = NULL;
   double *adMu        = NULL;

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

   adP = new double[cMu];
   if(adP==NULL)
   {
      hr = HH_OUTOFMEMORY;
      goto Error;
   }
   adAlphaTemp = new double[cMu];
   if(adAlphaTemp==NULL)
   {
      hr = HH_OUTOFMEMORY;
      goto Error;
   }

   // setup the return value
   PROTECT(rAns = allocVector(VECSXP, 3));
   // allocate theta, a vector
   PROTECT(rTheta = allocVector(REALSXP, cY));
   SET_VECTOR_ELT(rAns,0,rTheta);
   UNPROTECT(1); // rTheta
   // allocate Sigma, a vector
   PROTECT(rSigma2 = allocVector(REALSXP, cY));
   SET_VECTOR_ELT(rAns,1,rSigma2);
   UNPROTECT(1); // rSigma2
   // allocate Alpha, a vector, mass on mu[k]
   PROTECT(rAlpha = allocVector(REALSXP, cMu));
   SET_VECTOR_ELT(rAns,2,rAlpha);
   UNPROTECT(1); // rAlpha

   adTheta  = REAL(rTheta);
   adSigma2 = REAL(rSigma2);
   adMu     = REAL(radMu);
   adAlpha  = REAL(rAlpha);
   adY      = REAL(radY);

   // initialize mass points at y[k]
   for(i=0; i<cY; i++)
   {
      adTheta[i]  = adY[i];
      adSigma2[i] = REAL(radSigma2)[i];
   }
   for(k=0; k<cMu; k++)
   {
      adAlpha[k]  = 1.0/cMu;
   }

   // SBR begins here
   for(iSBR=0; iSBR<*INTEGER(rcIterations); iSBR++)
   {
      if(*INTEGER(rfVerbose)==1) Rprintf("%d\n",iSBR);
      for(k=0; k<cMu; k++)
      {
         adAlphaTemp[k] = 0.0;
      }
      for(i=0; i<cY; i++)
      {
         dNum = 0.0;
         dDen = 0.0;
         for(k=0; k<cMu; k++)
         {
            adP[k] = adAlpha[k]*dnorm(adY[i],adMu[k],sqrt(adSigma2[i]),0);
            dNum  += adP[k]*adMu[k];
            dDen  += adP[k];
         }
         for(k=0; k<cMu; k++)
         {
            adAlphaTemp[k] += adP[k]/dDen;
         }
         adTheta[i] = dNum/dDen;
      }

      dDen = 0.0;
      for(k=0; k<cMu; k++)
      {
         dDen += adAlphaTemp[k];
      }
      for(k=0; k<cMu; k++)
      {
         adAlpha[k] = adAlphaTemp[k]/dDen;
      }
   }

   UNPROTECT(1); // rAns

Cleanup:
   if(adP!=NULL)
   {
      delete [] adP;
      adP == NULL;
   }
   if(adAlphaTemp!=NULL)
   {
      delete [] adAlphaTemp;
      adAlphaTemp == NULL;
   }

   return rAns;
Error:
   goto Cleanup;
} // end GaussianNPML



} // end extern C



