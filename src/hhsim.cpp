// hhsim by Greg Ridgeway and Susan Paddock

#include <cmath>
#include "hhobject.h"
#include "gaussiangaussian.h"
#include "gaussiant.h"
//#include "gaussiansmix.h"
#include "gaussianmdp.h"


extern "C" {

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

typedef double* PDOUBLE;

SEXP test
(
    SEXP raY       // outcome or response
)
{
    unsigned long hr = 0;
    unsigned long i = 0;
    unsigned long cY = length(raY);

    SEXP rAns = NULL;

    // allocate the return vector
    PROTECT(rAns = allocVector(REALSXP, 10));
    for(i=0; i < 10; i++)
    {
       REAL(rAns)[i] = REAL(raY)[i];
    }
    UNPROTECT(1);

Cleanup:
    return rAns;
Error:
    goto Cleanup;
} // end test




SEXP GaussianGaussian
(
    SEXP radY,            // data
    SEXP radSigma2,       // known sigma2 vector
    SEXP radPriorParams,  // prior parameters s1, s2, a, A, t1, t2
    SEXP rcDraws          // number of MCMC draws
)
{
   unsigned long hr = HH_OK;

   SEXP rAns = NULL;
   SEXP rTheta = NULL;
   SEXP rTheta0 = NULL;
   SEXP rSigma = NULL;
   SEXP rSigma0 = NULL;
   SEXP rTau2 = NULL;
   SEXP rm = NULL;

   CGaussianGaussian gg;

   GetRNGstate();

   unsigned long i = 0;
   unsigned long iDraw = 0;
   unsigned long cY = length(radY);

   hr = gg.initialize(REAL(radY),
                      REAL(radSigma2),
                      length(radY),
                      REAL(radPriorParams));
   if(HH_FAILED(hr))
   {
      goto Error;
   }

   // setup the return value
   PROTECT(rAns = allocVector(VECSXP, 4));

   // allocate the thetas, a list of vectors
   PROTECT(rTheta = allocVector(VECSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,0,rTheta);
   UNPROTECT(1); // rTheta
   // allocate the Sigmas, a list of vectors
   PROTECT(rSigma = allocVector(VECSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,1,rSigma);
   UNPROTECT(1); // rSigma
   // allocate the tau2, a vector
   PROTECT(rTau2 = allocVector(REALSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,2,rTau2);
   UNPROTECT(1); // rTau2
   // allocate the m, a vector
   PROTECT(rm = allocVector(REALSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,3,rm);
   UNPROTECT(1); // rm

   for(iDraw=0; iDraw < *INTEGER(rcDraws); iDraw++)
   {
      hr = gg.sample();
      if(HH_FAILED(hr))
      {
         goto Error;
      }

      PROTECT(rTheta0 = allocVector(REALSXP, cY));
      SET_VECTOR_ELT(rTheta,iDraw,rTheta0);
      UNPROTECT(1); // rTheta0

      hr = gg.getTheta(REAL(rTheta0));
      if(HH_FAILED(hr))
      {
         goto Error;
      }

      PROTECT(rSigma0 = allocVector(REALSXP, cY));
      SET_VECTOR_ELT(rSigma,iDraw,rSigma0);
      UNPROTECT(1); // rSigma0

      hr = gg.getSigma(REAL(rSigma0));
      if(HH_FAILED(hr))
      {
         goto Error;
      }

      REAL(rTau2)[iDraw] = gg.dTau2;
      REAL(rm)[iDraw] = gg.dm;
   }

   UNPROTECT(1); // rAns

Cleanup:
   PutRNGstate();
   return rAns;
Error:
   goto Cleanup;
} // end GaussianGaussian



SEXP GaussianT
(
    SEXP radY,            // data
    SEXP radSigma2,
    SEXP radPriorParams,  // prior parameters s1, s2, a, A, t1, t2, nu
    SEXP rcDraws          // number of MCMC draws
)
{
   unsigned long hr = HH_OK;

   SEXP rAns = NULL;
   SEXP rTheta = NULL;
   SEXP rTheta0 = NULL;
   SEXP rSigma = NULL;
   SEXP rSigma0 = NULL;
   SEXP rTau2 = NULL;
   SEXP rm = NULL;

   CGaussianT gt;

   unsigned long i = 0;
   unsigned long iDraw = 0;
   unsigned long cY = length(radY);

   GetRNGstate();

   // initialize the GaussianT object
   hr = gt.initialize(REAL(radY),REAL(radSigma2),
                      length(radY),
                      REAL(radPriorParams));
   if(HH_FAILED(hr))
   {
      goto Error;
   }

   // setup the return value
   PROTECT(rAns = allocVector(VECSXP, 4));

   // allocate the thetas, a list of vectors
   PROTECT(rTheta = allocVector(VECSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,0,rTheta);
   UNPROTECT(1); // rTheta
   // allocate the Sigmas, a list of vectors
   PROTECT(rSigma = allocVector(VECSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,1,rSigma);
   UNPROTECT(1); // rSigma
   // allocate the tau2, a vector
   PROTECT(rTau2 = allocVector(REALSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,2,rTau2);
   UNPROTECT(1); // rTau2
   // allocate the m, a vector
   PROTECT(rm = allocVector(REALSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,3,rm);
   UNPROTECT(1); // rm

   for(iDraw=0; iDraw < *INTEGER(rcDraws); iDraw++)
   {
      hr = gt.sample();
      if(HH_FAILED(hr))
      {
         goto Error;
      }

      // transfer theta to R return structure
      PROTECT(rTheta0 = allocVector(REALSXP, cY));
      SET_VECTOR_ELT(rTheta,iDraw,rTheta0);
      UNPROTECT(1); // rTheta0

      hr = gt.getTheta(REAL(rTheta0));
      if(HH_FAILED(hr))
      {
         goto Error;
      }

      // transfer sigma to R return structure
      PROTECT(rSigma0 = allocVector(REALSXP, cY));
      SET_VECTOR_ELT(rSigma,iDraw,rSigma0);
      UNPROTECT(1); // rSigma0

      hr = gt.getSigma(REAL(rSigma0));
      if(HH_FAILED(hr))
      {
         goto Error;
      }

      REAL(rTau2)[iDraw] = gt.dTau2;
      REAL(rm)[iDraw] = gt.dm;
   }

   UNPROTECT(1); // rAns

Cleanup:
   PutRNGstate();
   return rAns;
Error:
   goto Cleanup;
}



SEXP GaussianMDP
(
   SEXP radY,      // data
   SEXP radSigma2, // known sigma2 vector
   SEXP riK,
   SEXP raiConfig,
   SEXP radPhi,
   SEXP radPriorParams, // prior parameters m, s2, par1, par2, alpha, w, W
   SEXP rcDraws         // number of MCMC draws
)
{
   unsigned long hr = HH_OK;

   SEXP rAns = NULL;
   SEXP rTheta = NULL;
   SEXP rTheta0 = NULL;

   SEXP rm = NULL;
   SEXP rs2 = NULL;
   SEXP ralpha = NULL;
   SEXP rnconfig = NULL;

   CGaussianMDP gmdp;

   unsigned long i = 0;
   unsigned long iDraw = 0;
   unsigned long cY = length(radY);

   GetRNGstate();

   hr = gmdp.initialize(REAL(radY),
                        REAL(radSigma2),
                        cY,
                        INTEGER(riK)[0],
                        INTEGER(raiConfig),
                        REAL(radPhi),
                        REAL(radPriorParams));
   if(HH_FAILED(hr))
   {
      goto Error;
   }

   // setup the return value

   PROTECT(rAns = allocVector(VECSXP, 5));

   // allocate the thetas, a list of vectors
   PROTECT(rTheta = allocVector(VECSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,0,rTheta);
   UNPROTECT(1); // rTheta
   // allocate the m, a vector
   PROTECT(rm = allocVector(REALSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,1,rm);
   UNPROTECT(1); // rm
   // allocate the s2, a vector
   PROTECT(rs2 = allocVector(REALSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,2,rs2);
   UNPROTECT(1); // rs2
   // allocate the alpha, a vector
   PROTECT(ralpha = allocVector(REALSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,3,ralpha);
   UNPROTECT(1); // ralpha
   PROTECT(rnconfig = allocVector(INTSXP, *INTEGER(rcDraws)));
   SET_VECTOR_ELT(rAns,4,rnconfig);
   UNPROTECT(1); // rnconfig


   for(iDraw=0; iDraw < *INTEGER(rcDraws); iDraw++)
   {
      hr = gmdp.sample();
      if(HH_FAILED(hr))
      {
         goto Error;
      }

      PROTECT(rTheta0 = allocVector(REALSXP, cY));
      SET_VECTOR_ELT(rTheta, iDraw, rTheta0);
      UNPROTECT(1); // rTheta0

      hr = gmdp.getTheta(REAL(rTheta0));
      if(HH_FAILED(hr))
      {
         goto Error;
      }

      REAL(rm)[iDraw] = gmdp.m;
      REAL(rs2)[iDraw] = gmdp.s2;
      REAL(ralpha)[iDraw] = gmdp.alpha;
      INTEGER(rnconfig)[iDraw] = gmdp.numconfig;

      // DEBUG also grab the number of configurations
   }

   UNPROTECT(1); // rAns

Cleanup:
   PutRNGstate();

   return rAns;
Error:
   goto Cleanup;
} // end GaussianMDP




void cumprob
(
   double dT,
   double *adVec,
   int iM,
   double &dProb
)
{
   int i = 0;

   while ( (i<iM) && (adVec[i] <= dT))  i++;
   dProb = (double)i/(double)iM;
}


/* --------------------------------------------------------------------- */
/*     GR estimation                                                     */
/* ********************************************************************* */
/*     Get quantiles - for ensemble estimation                           */
/* ***********************************************************************/
void grquantl
(
   int *piKmax,
   double *pdAlow,
   double *pdAhigh,
   double *adP,
   double *adPostmean,
   double *adEnsemble,
   double *adProbens,
   int *piNmcmc,
   int *piNiter
)
{
  /* this function looks OK */
   int i, j, k;
   double *adLow = NULL;
   double *adHigh = NULL;
   double *adProbmid = NULL;
   double *adMid = NULL;
   double dTemp;
   double **aadPmsort = NULL;

   adLow = new double[*piKmax];  /* assuming kmax=correct grid size */
   if(adLow==NULL)
   {
      Rprintf("memory alloc failed adLow\n");
      goto Error;
   }
   adHigh = new double[*piKmax];
   if(adHigh==NULL)
   {
      Rprintf("memory alloc failed adHigh\n");
      goto Error;
   }
   adProbmid = new double[*piKmax];
   if(adProbmid==NULL)
   {
      Rprintf("memory alloc failed adProbmid\n");
      goto Error;
   }
   adMid = new double[*piKmax];
   if(adMid==NULL)
   {
      Rprintf("memory alloc failed adMid\n");
      goto Error;
   }
   aadPmsort = new PDOUBLE[*piKmax];
   if(aadPmsort==NULL)
   {
      Rprintf("memory alloc failed aadPmsort\n");
      goto Error;
   }

/* sort mcmc output for each postmean */
   for(i = 0; i<*piKmax; i++)
   {
      aadPmsort[i] = new double[*piNmcmc];
      if(aadPmsort[i]==NULL)
      {
         Rprintf("memory alloc failed aadPmsort[%d]\n",i);
         goto Error;
      }
      for(j=0; j< (*piNmcmc); j++)
      {
         aadPmsort[i][j]=adPostmean[i*  (*piNmcmc)+j];
      }
      R_rsort(aadPmsort[i], *piNmcmc);
   }

/* *****Boundary as initial estimates uniformly ***********C */
   for (k = 0; k < (*piKmax); k++)
   {
      adLow[k] = *pdAlow;
      adHigh[k] = *pdAhigh;
   }

/* *****Iteration by interval halfing *********************C */
   for (k = 0; k <  (*piKmax); k++)
   {
      for (i = 0; i < (*piNiter); i++)
      {
         adMid[k] = (adLow[k] + adHigh[k]) / 2.;
         /* *****probability at midpoint of the interval ***********C */
         adProbmid[k] = 0.;
         for (j = 0; j <  (*piKmax); j++)
         {
            /* cumprob vs normal cdf with pmvar */
            cumprob(adMid[k], aadPmsort[j], (*piNmcmc), dTemp);

            adProbmid[k] += dTemp;
         }
         adProbmid[k] /= (double) (*piKmax);

         /* *****Cut interval into half ****************************C */
         if (adProbmid[k] <= adP[k])
         {
            adLow[k] = adMid[k];
         }
         else
         {
            adHigh[k] = adMid[k];
         }
      }
   }
   for (k = 0; k < (*piKmax); k++)
   {
      adEnsemble[k] = adMid[k];
      adProbens[k] = adProbmid[k];
   }


Cleanup:
   for(i=0; i<(*piKmax); i++)
   {
      if(aadPmsort[i]!=NULL)
      {
         delete [] aadPmsort[i];
         aadPmsort[i]=NULL;
      }
   }

   if(adLow!=NULL)
   {
       delete [] adLow;
       adLow=NULL;
   }
   if(adMid!=NULL)
   {
      delete [] adMid;
      adMid=NULL;
   }
   if(adHigh!=NULL)
   {
      delete [] adHigh;
      adHigh=NULL;
   }
   if(adProbmid!=NULL)
   {
      delete [] adProbmid;
      adProbmid=NULL;
   }
   if(aadPmsort!=NULL)
   {
      delete [] aadPmsort;
      aadPmsort=NULL;
   }

   return;

Error:
   goto Cleanup;
}


/* **********************************************************************C */
/*     Bayes risk (SEL and SSEL) for all estimates at each simulation   C */
/* **********************************************************************C */
void thetssel
(
   int *kmax,
   int *jmax,
   double *theta,
   double *est,
   double *weight,
   double *sel,
   double *ssel,
   double *wsel,
   double *wssel
)
{
   int j, k;
   /* est: jmax rows, kmax columns
      wsel, sel: jmax rows, kmax columns */

   for (j = 0; j <  jmax[0]; j++)
   {
      for (k = 0; k < kmax[0]; k++) {
         sel[j*kmax[0] + k] =
            (est[j*kmax[0]+k] - theta[k]) * (est[j*kmax[0]+k] - theta[k]);
         wsel[j*kmax[0]+k] = sel[j*kmax[0]+k] * weight[k];
      }
   }

   /* *****SSEL for each estimator **************************/
   for (j = 0; j < jmax[0]; j++)
   {
      ssel[j] = 0.0;
      wssel[j] = 0.0;
   }

   for (j = 0; j < jmax[0]; j++)
   {
      for (k = 0; k < kmax[0]; k++)
      {
         ssel[j] += sel[j*kmax[0]+k] / kmax[0];
         wssel[j] += wsel[j*kmax[0]+k];
      }
   }
}

/* **********************************************************************C */
/*     SEL and SSEL over simulations                                    C */
/* **********************************************************************C */
void thetsimu
(
   int *jmax,
   int *kmax,
   double *sel,
   double *ssel,
   double *simusel,
   double *simussel,
   double *wsel,
   double *wssel,
   double *simuwsel,
   double *simuwsse
)
{
   int j, k;

   for (j = 0; j < jmax[0]; j++)
   {
      for (k = 0; k < kmax[0] ; k++)
      {
         simusel[j*kmax[0]+k] += sel[j*kmax[0]+k];
         simuwsel[j*kmax[0]+k] += wsel[j*kmax[0]+k];
      }
   }
   for (j = 0; j <  jmax[0]; j++)
   {
      simussel[j] += ssel[j];
      simuwsse[j] += wssel[j];
   }
}




/* **********************************************************************C */
/*     Ensemble estimation: ISEL                                        C */
/* **********************************************************************C */
void ensrisk
(
   int *kmax,
   double *theta,
   double *est,
   double *isel
)
{
   double diff1, diff2, temp1, temp2;
   int k, k1, k2,kmax2;
   double epsilon,d1;
   double *thetast = NULL;
   double *sortest = NULL;
   double *comb = NULL;
   double *enssel = NULL;

   kmax2=kmax[0]*2;
   thetast = new double[kmax[0]];
   if(thetast==NULL)
   {
      Rprintf("memory alloc failed thetast\n");
      goto Error;
   }
   sortest = new double[kmax[0]];
   if(sortest==NULL)
   {
      Rprintf("memory alloc failed sortest\n");
      goto Error;
   }
   comb = new double[kmax2];
   if(comb==NULL)
   {
      Rprintf("memory alloc failed comb\n");
      goto Error;
   }
   enssel = new double[kmax2];
   if(enssel==NULL)
   {
      Rprintf("memory alloc failed enssel\n");
      goto Error;
   }

   for (k = 0; k < kmax[0]; k++)
   {
      thetast[k] = theta[k];
   }
   R_rsort(thetast, kmax[0]);

   for (k = 0; k < kmax[0]; k++) {
   sortest[k] = est[k];
   }
   R_rsort(sortest, kmax[0]);

   /* *****Create a combined vector of theta and estimate ***C */
   for (k = 0; k < kmax2; k++)
   {
      comb[k] = (k < kmax[0]) ? thetast[k] : sortest[k - kmax[0]];
   }
   /* *****Sort the combined vector *************************C */
   R_rsort(comb, kmax2);
   /* *****Keep track of k1, k2, and squared difference *****C */
   /* ***** k1 - true thetas; k2 - estimates ****************C */
   epsilon = 1e-6;
   k1 = 0;
   k2 = 0;

   for (k = 0; k < kmax2; k++)
   {
      // Get the current theta or estimate value

      temp1 = (k1 < kmax[0]) ? thetast[k1] : thetast[kmax[0]-1];
      temp2 = (k2 < kmax[0]) ? sortest[k2] : sortest[kmax[0]-1];

      // Calculate the difference between the pulled vector: check code
      d1 = comb[k] - temp1;
      diff1 = abs(d1);
      d1 = comb[k] - temp2;
      diff2 = abs(d1);
      // Identify where the current value from
      if (diff1 <= epsilon) k1++;
      if (diff2 <= epsilon) k2++;
      // Calculate vertical difference DEBUG integer division problems?
      enssel[k] = (k1 + 0. - k2) * (k1 + 0. - k2) / (kmax[0] * kmax[0]);
   }

   /* *****Compute ISEL: note sel(kmax2) is always zero *********** */
   isel[0] = 0.;
   for (k = 0; k < kmax2-1;k++)
   {
      isel[0] += (comb[k + 1] - comb[k]) * enssel[k];
   }

Cleanup:
   if(thetast != NULL) delete [] thetast;
   if(sortest != NULL) delete [] sortest;
   if(comb != NULL) delete [] comb;
   if(enssel != NULL) delete [] enssel;

   return;

Error:
   goto Cleanup;
}




} // end extern C



