/* header for random.c */

#ifndef RANDOM_H
#define RANDOM_H

#ifndef EXTERN_C
    #ifdef __cplusplus
        #define EXTERN_C    extern "C"
        #include <cstdlib> 
        #include <cmath>
        #include <ctime>
        #include <cassert>
        #include <cfloat>
    #else
        #define EXTERN_C    extern
        #include <stdlib.h> 
        #include <math.h>
        #include <time.h>
        #include <assert.h>
        #include <float.h>
    #endif
#endif // EXTERN_C


EXTERN_C double pnorm(double dX);
EXTERN_C double qnorm(double dP);

EXTERN_C double qtruncnorm(double dP,
                           double dL,
                           double dR);
EXTERN_C double rtruncnorm(double dL,
                           double dR);

EXTERN_C double         runi();
EXTERN_C double         runi2();
EXTERN_C double         runir(double dA,double dB);
EXTERN_C long           runii(long iA,long iB); 
EXTERN_C double         rnorm(double dMu,double dSd);
EXTERN_C double         rcnorm(double dMu,double dSd);
EXTERN_C double         rigaus(double dMu, double dLambda);
EXTERN_C unsigned long  rpois(double dMu);
EXTERN_C double         rexp(double dLambda);
EXTERN_C double         rcauchy(double dLoc,double dScale);
EXTERN_C double         rcauchy_icdf(double dLoc,double dScale);
EXTERN_C double         rgamma(double dAlpha);
EXTERN_C double         rgamma2(double dAlpha,double dBeta);
EXTERN_C double			rgamma_rs(double dAlpha);
EXTERN_C double         rbeta(double dAlpha,double dBeta);
EXTERN_C int            rdirichlet(double *pdP, double *pdAlpha,unsigned long k);
EXTERN_C int            rdanish(double *adP, double *adAlpha, double dGamma, int k);
EXTERN_C int            rmultnom(unsigned long *aiCount, 
                                 unsigned long n, 
                                 double *adPr, 
                                 unsigned long k);
EXTERN_C unsigned long  rmultnom1(double *adPr, unsigned long k);
EXTERN_C double         rstudent(double dT);
EXTERN_C double         rlogistic(double dLoc,double dScale);
EXTERN_C double         rlnorm(double dNorm, double dNorsd);
EXTERN_C int            rbin(int n,double p);
EXTERN_C double         rweibull(double dGamma);
EXTERN_C double         rf(double dT, double dU);
EXTERN_C double         rchisq(double dT); 
EXTERN_C int            rperm(unsigned long *piResult, int n);
EXTERN_C int            rperm_some(unsigned long *piResult, 
                                   unsigned long cLength, 
                                   unsigned long cPermute);

EXTERN_C void           setrandomseed(int s);
EXTERN_C void           randomseed();

#endif // RANDOM_H
