#include "random.h"

const double PI = 3.141592653589793;
const double EE = 2.718281828459046;


/*************************************/
void setrandomseed(int s)
{
    srand(s);
}


/*************************************/
void randomseed()
{
    srand((unsigned) time(NULL));
}


/*************************************/
// generates RAND_MAX unique random numbers between 0 and 1 exclusive.
inline double runi()
{
    return (rand()+1.0)/(RAND_MAX+2.0);
}


/*************************************/
// generates RAND_MAX*RAND_MAX unique random numbers between 0 and 1 exclusive.
double runi2()
{
    return (rand()*RAND_MAX+rand()+1.0)/(RAND_MAX*(RAND_MAX+1.0)+2.0);
}


/*************************************/
double runir(double dA,double dB)
{
    return (dB-dA)*runi()+dA;
}


/*************************************/
long runii(long iA,long iB)
{
    long lU = 0; 

    lU = (long)floor(runi2()*(iB-iA+1.0))+iA;

    return lU;
}


/*************************************/
double rnorm(double dMu,double dSd)
{
    static bool fCachedNormal = false;
    static double dCachedNormal = 0.0;
    double dE=0.0;
    double dV1=0.0;
    double dV2=0.0;
    double dW=0.0;
    if(!fCachedNormal)
    {
        do
        {
            dV1=2*runi()-1.0;
            dV2=2*runi()-1.0;
            dW=dV1*dV1+dV2*dV2;
        }while(dW>1.0);

        dE=sqrt((-2.0*log(dW))/dW);
        dCachedNormal=dV1*dE;
        fCachedNormal = true;;
        return dV2*dE*dSd+dMu;
    }
    else
    {
        fCachedNormal = false;
        return dCachedNormal*dSd+dMu;
    }
}


/*************************************/
double rcnorm(double dMu,double dSd)
/**************
Central limit based normal generator
Add 12 standard uniforms
**************/
{
    int i;
    double dResult=0.0;
    double dUnit=0.0;
    for(i=1; i<=12; i++)
    {
        dUnit = dUnit + runi();
    }
    dResult = dMu + dSd*( dUnit - 6.0);
    return(dResult);
}


/*************************************/
double rigaus(double dMu, double dLambda)
{
    double dU=0.0;
    double dV=0.0;
    double dW=0.0;
    double dC=0.0;
    double dX1=0.0;
    dU = rnorm(0.0,1.0);
    dV = dU*dU;
    dW = dMu*dV;
    dC = dMu / (2.0*dLambda);
    dX1 = dMu + dC*(dW - sqrt(dW * (4.0 * dLambda + dW)));
    return (runi() > dMu/(dMu + dX1) ? dMu*dMu/dX1 : dX1);
}


/*************************************/
unsigned long rpois(double dMu)
{ 
    double dP=1.0;
    unsigned long i=0;
    if(dMu<=0)
    {
        return 0;
    }
    dMu = exp(-dMu);
    do
    {
        dP = dP*runi();
        i++;
    }
    while(dP>=dMu);
    return i-1;
}


/*************************************/
double rexp(double dLambda)
{
    return -log(runi())/dLambda;
}


/*************************************/
/* Polar method                      */
double rcauchy(double dLoc,double dScale)
{
    static double dVert = 0.0;
    static double dHoriz = 0.0;

    do 
    { // tangent of random angle
        dVert = runi();
        dHoriz = 2.0*runi()-1.0;
    } while(dVert*dVert+dHoriz*dHoriz > 1.0);

    return dScale*dVert/dHoriz+dLoc;
}

/*************************************/
/* Slower inverse cdf method         */
double rcauchy_icdf(double dLoc,double dScale)
{
    return tan(runi()*PI)*dScale+dLoc;
}


/*************************************/
double rgamma(double dAlpha)
{
    static double dR1,dR2,dAA,dX,dW,dC1,dC2,dC3,dC4,dC5;
    double dReturnValue = -HUGE_VAL;
    if(dAlpha<=0.0)
    {
        dReturnValue = 0.0;
    }
    else if(dAlpha == 1.0)
    {
        dReturnValue = rexp(1.0);
    }
    else if(dAlpha<1.0)
    {
        dAA=(dAlpha+EE)/EE;
        do
        {
            dR1=runi();
            dR2=runi();
            if(dR1>1.0/dAA)
            {
                dX = -log(dAA*(1.0-dR1)/dAlpha);
                if(dR2<pow(dX,(dAlpha-1.0)))
                {
                    dReturnValue = dX;
                    break;
                }
            }
            else
            {
                dX = pow((dAA*dR1),(1.0/dAlpha));
                if(dR2<exp(-dX))
                {
                    dReturnValue = dX;
                    break;
                }
            }
        }while(dR2<2);
    }
    else  // dAlpha > 1.0
    {
        dC1=dAlpha-1.0;
        dC2=(dAlpha-1.0/(6.0*dAlpha))/dC1;
        dC3=2.0/dC1;
        dC4=dC3+2.0;
        dC5=1.0/sqrt(dAlpha);
        do
        {
            do
            {
                dR1=runi();
                dR2=runi();
                if(dAlpha>2.5)
                {
                    dR1=dR2+dC5*(1.0-1.86*dR1);
                }
            }while(dR1<=0 || dR1 >= 1);
            dW=dC2*dR2/dR1;
            if(dC3*dR1+dW+1/dW <= dC4)
            {
                dReturnValue = dC1*dW;
                break;
            }
            if(dC3*log(dR1)-log(dW)+dW<1)
            {
                dReturnValue = dC1*dW;
                break;
            }
        }while(dR2<2);
    }

    return dReturnValue;
}


double rgamma2(double dAlpha,double dBeta)
// mean = alpha/beta, var = alpha/beta^2
{
    return rgamma(dAlpha)/dBeta;
}


/*************************************/
double rbeta(double dAlpha,double dBeta)
{
   double dR1;
   if(dAlpha <=0.0 || dBeta <= 0.0) 
   {
       return 0.0;
   }
   dR1=rgamma(dAlpha);
   return dR1/(dR1+rgamma(dBeta));
}


/*************************************/
int rdirichlet(double *adP, double *adAlpha,unsigned long k)
{
    unsigned long j;
    double dTotal = 0.0;

    assert(adP!=NULL);
    assert(adAlpha!=NULL);

    for(j=0; j<k; j++)
    {
        assert(adAlpha[j]>0); 
        adP[j] = rgamma(adAlpha[j]);
        dTotal += adP[j];
    }

    for(j=0; j<k; j++)
    {
        adP[j] = adP[j]/dTotal;
    }

    return 1;
}


/*************************************/
int rdanish(double *adP, double *adAlpha, double dGamma, int k)
/*********************************************************************
For alpha[1,...,k] and  gamma = lambda*prod(alpha)^-1 (gamma is invariant)
then lambda_i = gamma*(alpha_i)^2.  Finally, y_i~ IG(alpha_i, lambda_i)
and y_i/sum(y_i) ~ danish(alpha, lambda)
*********************************************************************/
{
    int i;
    double dLambda;
    double dTotal=0.0;

    assert(adP!=NULL);
    assert(adAlpha!=NULL);

    for(i=0; i<k; i++)
    {
        assert(adAlpha[i]>0);
        dLambda = dGamma*adAlpha[i]*adAlpha[i];
        adP[i] = rigaus(adAlpha[i], dLambda);
        dTotal += adP[i];
    }
    for(i=0; i<k; i++)
    {
        adP[i] = adP[i]/dTotal;
    }

    return 1;
}


/*************************************/
int rmultnom(unsigned long *aiCount, 
             unsigned long n, 
             double *adPr, 
             unsigned long k)
{
    unsigned long i;
    unsigned long j;
    unsigned long ulSum = 0;
    double dPick;
    double dProb;

    assert(adPr!=NULL);

    for(i=0; i<k; i++)
    {
        aiCount[i] = 0;
    }

    for(j=0; j<n; j++)
    {
        dProb = 0.0;
        dPick = runi();
        for(i=0; i<k-1; i++)
        {
            dProb += adPr[i];
            if (dPick < dProb)
            {
                aiCount[i]++;
                ulSum++;
                break;
            }
        }
    }
    aiCount[k-1] = n - ulSum;
    return 1;
}


/*************************************/
unsigned long rmultnom1
(
    double *adPr, 
    unsigned long k
)
{
    double dP = runi();
    double dSum = 0.0; 
    unsigned long i = 0;

    assert(adPr!=NULL);

    for(i=0; i<k; i++)
    {
        dSum = dSum + adPr[i];
        if(dSum > dP)
        {
            return i;
        }
    }

    return k-1;
}


/*************************************/
double rlogistic(double dLoc,double dScale)
{
    return 1.0/(1.0-runi())-1.0;
}


/*************************************/
double rlnorm(double dNorm, double dNorsd)
{
    return exp(rnorm(dNorm,dNorsd));
}


/*************************************/
int rbin(int n,double p)
{
   int j=0;
   int k=0;
   for(k=0; k<n; k++)
   {
       if(runi()<p)
       {
           j++;
       }
   }
   return j;
}


/*************************************/
double rweibull(double dGamma)
{
   if(dGamma<=0.0)
   {
       return 0.0;
   }
   return pow(rexp(1.0),(1/dGamma));
}


/*************************************/
double rchisq(double dT)
{
    return rgamma(dT/2.0)*2.0;
}


/*************************************/
double rf(double dT, double dU)
{
    if((dT<=0) || (dU <= 0))
    {
        return 0.0;
    }
    return rchisq(dT)*dU/(dT*rchisq(dU));
}


/*************************************/
double rstudent(double dT)
{
   if(dT<=0.0)
   {
       return 0.0;
   }
   return rnorm(0.0,1.0)/sqrt(rchisq(dT)/dT);
}


/*************************************/
int rperm(unsigned long *aiResult, int n)
{
    int i;
    int iSample;
    int iSpare;
    int iNew;

    assert(aiResult!=NULL);

    for(i=0; i<n; i++)
    {
        iNew = runii(i,n-1);
        iSample = aiResult[iNew];
        iSpare = aiResult[i];
        aiResult[i] = iSample;
        aiResult[iNew] = iSpare;
    }

    return 1;
}


/*************************************/
int rperm_some
(
    unsigned long *aiResult, 
    unsigned long cLength, 
    unsigned long cPermute
)
{
    unsigned long i=0;
    unsigned long ulTemp=0;
    unsigned long iNew=0;

    assert(aiResult!=NULL);

    for(i=0; i<cPermute; i++)
    {
        iNew = runii(i,cLength-1);
        ulTemp = aiResult[iNew];
        aiResult[iNew] = aiResult[i];
        aiResult[i] = ulTemp;
    }

    return 0;
}


double rgamma_rs(double dAlpha)
{
    double dMu,dSigma,dTangent,dCauchy;
    double dRatio,dVert,dHoriz;
    do 
    { // repeat until acceptance
        do 
        { // generate positive Cauchy
            do 
            { // tangent of random angle
                dVert = runi();
                dHoriz = 2.0*runi()-1.0;
            } while(dVert*dVert+dHoriz*dHoriz > 1.0);
            dTangent = dVert/dHoriz;
            dMu = dAlpha-1.0;
            dSigma = sqrt(2.0*dAlpha-1.0);
            dCauchy = dSigma*dTangent+dMu;
        } while (dCauchy <= 0.0);
        dRatio = (1.0+dTangent*dTangent)*
                 exp(dMu*log(dCauchy/dMu)-dSigma*dTangent);
    } while (runi() > dRatio);
    return dCauchy;
}


double rtruncnorm
// uses inverse cdf method
(
    double dL,
    double dR
)
{
    return qtruncnorm(runi(),dL,dR);
}





// pdf, cdf, etc
double pnorm
(
    double u
)
/*
 * The standard normal CDF, for one random variable.
 *
 *   Author:  W. J. Cody
 *   URL:   http://www.netlib.org/specfun/erf
 *
 * This is the erfc() routine only, adapted by the
 * transform stdnormal_cdf(u)=(erfc(-u/sqrt(2))/2;
 */
{
    const double a[5] = {1.161110663653770e-002,3.951404679838207e-001,
                         2.846603853776254e+001,1.887426188426510e+002,
                         3.209377589138469e+003};
    const double b[5] = {1.767766952966369e-001,8.344316438579620e+000,
                         1.725514762600375e+002,1.813893686502485e+003,
                         8.044716608901563e+003};
    const double c[9] = {2.15311535474403846e-8,5.64188496988670089e-1,
                         8.88314979438837594e+00,6.61191906371416295e+01,
                         2.98635138197400131e+02,8.81952221241769090e+02,
                         1.71204761263407058e+03,2.05107837782607147e+03,
                         1.23033935479799725e+03};
    const double d[9] = {1.00000000000000000e00,1.57449261107098347e01,
                         1.17693950891312499e02,5.37181101862009858e02,
                         1.62138957456669019e03,3.29079923573345963e03,
                         4.36261909014324716e03,3.43936767414372164e03,
                         1.23033935480374942e03};
    const double p[6] = {1.63153871373020978e-2,3.05326634961232344e-1,
                         3.60344899949804439e-1,1.25781726111229246e-1,
                         1.60837851487422766e-2,6.58749161529837803e-4};
    const double q[6] = {1.00000000000000000e00,2.56852019228982242e00,
                         1.87295284992346047e00,5.27905102951428412e-1,
                         6.05183413124413191e-2,2.33520497626869185e-3};
    const double k_sqrt2 = 1.4142135623731;
    const double k_invsqrt2 = 0.70710678118655;
    const double k_sqrtpi = 1.7724538509055;
    register double y, z;

    if(u==HUGE_VAL) return 1.0;
    if(u==-HUGE_VAL) return 0.0;

    y = fabs(u);

    if(y <= 0.46875*k_sqrt2) 
    {
        /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
        z = y*y;
        y = u*((((a[0]*z+a[1])*z+a[2])*z+a[3])*z+a[4])
            /((((b[0]*z+b[1])*z+b[2])*z+b[3])*z+b[4]);
        return 0.5+y;
    }
    z = 0.5*exp(-0.5*y*y);
    if(y <= 4.0) 
    {
        /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
        y = y*k_invsqrt2;
        y = ((((((((c[0]*y+c[1])*y+c[2])*y+c[3])*y+c[4])*y+c[5])*y+c[6])*y+c[7])*y+c[8])
            /((((((((d[0]*y+d[1])*y+d[2])*y+d[3])*y+d[4])*y+d[5])*y+d[6])*y+d[7])*y+d[8]);
        y = z*y;
    } 
    else 
    {
        /* evaluate erfc() for |u| > sqrt(2)*4.0 */
        z = z*k_sqrt2/y;
        y = 2.0/(y*y);
        y = y*(((((p[0]*y+p[1])*y+p[2])*y+p[3])*y+p[4])*y+p[5])
            /(((((q[0]*y+q[1])*y+q[2])*y+q[3])*y+q[4])*y+q[5]);
        y = z*(k_sqrtpi-y);
    }
    return (u < 0.0 ? y : 1-y);
};



double qnorm
//   Author:      Peter J. Acklam <jacklam@math.uio.no>
//   URL:         http://www.math.uio.no/~jacklam
(
    double dP
)
{
    const double a[6] = {-3.969683028665376e+01,  2.209460984245205e+02,
                         -2.759285104469687e+02,  1.383577518672690e+02,
                         -3.066479806614716e+01,  2.506628277459239e+00};
    const double b[5] = {-5.447609879822406e+01,  1.615858368580409e+02,
                         -1.556989798598866e+02,  6.680131188771972e+01,
                         -1.328068155288572e+01};
    const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                         -2.400758277161838e+00, -2.549732539343734e+00,
                         4.374664141464968e+00,  2.938163982698783e+00};
    const double d[4] = {7.784695709041462e-03,  3.224671290700398e-01,
                         2.445134137142996e+00,  3.754408661907416e+00};

    register double q, t, u;

    if (dP == 0.0)
        return -HUGE_VAL;
    if (dP == 1.0)
        return HUGE_VAL;
    q = (dP < 0.5 ? dP : 1.0-dP);
    if (q > 0.02425) 
    {
        /* Rational approximation for central region. */
        u = q-0.5;
        t = u*u;
        u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
            /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1.0);
    }
    else 
    {
        /* Rational approximation for tail region. */
        t = sqrt(-2.0*log(q));
        u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
            /((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1.0);
    }
 /* The relative error of the approximation has absolute value less
    than 1.15e-9.  One iteration of Halley's rational method (third
    order) gives full machine precision... */
    //t = pnorm(u)-q;    /* error */
    //t = 2.506628*t*exp(u*u/2);   /* f(u)/df(u) */
    //u = u-t/(1+u*t/2);     /* Halley's method */

    return (dP > 0.5 ? -u : u);
}



double qtruncnorm
(
    double dP,
    double dL,
    double dR
)
{
    static double dPLeft;

    if(dP==0.0) return dL;
    if(dP==1.0) return dR;

    dPLeft = pnorm(dL);
    return qnorm(dP*(pnorm(dR)-dPLeft) + dPLeft);
}


