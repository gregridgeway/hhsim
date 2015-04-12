
#include "gaussiangaussian.h"

CGaussianGaussian::CGaussianGaussian()
{
}

CGaussianGaussian::~CGaussianGaussian()
{
}

HHRESULT CGaussianGaussian::initialize
(
    double *adY,            // data
    double *adSigma2,       // fixed sigma2 params
    unsigned long cY,
    double *adPriorParams  // prior parameters s1, s2, a, A, t1, t2
)
{
   HHRESULT hr = HH_OK;
   fInitialized = true;

   if(adY==NULL)
   {
      hr = HH_INVALIDARG;
      goto Error;
   }
   if(adSigma2==NULL)
   {
      hr = HH_INVALIDARG;
      goto Error;
   }
   if(adPriorParams==NULL)
   {
      hr = HH_INVALIDARG;
      goto Error;
   }

   ds1 = adPriorParams[0];
   ds2 = adPriorParams[1];
   da  = adPriorParams[2];
   dA  = adPriorParams[3];
   dt1 = adPriorParams[4];
   dt2 = adPriorParams[5];

   this->cY = cY;

   vecY.ReSize(cY);
   vecSigma2.ReSize(cY);
   vecTheta.ReSize(cY);
   dmatSigma.ReSize(cY);

   // must initialize by drawing these from their prior distributions
   dm = rnorm(da,pow(dA,0.5));
   dTau2 = 1/rgamma2(dt1,dt2);

   symmatD.ReSize(cY);
   dmatTau2.ReSize(cY);
   vecOne.ReSize(cY);
   vecd.ReSize(cY);
   vecMu.ReSize(cY);
   ltmatSqrtD.ReSize(cY);

   dTemp1 = 0.0;
   dTemp2 = 0.0;

   vecY << adY;
   vecSigma2 << adSigma2;
   vecTheta = 0.0;
   for(i=1;i<=cY;i++)
   {
      dmatSigma(i,i) = vecSigma2(i);
   }
   vecOne = 1.0;

Cleanup:
    return hr;
Error:
    goto Cleanup;
} // end test




HHRESULT CGaussianGaussian::sample()
{
   HHRESULT hr = HH_OK;

   dmatTau2 = dTau2;

   // draw new Theta
   dTemp1 = 1./((vecOne.t()*dmatTau2.i()*vecOne + 1/dA).AsScalar());
   symmatD << (dmatSigma.i() + dmatTau2.i() - dmatTau2.i()*
               vecOne*dTemp1*vecOne.t()*dmatTau2.i()).i();

   vecd = dmatSigma.i()*vecY + da*(dmatTau2 + dA*vecOne*vecOne.t()).i()*vecOne;
   vecMu = symmatD*vecd;
   ltmatSqrtD = Cholesky(symmatD);
   hr = rmvnorm(vecMu, ltmatSqrtD, vecTheta);

   // draw new Sigma
   // 2-2-04: fix for initial simulation studies
   //   for(i=1; i<=cY; i++)
   //   {
   //      dTemp1 = vecY(i)-vecTheta(i);
   //      dTemp2 = rgamma2(ds1+0.5, ds2+0.5*dTemp1*dTemp1);
   //      dmatSigma(i,i) = 1.0/dTemp2;
   //   }

   // draw new m (mean of thetas)
   dTemp1 = 1./((vecOne.t()*dmatTau2.i()*vecOne + 1/dA).AsScalar());
   dTemp2 = (vecOne.t()*dmatTau2.i()*vecTheta).AsScalar() + da/dA;
   dm =  rnorm(dTemp1*dTemp2, pow(dTemp1,0.5));

   // draw new tau2
   dTemp1 = rgamma2(dt1 + 0.5*cY,
                   dt2 + 0.5*(vecTheta-dm).SumSquare());
   dTau2 = 1.0/dTemp1;

Cleanup:
   return hr;
Error:
    goto Cleanup;
}



HHRESULT CGaussianGaussian::rmvnorm
(
    ColumnVector &cvMu,
    LowerTriangularMatrix &matSqrtSigma,
    ColumnVector &cvX
)
{
   HHRESULT hr = HH_OK;

   unsigned long i=0;
   for(i=0; i<cvMu.Nrows(); i++)
   {
      cvX[i] = rnorm(0.0,1.0); // NOTE: this takes the sd, not the variance
   }
   cvX = matSqrtSigma*cvX + cvMu;

   return hr;
}


HHRESULT CGaussianGaussian::getSigma
(
   double *adSigma
)
{
   HHRESULT hr = HH_OK;

   unsigned long i = 0;
   for(i=0; i<cY; i++)
   {
      adSigma[i] = dmatSigma(i+1);
   }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}
