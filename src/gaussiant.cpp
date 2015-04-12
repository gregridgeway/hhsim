
#include "gaussiant.h"

CGaussianT::CGaussianT()
{
}

CGaussianT::~CGaussianT()
{
}

HHRESULT CGaussianT::initialize
(
    double *adY,            // data
    double *adSigma2,
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
   dnu = adPriorParams[6];

   this->cY = cY;

   vecY.ReSize(cY);
   vecSigma2.ReSize(cY);
   vecTheta.ReSize(cY);
   vecY << adY;
   vecTheta = 0.0;
   vecSigma2 << adSigma2;

   vecLambda.ReSize(cY); 
   vecLambda = 1.0; //maybe later will want to read in initial values
   dmatIdentity.ReSize(cY); 
   dmatIdentity=1.0;
   dmatLambda.ReSize(cY);
//   dmatLambda << dmatIdentity*vecLambda;  comment for now
   dmatLambda = 1.0; //for now - to test - it works

//   for(i=1;i<=cY;i++) {
//          Rprintf("%f ",dmatLambda(i,i)); // these initial vals look wrong if Lambda is not set to 1 initially
//    }

   dmatSigma.ReSize(cY);
   dm = 0.0;
   dTau2 = 1.0*(dnu-2)/dnu; // for simul study, want var of deviate to be 1

   symmatD.ReSize(cY);
   dmatTau2.ReSize(cY);
   vecOne.ReSize(cY);
   vecd.ReSize(cY);
   vecMu.ReSize(cY);
   ltmatSqrtD.ReSize(cY);

   dTemp1 = 0.0;
   dTemp2 = 0.0;
   vecOne = 1.0;

// must initialize by drawing these from their prior distributions
   dm = rnorm(da,pow(dA,0.5));
   dTau2 = 1/rgamma2(dt1,dt2);

   for(i=1;i<=cY;i++)
   {
      dmatSigma(i,i) = vecSigma2(i);
   }

Cleanup:
    return hr;
Error:
    goto Cleanup;
} // end initialize




HHRESULT CGaussianT::sample()
{
   HHRESULT hr = HH_OK;

   dmatTau2 = dTau2;

   dmatLambdainv = dmatLambda.i(); //smp: to create a mtx from Lambda vector
   dmatTau2Lambdainv = dTau2*dmatLambdainv;

   // draw new Theta
   dTemp1 = 1./((vecOne.t()*dmatTau2Lambdainv.i()*vecOne + 1/dA).AsScalar());
   symmatD << (dmatSigma.i() + dmatTau2Lambdainv.i() -
               dmatTau2Lambdainv.i()*vecOne*dTemp1*vecOne.t()*dmatTau2Lambdainv.i()).i();
   vecd = dmatSigma.i()*vecY + da*(dmatTau2Lambdainv + dA*vecOne*vecOne.t()).i()*vecOne;

   vecMu = symmatD*vecd;
   ltmatSqrtD = Cholesky(symmatD);
   hr = rmvnorm(vecMu, ltmatSqrtD, vecTheta);

   // draw new Sigma (same as for the t-distribution)
   // this fc applies to drawing a separate sigma2 for each obs
   // need to augment for the case of assumed equal sigma2
   // not critical now (1-24-04) because our simul studies will fix sigma2
   // anyway
//   for(i=1; i<=cY; i++)
//   {
//      dTemp1 = vecY(i)-vecTheta(i);
//      dTemp2 = rgamma2(ds1+0.5, ds2+0.5*dTemp1*dTemp1);
//      dmatSigma(i,i) = 1.0/dTemp2;
//   }

   // draw new m (mean of thetas)
   dTemp1 = 1./((vecOne.t()*dmatTau2Lambdainv.i()*vecOne + 1/dA).AsScalar());
   dTemp2 = (vecOne.t()*dmatTau2Lambdainv.i()*vecTheta).AsScalar() + da/dA;
   dm = rnorm(dTemp1*dTemp2, pow(dTemp1,0.5));

   // draw new tau2
//   dTemp1 = 0.0;
   for(i=1;i<=cY; i++)
   {
      dTemp2 = vecTheta(i)-dm;
      dTemp1 += dmatLambda(i,i)*dTemp2*dTemp2;
   }
   dTau2 = 1/rgamma2(dt1 + 0.5*cY, 0.5*dTemp1 + dt2);

   // draw new lambda parameters
   for(i=1;i<=cY;i++)
   {
      dTemp1 = vecTheta(i) - dm;
      dTemp2 = 0.5*(dnu + dTemp1*dTemp1/dTau2);
      dmatLambda(i,i) = rgamma2((dnu+1)/2,dTemp2);
   }

Cleanup:
   return hr;
Error:
    goto Cleanup;
}


HHRESULT CGaussianT::rmvnorm
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
      cvX[i] = rnorm(0.0,1.0);
   }
   cvX = matSqrtSigma*cvX + cvMu;

   return hr;
}


HHRESULT CGaussianT::getSigma
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
