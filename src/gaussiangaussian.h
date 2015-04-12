
#ifndef GAUSSIANGAUSSIAN_H
#define GAUSSIANGAUSSIAN_H

#include "hhobject.h"

class CGaussianGaussian : public CHHObject
{
public:

   CGaussianGaussian();

   virtual ~CGaussianGaussian();

   HHRESULT initialize(double *adY,            // data
                       double *adSigma2, //  known sigma2 values  
                       unsigned long cY,
                       double *adPriorParams); 

   HHRESULT sample();
   HHRESULT rmvnorm(ColumnVector &cvMu,
                    LowerTriangularMatrix &matSqrtSigma,
                    ColumnVector &cvX);
   HHRESULT getSigma(double *adSigma);

   ColumnVector vecY;
   ColumnVector vecSigma2;

   DiagonalMatrix dmatSigma;
   double dm;
   double dTau2;

   SymmetricMatrix symmatD;
   DiagonalMatrix dmatTau2;
   ColumnVector vecOne;
   ColumnVector vecd;
   ColumnVector vecMu;
   LowerTriangularMatrix ltmatSqrtD;
   double dTemp1;
   double dTemp2;

private:
   double ds1;
   double ds2;
   double da;
   double dA;
   double dt1;
   double dt2;

   unsigned long i;
};

#endif // GAUSSIANGAUSSIAN_H
