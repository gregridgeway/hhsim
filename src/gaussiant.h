
#ifndef GAUSSIANT_H
#define GAUSSIANT_H

#include "hhobject.h"

class CGaussianT : public CHHObject
{
public:

   CGaussianT();

   virtual ~CGaussianT();

   HHRESULT initialize(double *adY,            // data
		   	double *adSigma2,
                       unsigned long cY,
                       double *adPriorParams);  // prior parameters s1, s2, a, A, t1, t2

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
   DiagonalMatrix dmatTau2Lambdainv;
   ColumnVector vecOne;
   ColumnVector vecd;
   ColumnVector vecMu;
   LowerTriangularMatrix ltmatSqrtD;
   double dTemp1;
   double dTemp2;

   ColumnVector vecLambda;
   DiagonalMatrix dmatIdentity;
   DiagonalMatrix dmatLambda;
   DiagonalMatrix dmatLambdainv;
   DiagonalMatrix datTau2Lambdainv; 

   double ds1;
   double ds2;
   double da;
   double dA;
   double dt1;
   double dt2;
   double dnu;

   unsigned long i;
   unsigned long j;
};

#endif // GAUSSIANT_H
