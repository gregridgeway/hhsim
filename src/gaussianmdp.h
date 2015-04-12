
#ifndef GAUSSIANMDP_H
#define GAUSSIANMDP_H

#include "hhobject.h"
#include <Rmath.h>

class CGaussianMDP : public CHHObject
{
public:

   CGaussianMDP();

   virtual ~CGaussianMDP();

   HHRESULT initialize(double *Y,
                       double *sigma2,
                       unsigned long cY,
                       int K,
                       int *config,
                       double *phi,
                       double *adPriorParams);  // m, s2, alpha, w, W
   HHRESULT sample();
   HHRESULT getTheta(double *adTheta);

   double m;
   double s2;
   double alpha;
   int numconfig;

private:
   HHRESULT sample_config(int *&config,
                          int obs,
                          double *sigma2,
                          int n,
                          double *y,
                          double *phi,
                          double alpha);
   double sample_phi0(int n,
                      double *ysamp,
                      double *sigmasamp,
                      double s2,
                      double m);
   HHRESULT sample_phi(int *config,
                       double *y,
                       double *sigma2,
                       double s2,
                       double mn,
                       int n,
                       double *&phi,
                       int K);
   HHRESULT sample_theta(int *config,
                         double *phi,
                         int n,
                         double *&theta);
   HHRESULT sample_m(double s2,
                     double *phi,
                     int sizephi,
                     double &m);
   HHRESULT sample_s2(double &s2,
                      double *phi,
                      double m,
                      int sizephi,
                      double w,
                      double W);
   HHRESULT sample_alpha(double par1,
                         double par2,
                         int n,
                         int k,
                         double &alpha);
   int multinomial(int ncell, 
                   double * nvec);

   double par1;
   double par2;
   double w;
   double W;
   
   double *theta;

   int K;
   int *config;
   double *phi;
   
   double *Y;            // data
   double *sigma2;
   
   int *configtmp;
   int *nconfig;
   double *prob;
   double *ysamp;
   double *sigmasamp;

   unsigned long i;
};

#endif // GAUSSIANMDP_H
