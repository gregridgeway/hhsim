
#ifndef HHOBJECT_H
#define HHOBJECT_H

#include "buildinfo.h"
#include "newmat.h"
#include "newmatap.h"   // need matrix applications for Cholesky
#include "random.h"

class CHHObject
{

public:

   CHHObject();
   virtual ~CHHObject();

   virtual HHRESULT sample() = 0;
   HHRESULT getTheta(double *adTheta);

   bool fInitialized;
   ColumnVector vecTheta;
   unsigned long cY;
};
typedef CHHObject *PCHHObject;

#endif // HHOBJECT_H



