#include "hhobject.h"

CHHObject::CHHObject()
{
   fInitialized = false;
}

CHHObject::~CHHObject()
{
}


HHRESULT CHHObject::getTheta
(
   double *adTheta
)
{
   HHRESULT hr = HH_OK;

   unsigned long i = 0;
   for(i=0; i<cY; i++)
   {
      adTheta[i] = vecTheta(i+1);
   }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}
