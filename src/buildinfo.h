
#ifndef BUILDINFO_H
#define BUILDINFO_H

    #include <R.h>

    #ifdef IEEE_754
    #undef ISNAN
    #define ISNAN(x)   R_IsNaNorNA(x)
    #endif

    #define HH_FAILED(hr) ((unsigned long)hr != 0)
    typedef unsigned long HHRESULT;
    #define HH_OK 0
    #define HH_FAIL 1
    #define HH_INVALIDARG 2
    #define HH_OUTOFMEMORY 3
    #define HH_INVALID_DATA 4
    #define HH_NOTIMPL 5


#endif // BUILDINFO_H
