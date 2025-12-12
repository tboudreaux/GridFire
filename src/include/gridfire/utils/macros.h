#pragma once


#if defined(GF_USE_OPENMP)
    #define GF_OMP_PRAGMA(x) _Pragma(#x)
    #define GF_OMP(omp_args, _) GF_OMP_PRAGMA(omp omp_args)
#else
    #define GF_OMP(_,fallback_args) fallback_args
#endif