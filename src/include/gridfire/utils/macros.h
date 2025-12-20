#pragma once


#if defined(GF_USE_OPENMP)
    #define GF_OMP_PRAGMA(x) _Pragma(#x)
    #define GF_OMP(omp_args, extra) GF_OMP_PRAGMA(omp omp_args) extra
    #define GF_OMP_MAX_THREADS omp_get_max_threads()
    #define GF_OMP_THREAD_NUM omp_get_thread_num()
#else
    #define GF_OMP(_,fallback_args) fallback_args
    #define GF_OMP_MAX_THREADS 1
    #define GF_OMP_THREAD_NUM 0
#endif