#pragma once


#if defined(GF_USE_OPENMP)
    #include <omp.h>
    #include <algorithm>
    #include "cppad/configure.hpp"

    #define GF_OMP_PRAGMA(x) _Pragma(#x)

    namespace gridfire::omp {
        inline int capped_max_threads() {
            return std::min<int>(omp_get_max_threads(),
                                  static_cast<int>(CPPAD_MAX_NUM_THREADS));
        }
    }

    #define GF_OMP_NUM_THREADS (gridfire::omp::capped_max_threads())
    #define GF_OMP(omp_args, extra) \
    GF_OMP_PRAGMA(omp omp_args num_threads(GF_OMP_NUM_THREADS)) extra

    #define GF_OMP_MAX_THREADS (gridfire::omp::capped_max_threads())
    #define GF_OMP_THREAD_NUM  omp_get_thread_num()
#else
    #define GF_OMP(_, fallback_args) fallback_args
    #define GF_OMP_MAX_THREADS 1
    #define GF_OMP_THREAD_NUM 0
#endif