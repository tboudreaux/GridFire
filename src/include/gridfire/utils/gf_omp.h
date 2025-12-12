#pragma once
#include "fourdst/logging/logging.h"
#include "quill/LogMacros.h"

#if defined(GF_USE_OPENMP)

#include <omp.h>

namespace gridfire::omp {
    static bool s_par_mode_initialized = false;

    inline unsigned long get_thread_id() {
        return static_cast<unsigned long>(omp_get_thread_num());
    }

    inline bool in_parallel() {
        return omp_in_parallel() != 0;
    }

    inline void init_parallel_mode() {
        if (s_par_mode_initialized) {
            return; // Only initialize once
        }
        quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        LOG_INFO(logger, "Initializing OpenMP parallel mode with {} threads", static_cast<unsigned long>(omp_get_max_threads()));
        CppAD::thread_alloc::parallel_setup(
            static_cast<size_t>(omp_get_max_threads()), // Max threads
            []() -> bool { return in_parallel(); },                              // Function to get thread ID
            []() -> size_t { return get_thread_id(); }                                 // Function to check parallel state
        );

        CppAD::thread_alloc::hold_memory(true);
        s_par_mode_initialized = true;
    }
}

#define GF_PAR_INIT() gridfire::omp::init_parallel_mode();

#else

namespace gridfire::omp {
    inline void log_not_in_parallel_mode() {
        quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        LOG_INFO(logger, "This is not an error! Note: OpenMP parallel mode is not enabled. GF_USE_OPENMP is not defined. Pass -DGF_USE_OPENMP when compiling to enable OpenMP support. When using meson use the option -Dopenmp_support=true");
    }
}

#define GF_PAR_INIT() gridfire::omp::log_not_in_parallel_mode();

#endif