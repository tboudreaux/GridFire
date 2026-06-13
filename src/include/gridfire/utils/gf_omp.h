#pragma once
#include "fourdst/logging/logging.h"
#include "quill/LogMacros.h"

#if defined(GF_USE_OPENMP)

#include <omp.h>

namespace gridfire::omp {
    /**
     * @brief Namespace containing utilities for initializing and managing parallel execution with OpenMP and CppAD.
     *
     * @note GF_PAR_INIT should be called at the start of your program regardless of if you run in parallel or serial
     * mode. When GridFire has been compiled without openMP support, GF_PAR_INIT will simply log a message that you are not in parallel mode and return.
     * However, if in the future you wish to relink against a version which has been compiled with parallel support
     * missing GF_PAR_INIT may lead to a silent failure.
     *
     * @note An end user should only ever need to call the GF_PAR_INIT macro. i.e. never call any of the actual
     * functions in this header directly.
     */
    inline bool s_par_mode_initialized = false;

    inline unsigned long get_thread_id() {
        return static_cast<unsigned long>(omp_get_thread_num());
    }

    inline bool in_parallel() {
        return omp_in_parallel() != 0;
    }

    inline void init_parallel_mode() {
        if (s_par_mode_initialized) return;
        [[maybe_unused]] quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");

        auto n_threads = static_cast<size_t>(omp_get_max_threads());
        if (n_threads > CPPAD_MAX_NUM_THREADS) {
            LOG_CRITICAL(logger,
                "OpenMP reports {} threads but CppAD was built with CPPAD_MAX_NUM_THREADS={}; "
                "clamping OpenMP to {} threads.",
                n_threads, static_cast<size_t>(CPPAD_MAX_NUM_THREADS),
                static_cast<size_t>(CPPAD_MAX_NUM_THREADS));
            throw std::runtime_error(std::format(
                "OpenMP reports {} threads but CppAD was built with CPPAD_MAX_NUM_THREADS={}; clamping default OpenMP number of threads to {}. Rebuild CppAD with a higher CPPAD_MAX_NUM_THREADS if you need more threads. Alternative, set the environmental variable OMP_NUM_THREADS to a value less than or equal to {} to avoid this error.",
                n_threads, static_cast<size_t>(CPPAD_MAX_NUM_THREADS),
                static_cast<size_t>(CPPAD_MAX_NUM_THREADS),
                static_cast<size_t>(CPPAD_MAX_NUM_THREADS)));
        }

        LOG_INFO(logger, "Initializing OpenMP parallel mode with {} threads", n_threads);
        CppAD::thread_alloc::parallel_setup(
            n_threads,
            []() -> bool { return in_parallel(); },
            []() -> size_t { return get_thread_id(); });
        CppAD::thread_alloc::hold_memory(true);
        CppAD::parallel_ad<double>();
        s_par_mode_initialized = true;
    }
}

#define GF_PAR_INIT() gridfire::omp::init_parallel_mode();

#else

namespace gridfire::omp {
    inline void log_not_in_parallel_mode() {
        [[maybe_unused]] quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        LOG_INFO(logger, "This is not an error! Note: OpenMP parallel mode is not enabled. GF_USE_OPENMP is not defined. Pass -DGF_USE_OPENMP when compiling to enable OpenMP support. When using meson use the option -Dopenmp_support=true");
    }
}

#define GF_PAR_INIT() gridfire::omp::log_not_in_parallel_mode();

#endif