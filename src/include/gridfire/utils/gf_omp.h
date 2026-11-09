#pragma once
#include "fourdst/logging/logging.h"
#include "quill/LogMacros.h"

#include <cstddef>

#if defined(GF_USE_OPENMP)

namespace gridfire::omp {
    /**
     * @brief Thread handling for CppAD when GridFire is built with OpenMP.
     *
     * CppAD keeps per-thread memory and tapes and needs a number for the calling thread. GridFire gives every thread
     * its own number the first time it reaches CppAD, instead of asking an OpenMP runtime: a host program may run
     * GridFire from threads created by a different OpenMP runtime (GCC's libgomp, while GridFire carries LLVM's libomp)
     * or from plain threads, and those all look like thread 0 to GridFire's runtime. Numbers are reused when a thread
     * exits. The state lives in libgridfire (not in this header) so every library and program sees the same instance.
     *
     * @note An end user only calls GF_PAR_INIT (outside any parallel region, before using GridFire from several
     * threads). The C API does this in gf_init.
     */

    /// Set up CppAD for use from several threads; the calling thread becomes thread 0. Must be called outside any
    /// parallel region; later calls do nothing.
    void init_parallel_mode();

    /// Whether init_parallel_mode has run.
    [[nodiscard]] bool parallel_mode_initialized();

    /// This thread's CppAD number (assigned on first use; 0 for the thread that called init_parallel_mode).
    [[nodiscard]] size_t thread_id();

    /// Whether CppAD should treat the calling thread as running in parallel: inside an OpenMP region of GridFire's
    /// runtime, or any thread other than thread 0.
    [[nodiscard]] bool in_parallel();
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
