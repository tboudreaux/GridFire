#include "gridfire/utils/gf_omp.h"

#if defined(GF_USE_OPENMP)

#include "cppad/cppad.hpp"

#include <omp.h>

#include <cstdio>
#include <cstdlib>
#include <format>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace gridfire::omp {
    namespace {
        constexpr size_t UNASSIGNED = std::numeric_limits<size_t>::max();

        // Never destroyed: runtime threads (OpenMP pools) release their numbers during process exit, after static
        // objects would have been destroyed
        std::mutex& ids_mutex() {
            static auto* m = new std::mutex;
            return *m;
        }
        std::vector<size_t>& free_ids() {  // numbers released by threads that exited (guarded by ids_mutex)
            static auto* ids = new std::vector<size_t>;
            return *ids;
        }
        size_t next_id = 1;            // next never-used number; 0 belongs to the initialising thread
        bool initialized = false;

        size_t acquire_id() {
            std::lock_guard lock(ids_mutex());
            if (!free_ids().empty()) {
                const size_t id = free_ids().back();
                free_ids().pop_back();
                return id;
            }
            if (next_id >= CPPAD_MAX_NUM_THREADS) {
                // Called from inside CppAD, so report and stop rather than throw through its frames
                std::fprintf(stderr, "GridFire: more than %zu threads are using GridFire at once (CPPAD_MAX_NUM_THREADS); "
                             "rebuild CppAD with a larger CPPAD_MAX_NUM_THREADS.\n", static_cast<size_t>(CPPAD_MAX_NUM_THREADS));
                std::abort();
            }
            return next_id++;
        }

        /// The calling thread's number, returned to the pool when the thread exits
        struct ThreadSlot {
            size_t id = UNASSIGNED;
            ~ThreadSlot() {
                if (id != UNASSIGNED && id != 0) {
                    std::lock_guard lock(ids_mutex());
                    free_ids().push_back(id);
                }
            }
        };
        thread_local ThreadSlot slot;

        bool cppad_in_parallel() { return in_parallel(); }
        size_t cppad_thread_num() { return thread_id(); }
    }

    size_t thread_id() {
        if (slot.id == UNASSIGNED) {
            slot.id = acquire_id();
        }
        return slot.id;
    }

    bool in_parallel() {
        return omp_in_parallel() != 0 || thread_id() != 0;
    }

    bool parallel_mode_initialized() {
        return initialized;
    }

    void init_parallel_mode() {
        if (initialized) return;
        [[maybe_unused]] quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        if (omp_in_parallel()) {
            throw std::logic_error("gridfire::omp::init_parallel_mode must be called outside any parallel region.");
        }
        if (static_cast<size_t>(omp_get_max_threads()) > CPPAD_MAX_NUM_THREADS) {
            throw std::runtime_error(std::format(
                "OpenMP reports {} threads but CppAD was built with CPPAD_MAX_NUM_THREADS={}. Set OMP_NUM_THREADS to at "
                "most {}, or rebuild CppAD with a larger CPPAD_MAX_NUM_THREADS.",
                omp_get_max_threads(), static_cast<size_t>(CPPAD_MAX_NUM_THREADS), static_cast<size_t>(CPPAD_MAX_NUM_THREADS)));
        }

        slot.id = 0;
        LOG_INFO(logger, "Initializing CppAD for up to {} threads", static_cast<size_t>(CPPAD_MAX_NUM_THREADS));
        CppAD::thread_alloc::parallel_setup(CPPAD_MAX_NUM_THREADS, cppad_in_parallel, cppad_thread_num);
        CppAD::thread_alloc::hold_memory(true);
        CppAD::parallel_ad<double>();
        initialized = true;
    }
}

#endif
