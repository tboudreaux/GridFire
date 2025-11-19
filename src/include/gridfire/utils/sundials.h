#pragma once
#include <unordered_map>
#include <nvector/nvector_serial.h>

#include "gridfire/exceptions/error_solver.h"
#include "sundials/sundials_nvector.h"

namespace gridfire::utils {
    inline std::unordered_map<int, std::string> cvode_ret_code_map {
        {0, "CV_SUCCESS: The solver succeeded."},
        {1, "CV_TSTOP_RETURN: The solver reached the specified stopping time."},
        {2, "CV_ROOT_RETURN: A root was found."},
        {-99, "CV_WARNING: CVODE succeeded but in an unusual manner"},
        {-1, "CV_TOO_MUCH_WORK: The solver took too many internal steps."},
        {-2, "CV_TOO_MUCH_ACC: The solver could not satisfy the accuracy requested."},
        {-3, "CV_ERR_FAILURE: The solver encountered a non-recoverable error."},
        {-4, "CV_CONV_FAILURE: The solver failed to converge."},
        {-5, "CV_LINIT_FAIL: The linear solver's initialization function failed."},
        {-6, "CV_LSETUP_FAIL: The linear solver's setup function failed."},
        {-7, "CV_LSOLVE_FAIL: The linear solver's solve function failed."},
        {-8, "CV_RHSFUNC_FAIL: The right-hand side function failed in an unrecoverable manner."},
        {-9, "CV_FIRST_RHSFUNC_ERR: The right-hand side function failed at the first call."},
        {-10, "CV_REPTD_RHSFUNC_ERR: The right-hand side function repeatedly failed recoverable."},
        {-11, "CV_UNREC_RHSFUNC_ERR: The right-hand side function failed unrecoverably."},
        {-12, "CV_RTFUNC_FAIL: The rootfinding function failed in an unrecoverable manner."},
        {-13, "CV_NLS_INIT_FAIL: The nonlinear solver's initialization function failed."},
        {-14, "CV_NLS_SETUP_FAIL: The nonlinear solver's setup function failed."},
        {-15, "CV_CONSTR_FAIL : The inequality constraint was violated and the solver was unable to recover."},
        {-16, "CV_NLS_FAIL: The nonlinear solver's solve function failed."},
        {-20, "CV_MEM_FAIL: Memory allocation failed."},
        {-21, "CV_MEM_NULL: The CVODE memory structure is NULL."},
        {-22, "CV_ILL_INPUT: An illegal input was detected."},
        {-23, "CV_NO_MALLOC: The CVODE memory structure has not been allocated."},
        {-24, "CV_BAD_K: The value of k is invalid."},
        {-25, "CV_BAD_T: The value of t is invalid."},
        {-26, "CV_BAD_DKY: The value of dky is invalid."},
        {-27, "CV_TOO_CLOSE: The time points are too close together."},
        {-28, "CV_VECTOROP_ERR: A vector operation failed."},
        {-29, "CV_PROJ_MEM_NULL: The projection memory structure is NULL."},
        {-30, "CV_PROJFUNC_FAIL: The projection function failed in an unrecoverable manner."},
        {-31, "CV_REPTD_PROJFUNC_ERR: The projection function has repeated recoverable errors."}
    };

    inline void check_cvode_flag(const int flag, const std::string& func_name) {
        if (flag < 0) {
            if (!cvode_ret_code_map.contains(flag)) {
                throw exceptions::CVODESolverFailureError("CVODE error in " + func_name + ": Unknown error code: " + std::to_string(flag));
            }
            throw exceptions::CVODESolverFailureError("CVODE error in " + func_name + ": " + cvode_ret_code_map.at(flag));
        }
    }

    inline N_Vector init_sun_vector(uint64_t size, SUNContext sun_ctx) {
        #ifdef SUNDIALS_HAVE_OPENMP
            const N_Vector vec = N_VNew_OpenMP(size, 0, sun_ctx);
        #elif SUNDIALS_HAVE_PTHREADS
            const N_Vector vec = N_VNew_Pthreads(size, sun_ctx);
        #else
            const N_Vector vec = N_VNew_Serial(static_cast<long long>(size), sun_ctx);
        #endif

        check_cvode_flag(vec == nullptr ? -1 : 0, "N_VNew");
        return vec;
    }

}
