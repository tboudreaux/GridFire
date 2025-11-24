#pragma once
#include <unordered_map>
#include <nvector/nvector_serial.h>

#include "gridfire/exceptions/error_solver.h"
#include "sundials/sundials_nvector.h"

namespace gridfire::utils {
    enum class SUNDIALS_RET_CODE_TYPES {
        CVODE,
        KINSOL
    };
    static inline std::unordered_map<int, std::string> cvode_ret_code_map {
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
    static inline std::unordered_map<int, std::string> kinsol_ret_code_map {
        {0, "KIN_SUCCESS: The solver succeeded."},
        {1, "KIN_INITIAL_GUESS_OKAY: The guess, u=u0, satisfied the system F(u) = 0 within the tolerances specified"},
        {2, "KIN_STEP_LT_STPTOL: KINSOL stopped based on scaled step length. This means that the current iterate may be an approximate solution of the given nonlinear system, but it is also quite possible that the algorithm is “stalled” (making insufficient progress) near an invalid solution, or that the scalar, scsteptol, is too large."},
        {99, "KIN_WARNING: KINSOL succeeded but in an unusual way"},
        {-1, "KIN_MEM_NULL: The KINSOL memory pointer is NULL."},
        {-2, "KIN_ILL_INPUT: An illegal value was specified for an input argument."},
        {-3, "KIN_NO_MALLOC: The KINSOL memory structure has not been allocated."},
        {-4, "KIN_MEM_FAIL: A memory allocation failed."},
        {-5, "KIN_LINESEARCH_NONCONV: The line search algorithm was unable to find an iterate sufficiently distinct from the current iterate."},
        {-6, "KIN_MAXITER_REACHED: The maximum number of iterations was reached before convergence."},
        {-7, "KIN_MXNEWT_5X_EXCEEDED: Five consecutive steps have been taken that satisfy a scaled step length test."},
        {-8, "KIN_LINESEARCH_BCFAIL: The line search algorithm was unable to satisfy the beta-condition for nbcf fails iterations."},
        {-9, "KIN_LINSOLV_NO_RECOVERY: The linear solver's solve function failed recoverably, but the Jacobian data is already current."},
        {-10, "KIN_LINIT_FAIL: The linear solver's init routine failed."},
        {-11, "KIN_LSETUP_FAIL: The linear solver's setup function failed in an unrecoverable manner."},
        {-12, "KIN_LSOLVE_FAIL: The linear solver's solve function failed in an unrecoverable manner."},
        {-13, "KIN_SYSFUNC_FAIL: The system function failed in an unrecoverable manner."},
        {-14, "KIN_FIRST_SYSFUNC_ERR: The system function failed at the first call."},
        {-15, "KIN_REPTD_SYSFUNC_ERR: Unable to correct repeated recoverable system function errors."},
        {-16, "KIN_VECTOROP_ERR: A vector operation failed."},
        {-17, "KIN_CONTEXT_ERR: A context error occurred."},
        {-18, "KIN_DAMPING_FN_ERR: The damping function failed."},
        {-19, "KIN_DEPTH_FN_ERR: The depth function failed."}
    };

    inline const std::unordered_map<int, std::string>& sundials_retcode_map(const SUNDIALS_RET_CODE_TYPES type) {
        switch (type) {
            case SUNDIALS_RET_CODE_TYPES::CVODE:
                return cvode_ret_code_map;
            case SUNDIALS_RET_CODE_TYPES::KINSOL:
                return kinsol_ret_code_map;
            default:
                throw exceptions::CVODESolverFailureError("Unknown SUNDIALS return code type.");
        }

    }

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

    inline void check_sundials_flag(const int flag, const std::string& func_name, const SUNDIALS_RET_CODE_TYPES type) {
        if (flag < 0) {
            const auto& ret_code_map = sundials_retcode_map(type);
            if (!ret_code_map.contains(flag)) {
                switch (type) {
                    case (SUNDIALS_RET_CODE_TYPES::CVODE):
                        throw exceptions::CVODESolverFailureError("CVODE error in " + func_name + ": Unknown error code: " + std::to_string(flag));
                    case (SUNDIALS_RET_CODE_TYPES::KINSOL):
                        throw exceptions::KINSolSolverFailureError("KINSol error in " + func_name + ": Unknown error code: " + std::to_string(flag));
                    default:
                        throw exceptions::CVODESolverFailureError("SUNDIALS error in " + func_name + ": Unknown error code: " + std::to_string(flag));
                }
            }
            switch (type) {
                case (SUNDIALS_RET_CODE_TYPES::CVODE):
                    throw exceptions::CVODESolverFailureError("CVODE error in " + func_name + ": " + ret_code_map.at(flag));
                case (SUNDIALS_RET_CODE_TYPES::KINSOL):
                    throw exceptions::KINSolSolverFailureError("KINSol error in " + func_name + ": " + ret_code_map.at(flag));
                default:
                    throw exceptions::CVODESolverFailureError("SUNDIALS error in " + func_name + ": " + ret_code_map.at(flag));
            }
        }
    }
}
