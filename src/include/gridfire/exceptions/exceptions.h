#pragma once

#include "gridfire/exceptions/error_gridfire.h"
#include "gridfire/exceptions/error_engine.h"
#include "gridfire/exceptions/error_utils.h"
#include "gridfire/exceptions/error_debug.h"
#include "gridfire/exceptions/error_policy.h"
#include "gridfire/exceptions/error_reaction.h"
#include "gridfire/exceptions/error_solver.h"
#include "gridfire/exceptions/error_utils.h"

/**
 * @namespace gridfire::exceptions
 * @brief Namespace for GridFire exception classes. All custom exceptions
 * defined in GridFire inherit from std::exception and are organized
 * within this namespace for clarity and modularity.
 *
 * All GridFire exception classes derive from std::exception and then from gridfire::exceptions::GridFireError;
 */
namespace gridfire::exceptions{}