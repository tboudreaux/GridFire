#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/engine/scratchpads/scratchpad_abstract.h"
#include "gridfire/engine/scratchpads/types.h"
#include "gridfire/engine/scratchpads/engine_multiscale_scratchpad.h"

#include "fourdst/atomic/atomicSpecies.h"

#include <vector>
#include <memory>
#include <unordered_map>

#include "sundials/sundials_context.h"


namespace gridfire::engine::scratch {
     bool MultiscalePartitioningEngineViewScratchPad::is_initialized() const{ return has_initialized; }

     void MultiscalePartitioningEngineViewScratchPad::initialize() {
          if (has_initialized) return;

          const int flag = SUNContext_Create(SUN_COMM_NULL, &sun_ctx);
          if (flag != 0) {
               throw std::runtime_error("Failed to create SUNContext in MultiscalePartitioningEngineViewScratchPad.");
          }

          has_initialized = true;
     }

     MultiscalePartitioningEngineViewScratchPad::~MultiscalePartitioningEngineViewScratchPad() {
          qse_solvers.clear();
          if (sun_ctx != nullptr) {
               SUNContext_Free(&sun_ctx);
               sun_ctx = nullptr;
          }
     }

     std::unique_ptr<AbstractScratchPad> MultiscalePartitioningEngineViewScratchPad::clone() const {
          auto clone_pad = std::make_unique<MultiscalePartitioningEngineViewScratchPad>();
          clone_pad->qse_groups = this->qse_groups;
          clone_pad->dynamic_species = this->dynamic_species;
          clone_pad->algebraic_species = this->algebraic_species;
          clone_pad->composition_cache = this->composition_cache;

          clone_pad->initialize();
          clone_pad->qse_solvers.reserve(this->qse_solvers.size());
          for (const auto& solver : qse_solvers) {
               clone_pad->qse_solvers.push_back(solver->clone(clone_pad->sun_ctx)); // Must rebind context to new SUNContext
          }
          return clone_pad;
     }

}