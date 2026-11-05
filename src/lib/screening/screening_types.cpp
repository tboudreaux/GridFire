#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/screening/screening_weak.h"
#include "gridfire/screening/screening_bare.h"
#include "gridfire/screening/screening_CD09.h"
#include "gridfire/screening/screening_blended.h"

#include <memory>
#include <stdexcept>

namespace gridfire::screening {
    std::unique_ptr<ScreeningModel> selectScreeningModel(const ScreeningType type, const ScreeningBlendWindow& window) {
      switch (type) {
          case ScreeningType::WEAK:
              return std::make_unique<WeakScreeningModel>();
          case ScreeningType::BARE:
              return std::make_unique<BareScreeningModel>();
          case ScreeningType::CD09:
              return std::make_unique<CD09ScreeningModel>();
          case ScreeningType::BLENDED:
              return std::make_unique<BlendedScreeningModel>(
                  std::make_unique<WeakScreeningModel>(), std::make_unique<CD09ScreeningModel>(), window);
      }
      throw std::runtime_error("Unknown ScreeningType");
    }
}
