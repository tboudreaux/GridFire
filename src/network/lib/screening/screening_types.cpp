#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/screening/screening_weak.h"
#include "gridfire/screening/screening_bare.h"

#include <memory>

namespace gridfire::screening {
    std::unique_ptr<ScreeningModel> selectScreeningModel(const ScreeningType type) {
        switch (type) {
            case ScreeningType::WEAK:
                return std::make_unique<WeakScreeningModel>();
            case ScreeningType::BARE:
                return std::make_unique<BareScreeningModel>();
            default:
                return std::make_unique<BareScreeningModel>();
        }
    }
}
