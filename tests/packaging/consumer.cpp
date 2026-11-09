// A downstream C++ program built only from what an installed GridFire reports: builds the main-sequence network and
// evaluates the solar-core energy generation rate through the C++ API.
#include "gridfire/gridfire.h"

#include "fourdst/composition/composition.h"

#include <cstdio>
#include <string>
#include <vector>

int main() {
    const std::vector<std::string> names = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
    const std::vector<double> Y = {0.7026166, 9.7479e-06, 0.0689551, 0.00025, 7.8554e-05, 0.00060144, 8.103e-05, 2.1513e-05};
    fourdst::composition::Composition comp;
    for (size_t i = 0; i < names.size(); ++i) {
        comp.registerSymbol(names[i]);
        comp.setMolarAbundance(names[i], Y[i]);
    }
    gridfire::policy::MainSequencePolicy policy(comp);
    auto [engine, ctx] = policy.construct();
    const auto rhs = engine.calculateRHSAndEnergy(*ctx, comp, 0.015, 150.0, false);
    if (!rhs || !(rhs->nuclearEnergyGenerationRate > 1.0 && rhs->nuclearEnergyGenerationRate < 1e3)) {
        std::printf("FAIL: no plausible solar-core eps_nuc\n");
        return 1;
    }
    std::printf("C++ consumer OK: %zu reactions, eps_nuc = %.6e erg/g/s\n",
                engine.getNetworkReactions(*ctx).size(), rhs->nuclearEnergyGenerationRate);
    return 0;
}
