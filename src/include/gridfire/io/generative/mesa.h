#pragma once
#include "fourdst/atomic/atomicSpecies.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine_abstract.h"

#include <format>

namespace gridfire::io::generative {
    std::string get_mesa_iso_name(const fourdst::atomic::Species& species);

    bool is_proton(const fourdst::atomic::Species& species);
    bool is_alpha(const fourdst::atomic::Species& species);
    bool is_neutron(const fourdst::atomic::Species& species);

    std::string get_mesa_reaction_name(const reaction::Reaction& reaction);
    std::string export_engine_to_mesa_net(const engine::DynamicEngine& engine, engine::scratch::StateBlob& ctx, bool skip_weak);
}
