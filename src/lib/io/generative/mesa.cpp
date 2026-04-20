#include "gridfire/io/generative/mesa.h"

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/reaction/reaction.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "gridfire/utils/config.h"

#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cctype>

namespace gridfire::io::generative {
    std::string get_mesa_iso_name(const fourdst::atomic::Species& species) {
        auto name = std::string(species.name());
        std::ranges::transform(name, name.begin(), ::tolower);

        name.erase(std::ranges::remove(name, '-').begin(), name.end());

        if (name == "p") return "h1";
        if (name == "n" || name == "n1") return "neut";
        if (name == "d") return "h2";
        if (name == "t") return "h3";
        if (name == "a") return "he4";

        return name;
    }

    bool is_proton(const fourdst::atomic::Species& s) { return get_mesa_iso_name(s) == "h1"; }
    bool is_alpha(const fourdst::atomic::Species& s) { return get_mesa_iso_name(s) == "he4"; }
    bool is_neutron(const fourdst::atomic::Species& s) { return get_mesa_iso_name(s) == "neut"; }


    std::string get_mesa_reaction_name(const reaction::Reaction& reaction) {
        std::vector<fourdst::atomic::Species> react_sorted = reaction.reactants();
        std::vector<fourdst::atomic::Species> prod_sorted = reaction.products();

        auto sort_species = [](std::vector<fourdst::atomic::Species>& list) {
            std::ranges::sort(list, [](const auto& a, const auto& b) {
                if (a.z() != b.z()) return a.z() < b.z();
                return a.a() < b.a();
            });
        };

        sort_species(react_sorted);
        sort_species(prod_sorted);

        if (react_sorted.size() == 1 && prod_sorted.size() == 1) {
            if (reaction.type() == reaction::ReactionType::WEAK ||
                reaction.type() == reaction::ReactionType::REACLIB_WEAK ||
                reaction.type() == reaction::ReactionType::LOGICAL_REACLIB_WEAK) {
                return "r_" + get_mesa_iso_name(react_sorted[0]) + "_wk_" + get_mesa_iso_name(prod_sorted[0]);
            }
        }

        if (react_sorted.size() == 2 && prod_sorted.size() == 1) {
            std::string x, cap;
            if (is_proton(react_sorted[0]) || is_proton(react_sorted[1])) {
                cap = "pg";
                x = is_proton(react_sorted[0]) ? get_mesa_iso_name(react_sorted[1]) : get_mesa_iso_name(react_sorted[0]);
            }
            else if (is_alpha(react_sorted[0]) || is_alpha(react_sorted[1])) {
                cap = "ag";
                x = is_alpha(react_sorted[0]) ? get_mesa_iso_name(react_sorted[1]) : get_mesa_iso_name(react_sorted[0]);
            }
            else if (is_neutron(react_sorted[0]) || is_neutron(react_sorted[1])) {
                cap = "ng";
                x = is_neutron(react_sorted[0]) ? get_mesa_iso_name(react_sorted[1]) : get_mesa_iso_name(react_sorted[0]);
            }

            if (!cap.empty()) return "r_" + x + "_" + cap + "_" + get_mesa_iso_name(prod_sorted[0]);
        }

        if (react_sorted.size() == 1 && prod_sorted.size() == 2) {
            std::string x, em;
            if (is_proton(prod_sorted[0]) || is_proton(prod_sorted[1])) {
                em = "gp";
                x = is_proton(prod_sorted[0]) ? get_mesa_iso_name(prod_sorted[1]) : get_mesa_iso_name(prod_sorted[0]);
            }
            else if (is_alpha(prod_sorted[0]) || is_alpha(prod_sorted[1])) {
                em = "ga";
                x = is_alpha(prod_sorted[0]) ? get_mesa_iso_name(prod_sorted[1]) : get_mesa_iso_name(prod_sorted[0]);
            }
            else if (is_neutron(prod_sorted[0]) || is_neutron(prod_sorted[1])) {
                em = "gn";
                x = is_neutron(prod_sorted[0]) ? get_mesa_iso_name(prod_sorted[1]) : get_mesa_iso_name(prod_sorted[0]);
            }

            if (!em.empty()) return "r_" + get_mesa_iso_name(react_sorted[0]) + "_" + em + "_" + x;
        }

        if (react_sorted.size() == 2 && prod_sorted.size() == 2) {
            int r_p = -1, r_a = -1, r_n = -1;
            int p_p = -1, p_a = -1, p_n = -1;

            for(int i=0; i<2; ++i) {
                if(is_proton(react_sorted[i])) r_p = i;
                if(is_alpha(react_sorted[i])) r_a = i;
                if(is_neutron(react_sorted[i])) r_n = i;

                if(is_proton(prod_sorted[i])) p_p = i;
                if(is_alpha(prod_sorted[i])) p_a = i;
                if(is_neutron(prod_sorted[i])) p_n = i;
            }

            std::string x, y, exc;
            if (r_a != -1 && p_p != -1) { exc = "ap"; x = get_mesa_iso_name(react_sorted[1-r_a]); y = get_mesa_iso_name(prod_sorted[1-p_p]); }
            else if (r_p != -1 && p_a != -1) { exc = "pa"; x = get_mesa_iso_name(react_sorted[1-r_p]); y = get_mesa_iso_name(prod_sorted[1-p_a]); }
            else if (r_n != -1 && p_p != -1) { exc = "np"; x = get_mesa_iso_name(react_sorted[1-r_n]); y = get_mesa_iso_name(prod_sorted[1-p_p]); }
            else if (r_p != -1 && p_n != -1) { exc = "pn"; x = get_mesa_iso_name(react_sorted[1-r_p]); y = get_mesa_iso_name(prod_sorted[1-p_n]); }
            else if (r_n != -1 && p_a != -1) { exc = "na"; x = get_mesa_iso_name(react_sorted[1-r_n]); y = get_mesa_iso_name(prod_sorted[1-p_a]); }
            else if (r_a != -1 && p_n != -1) { exc = "an"; x = get_mesa_iso_name(react_sorted[1-r_a]); y = get_mesa_iso_name(prod_sorted[1-p_n]); }

            if (!exc.empty()) return "r_" + x + "_" + exc + "_" + y;
        }

        std::string fallback = "r";
        for (const auto& s : react_sorted) fallback += "_" + get_mesa_iso_name(s);
        fallback += "_to";
        for (const auto& s : prod_sorted) fallback += "_" + get_mesa_iso_name(s);

        return fallback;
    }

    std::string export_engine_to_mesa_net(const engine::DynamicEngine& engine, engine::scratch::StateBlob& ctx, bool skip_weak) {
        std::stringstream ss;
        ss << "! Auto-generated MESA .net file from GridFire\n";
        ss << "! Generated by GridFire version: " << version().toString() << "\n";
        ss << "! Generated on " << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << "\n\n";

        ss << "add_isos(\n";
        for (const auto& species : engine.getNetworkSpecies(ctx)) {
            ss << "    " << get_mesa_iso_name(species) << "\n";
        }
        ss << ")\n\n";

        ss << "add_reactions(\n";
        const auto& reactions = engine.getNetworkReactions(ctx);

        for (const auto& reaction_ptr : reactions) {
            if (skip_weak && (reaction_ptr->type() == reaction::ReactionType::WEAK ||
                reaction_ptr->type() == reaction::ReactionType::REACLIB_WEAK ||
                reaction_ptr->type() == reaction::ReactionType::LOGICAL_REACLIB_WEAK)) {
                continue;
            }
            ss << "    " << get_mesa_reaction_name(*reaction_ptr) << "\n";
        }
        ss << ")\n";

        return ss.str();
    }

}
