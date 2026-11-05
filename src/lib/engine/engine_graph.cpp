#include "gridfire/engine/engine_graph.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/types/types.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/engine/procedures/priming.h"
#include "gridfire/partition/partition_ground.h"
#include "gridfire/engine/procedures/construction.h"
#include "gridfire/utils/hashing.h"
#include "gridfire/utils/table_format.h"

#include "gridfire/engine/scratchpads/engine_graph_scratchpad.h"
#include "gridfire/engine/scratchpads/blob.h"
#include "gridfire/engine/scratchpads/utils.h"

#include "fourdst/atomic/species.h"
#include "fourdst/atomic/atomicSpecies.h"

#include "quill/LogMacros.h"

// ReSharper disable once CppUnusedIncludeDirective
#include <algorithm>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <fstream>
#include <ranges>

#include "cppad/cppad.hpp"
#include "cppad/utility/sparse_rc.hpp"
#include "cppad/utility/sparse_rcv.hpp"
#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "gridfire/reaction/reaclib.h"
#include "gridfire/reaction/weak/neutrino_energies.h"


namespace {
    enum class REACLIB_WEAK_TYPES {
        BETA_PLUS_DECAY,
        BETA_MINUS_DECAY,
        ELECTRON_CAPTURE,
        POSITRON_CAPTURE,
        NONE
    };

    bool is_reaclib_weak(const gridfire::reaction::Reaction& r) {
        return r.type() == gridfire::reaction::ReactionType::REACLIB_WEAK ||
               r.type() == gridfire::reaction::ReactionType::LOGICAL_REACLIB_WEAK;
    }

    REACLIB_WEAK_TYPES get_weak_reaclib_reaction_type(const gridfire::reaction::Reaction& r) {
        if (!is_reaclib_weak(r)) {
            return REACLIB_WEAK_TYPES::NONE;
        }

        // Get the () part of the id
        const std::string_view id = r.id();
        const size_t open_paren_pos = id.find('(');
        const size_t close_paren_pos = id.find(')');
        if (open_paren_pos == std::string_view::npos || close_paren_pos == std::string_view::npos || close_paren_pos <= open_paren_pos) {
            throw gridfire::exceptions::ReactionParsingError("Invalid REACLIB weak reaction ID format.", std::string(id));
        }

        const std::string_view reaction_type_str = id.substr(open_paren_pos + 1, close_paren_pos - open_paren_pos - 1);

        // Find the comma and extract the part to the left and right of it
        const size_t reaction_type_pos = reaction_type_str.find(',');
        if (reaction_type_pos == std::string_view::npos) {
            throw gridfire::exceptions::ReactionParsingError("Invalid REACLIB weak reaction ID format: missing comma.", std::string(id));
        }

        const std::string_view projectiles_str = reaction_type_str.substr(0, reaction_type_pos);
        const std::string_view ejectiles_str = reaction_type_str.substr(reaction_type_pos + 1);

        // Check if the projectiles string has "e+" or "e-"
        const bool has_captured_positron = (projectiles_str.find("e+") != std::string_view::npos);
        const bool has_captured_electron = (projectiles_str.find("e-") != std::string_view::npos);

        const bool has_ejected_electron = (ejectiles_str.find("e-") != std::string_view::npos);
        const bool has_ejected_positron = (ejectiles_str.find("e+") != std::string_view::npos);

        // Assert that only one of the four possibilities is true
        const int true_count = static_cast<int>(has_captured_positron) +
                         static_cast<int>(has_captured_electron) +
                         static_cast<int>(has_ejected_electron) +
                         static_cast<int>(has_ejected_positron);
        if (true_count != 1) {
            throw gridfire::exceptions::ReactionParsingError("Invalid REACLIB weak reaction ID format: must have exactly one of e+, e- in projectiles or ejectiles.", std::string(id));
        }

        if (has_ejected_positron) {
            return REACLIB_WEAK_TYPES::BETA_PLUS_DECAY;
        } if (has_ejected_electron) {
            return REACLIB_WEAK_TYPES::BETA_MINUS_DECAY;
        } if (has_captured_electron) {
            return REACLIB_WEAK_TYPES::ELECTRON_CAPTURE;
        } if (has_captured_positron) {
            return REACLIB_WEAK_TYPES::POSITRON_CAPTURE;
        }

        return REACLIB_WEAK_TYPES::NONE;
    }
}

namespace gridfire::engine {
    namespace {
        /// The composition a network is built from, without the passive species (which must not seed reactions)
        fourdst::composition::Composition without_passive_species(
            const fourdst::composition::Composition& composition,
            const config::GridFireConfig& config
        ) {
            if (config.engine.passiveSpecies.empty()) {
                return composition;
            }
            fourdst::composition::Composition seeds;
            for (const auto& species : composition.getRegisteredSpecies()) {
                if (std::ranges::find(config.engine.passiveSpecies, std::string(species.name())) == config.engine.passiveSpecies.end()) {
                    seeds.registerSpecies(species);
                    seeds.setMolarAbundance(species, composition.getMolarAbundance(species));
                }
            }
            return seeds;
        }
    }

    GraphEngine::GraphEngine(
        const fourdst::composition::Composition &composition,
        const BuildDepthType buildDepth
    ): GraphEngine(composition, partition::GroundStatePartitionFunction(), buildDepth) {}

    GraphEngine::GraphEngine(
        const fourdst::composition::Composition &composition,
        const partition::PartitionFunction& partitionFunction,
        const BuildDepthType buildDepth
    ) : GraphEngine(composition, partitionFunction, buildDepth, NetworkConstructionFlags::DEFAULT){}

    GraphEngine::GraphEngine(
        const fourdst::composition::Composition &composition,
        const partition::PartitionFunction &partitionFunction,
        const BuildDepthType buildDepth,
        const NetworkConstructionFlags reactionTypes
    ) : GraphEngine(composition, partitionFunction, buildDepth, reactionTypes, config::GridFireConfig{}) {}

    GraphEngine::GraphEngine(
        const reaction::ReactionSet &reactions
    ) : GraphEngine(reactions, config::GridFireConfig{}) {}

    GraphEngine::GraphEngine(
        const reaction::ReactionSet &reactions,
        const partition::PartitionFunction &partitionFunction
    ) : GraphEngine(reactions, partitionFunction, config::GridFireConfig{}) {}

    GraphEngine::GraphEngine(
        const fourdst::composition::Composition &composition,
        const config::GridFireConfig &config
    ) : GraphEngine(composition, partition::GroundStatePartitionFunction(), NetworkBuildDepth::Full, NetworkConstructionFlags::DEFAULT, config) {}

    GraphEngine::GraphEngine(
        const fourdst::composition::Composition &composition,
        const partition::PartitionFunction &partitionFunction,
        const BuildDepthType buildDepth,
        const NetworkConstructionFlags reactionTypes,
        const config::GridFireConfig &config
    ) :
    m_weakRateInterpolator(rates::weak::UNIFIED_WEAK_DATA),
    m_reactions(build_nuclear_network(without_passive_species(composition, config), m_weakRateInterpolator, buildDepth, reactionTypes)),
    m_partitionFunction(partitionFunction.clone()),
    m_depth(buildDepth)
    {
        configure(config);
        syncInternalMaps();
    }

    GraphEngine::GraphEngine(
        const reaction::ReactionSet &reactions,
        const config::GridFireConfig &config
    ) :
    m_weakRateInterpolator(rates::weak::UNIFIED_WEAK_DATA),
    m_reactions(reactions)
    {
        configure(config);
        syncInternalMaps();
    }

    GraphEngine::GraphEngine(
        const reaction::ReactionSet &reactions,
        const partition::PartitionFunction &partitionFunction,
        const config::GridFireConfig &config
    ) :
    m_weakRateInterpolator(rates::weak::UNIFIED_WEAK_DATA),
    m_reactions(reactions),
    m_partitionFunction(partitionFunction.clone())
    {
        configure(config);
        syncInternalMaps();
    }

    void GraphEngine::configure(const config::GridFireConfig &config) {
        m_config.mutate([&config](auto& cfg) {
            cfg = config;
        });
        m_screeningType = m_config->engine.screening;
        m_screeningModel = screening::selectScreeningModel(m_screeningType, m_config->engine.screeningBlend);
        m_jacobianMethod = m_config->engine.jacobian;
    }

    void GraphEngine::addReaction(
        const reaction::Reaction& reaction
    ) {
        m_reactions.add_reaction(reaction);
        syncInternalMaps();
    }

    void GraphEngine::addReaction(
        const std::string& reaction_id
    ) {
        const auto& allReaclibReactions = reaclib::get_all_reaclib_reactions();
        const auto& reaction = allReaclibReactions.get(reaction_id);
        if (reaction.has_value()) {
            m_reactions.add_reaction(reaction.value()->clone());
        } else {
            throw exceptions::BadCollectionError(std::format("Unable to locate reaction with ID {} in reaclib set", reaction_id));
        }
    }

    std::unique_ptr<scratch::StateBlob> GraphEngine::constructStateBlob(const scratch::StateBlob *blob) const {
        if (blob) {
            throw exceptions::ScratchPadError("GraphEngine does not support accepting an external StateBlob. The state blob for GraphEngine must be constructed internally to ensure it contains the correct scratchpad states.");
        }
        auto i_blob = std::make_unique<scratch::StateBlob>();
        i_blob->enroll<engine::scratch::GraphEngineScratchPad>();
        auto* state = scratch::get_state<scratch::GraphEngineScratchPad, false>(*i_blob);
        state->initialize(*this);
        return i_blob;
    }

    std::expected<StepDerivatives<double>, EngineStatus> GraphEngine::calculateRHSAndEnergy(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        bool trust
    ) const {
        return calculateRHSAndEnergy(ctx, comp, T9, rho, m_reactions);
    }

    std::expected<StepDerivatives<double>, EngineStatus> GraphEngine::calculateRHSAndEnergy(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const reaction::ReactionSet &activeReactions
    ) const {
        auto* state = scratch::get_state<scratch::GraphEngineScratchPad, true>(ctx);
        LOG_TRACE_L3(m_logger, "Calculating RHS and Energy in GraphEngine at T9 = {}, rho = {}.", T9, rho);
        const double Ye = comp.getElectronAbundance();
        const std::vector<double> molarAbundances = comp.getMolarAbundanceVector();
        if (m_usePrecomputation) {
            const std::size_t state_hash = utils::hash_state(comp, T9, rho, activeReactions);
            if (state->stepDerivativesCache.contains(state_hash)) {
                const StepDerivatives<double>& cached = state->stepDerivativesCache.at(state_hash);
                state->most_recent_rhs_calculation = cached;
                return cached;
            }
            LOG_TRACE_L3(m_logger, "Using precomputation for reaction rates in GraphEngine calculateRHSAndEnergy.");
            std::vector<double> bare_rates;
            std::vector<double> bare_reverse_rates;
            compute_bare_rates(comp, T9, rho, activeReactions, bare_rates, bare_reverse_rates);

            LOG_TRACE_L3(m_logger, "Precomputed {} forward and {} reverse reaction rates for active reactions.", bare_rates.size(), bare_reverse_rates.size());

            // --- The public facing interface can always use the precomputed version since taping is done internally ---
            StepDerivatives<double> result =  calculateAllDerivativesUsingPrecomputation(ctx, comp, bare_rates, bare_reverse_rates, T9, rho, activeReactions);
            state->stepDerivativesCache.insert(std::make_pair(state_hash, result));
            state->most_recent_rhs_calculation = result;
            return result;
        } else {
            LOG_TRACE_L2(m_logger, "Not using precomputation for reaction rates in GraphEngine calculateRHSAndEnergy.");
            StepDerivatives<double> result = calculateAllDerivatives<double>(
                molarAbundances,
                T9,
                rho,
                Ye,
                0.0,
                [&comp](const fourdst::atomic::Species& species) -> std::optional<size_t> {
                    if (comp.contains(species)) {
                        return comp.getSpeciesIndex(species); // Return the index of the species in the composition
                    }
                    return std::nullopt; // Species not found in the composition
                },
                [&activeReactions](const reaction::Reaction& reaction) -> bool {
                    if (activeReactions.contains(reaction)) { return true; }
                    return false;
                }
            );
            state->most_recent_rhs_calculation = result;
            return result;
        }
    }

    void GraphEngine::compute_bare_rates(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const reaction::ReactionSet &activeReactions,
        std::vector<double> &bare_rates,
        std::vector<double> &bare_reverse_rates
    ) const {
        const double Ye = comp.getElectronAbundance();
        const std::vector<double> molarAbundances = comp.getMolarAbundanceVector();
        bare_rates.clear();
        bare_reverse_rates.clear();
        bare_rates.reserve(activeReactions.size());
        bare_reverse_rates.reserve(activeReactions.size());

        for (const auto& reaction: activeReactions) {
            assert(m_reactions.contains(*reaction)); // A bug which results in this failing indicates a serious internal inconsistency and should only be present during development.
            bare_rates.push_back(reaction->calculate_rate(T9, rho, Ye, 0.0, molarAbundances, m_indexToSpeciesMap));
            // Always push so bare_reverse_rates stays index-aligned with activeReactions (the kernel indexes it that way)
            bare_reverse_rates.push_back(has_reverse_rate(*reaction) ? calculateReverseRate(*reaction, T9, rho, comp) : 0.0);
        }
    }

    std::vector<double> GraphEngine::calculateEnergyChannels(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        return calculateEnergyChannels(ctx, comp, T9, rho, m_reactions);
    }

    std::vector<double> GraphEngine::calculateEnergyChannels(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const reaction::ReactionSet &activeReactions
    ) const {
        // Uses the same rate evaluation and flow kernel as calculateRHSAndEnergy so that the channels sum to eps
        std::vector<double> bare_rates;
        std::vector<double> bare_reverse_rates;
        compute_bare_rates(comp, T9, rho, activeReactions, bare_rates, bare_reverse_rates);

        // Network-ordered abundances for screening as well as the flows (comp may hold a subset, in its own order)
        const std::vector<double> local_abundances = networkAbundances(comp);
        const std::vector<double> screeningFactors = m_screeningModel->calculateScreeningFactors(
            activeReactions,
            m_networkSpecies,
            local_abundances,
            T9,
            rho
        );

        const PrecomputationKernelResults kernel = accumulate_flows_serial(
            ctx,
            local_abundances,
            screeningFactors,
            bare_rates,
            bare_reverse_rates,
            rho,
            activeReactions
        );

        std::vector<double> eps_by_channel(m_energy_channels.names.size(), 0.0);
        for (const auto& [flow, j] : std::views::zip(kernel.molar_flows, kernel.reaction_indices)) {
            const PrecomputedReaction& precomp = m_precomputed_reactions[j];
            eps_by_channel[precomp.channel_index] += flow * precomp.energy_per_flow;
        }
        return eps_by_channel;
    }

    std::optional<ReactionContributions> GraphEngine::calculateReactionContributions(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        // The same rate evaluation and flow kernel as calculateRHSAndEnergy, so the contributions sum to its totals
        std::vector<double> bare_rates;
        std::vector<double> bare_reverse_rates;
        compute_bare_rates(comp, T9, rho, m_reactions, bare_rates, bare_reverse_rates);
        // Network-ordered abundances for screening as well as the flows (comp may hold a subset, in its own order)
        const std::vector<double> local_abundances = networkAbundances(comp);
        const std::vector<double> screeningFactors = m_screeningModel->calculateScreeningFactors(
            m_reactions, m_networkSpecies, local_abundances, T9, rho
        );
        const PrecomputationKernelResults kernel = accumulate_flows_serial(
            ctx, local_abundances, screeningFactors, bare_rates, bare_reverse_rates, rho, m_reactions
        );

        const size_t n = m_reactions.size();
        ReactionContributions contributions{
            .molar_flow = std::vector<double>(n, 0.0),
            .energy_rate = std::vector<double>(n, 0.0),
            .neutrino_energy_loss_rate = std::vector<double>(n, 0.0),
            .neutrino_rate = std::vector<double>(n, 0.0),
        };
        for (const auto& [flow, j] : std::views::zip(kernel.molar_flows, kernel.reaction_indices)) {
            const PrecomputedReaction& precomp = m_precomputed_reactions[j];
            contributions.molar_flow[j] = flow;
            contributions.energy_rate[j] = flow * precomp.energy_per_flow;
            if (precomp.emits_neutrino) {
                contributions.neutrino_energy_loss_rate[j] = flow * precomp.neutrino_energy_per_flow;
                contributions.neutrino_rate[j] = flow * m_constants.Na;
            }
        }
        return contributions;
    }

    const std::vector<std::string>& GraphEngine::getEnergyChannelNames(
        scratch::StateBlob& ctx
    ) const {
        return m_energy_channels.names;
    }

    EnergyDerivatives GraphEngine::calculateEpsDerivatives(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        return calculateEpsDerivatives(ctx, comp, T9, rho, m_reactions);
    }

    EnergyDerivatives GraphEngine::calculateEpsDerivatives(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const reaction::ReactionSet &activeReactions
    ) const {
        if (m_jacobianMethod == JacobianMethod::ANALYTIC) {
            if (const auto analytic = calculateAnalyticEpsDerivatives(comp, T9, rho)) {
                return *analytic;
            }
        }
        return calculateAutoDiffEpsDerivatives(ctx, comp, T9, rho);
    }

    std::optional<EnergyDerivatives> GraphEngine::calculateAnalyticEpsDerivatives(
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        for (const auto& precomp : m_precomputed_reactions) {
            if (!precomp.analytic_rate_derivative || (m_useReverseReactions && precomp.reverse_symmetry_factor != 0.0)) {
                return std::nullopt;
            }
        }
        const std::vector<double> Y = networkAbundances(comp);
        const std::vector<double> screeningFactors = m_screeningModel->calculateScreeningFactors(
            m_reactions, m_networkSpecies, Y, T9, rho
        );
        const auto dlnf = m_screeningModel->calculateScreeningStateLogDerivatives(
            m_reactions, m_networkSpecies, Y, T9, rho
        );
        if (!dlnf) {
            return std::nullopt;
        }
        std::vector<double> bare_rates;
        std::vector<double> bare_reverse_rates;
        compute_bare_rates(comp, T9, rho, m_reactions, bare_rates, bare_reverse_rates);

        // Unused by REACLIB rates, but the interface takes one
        static const fourdst::composition::Composition no_composition;

        // eps = sum_r F_r E_r at fixed composition, so d eps / d ln x = sum_r F_r E_r d ln F_r / d ln x
        double dEps_dlnT = 0.0;
        double dEps_dlnRho = 0.0;
        for (size_t r = 0; r < m_reactions.size(); ++r) {
            const reaction::Reaction& reaction = m_reactions[r];
            const PrecomputedReaction& precomp = m_precomputed_reactions[r];

            // F = f k sym rho^(n-1) prod_u Y_u^p_u, exactly as the RHS forms it
            const size_t numReactants = reaction.reactants().size();
            const double rhoPower = numReactants > 1 ? static_cast<double>(numReactants) - 1.0 : 0.0;
            double flow = screeningFactors[r] * bare_rates[r] * precomp.symmetry_factor * std::pow(rho, rhoPower);
            for (size_t u = 0; u < precomp.unique_reactant_indices.size(); ++u) {
                for (int p = 0; p < precomp.reactant_powers[u]; ++p) {
                    flow *= Y[precomp.unique_reactant_indices[u]];
                }
            }
            const double weight = flow * precomp.energy_per_flow;
            if (weight == 0.0) {
                continue;
            }

            const double dlnk_dlnT = T9 * reaction.calculate_log_rate_partial_deriv_wrt_T9(T9, rho, 0.0, 0.0, no_composition);
            const double dlnk_dlnRho = precomp.electron_capture ? 1.0 : 0.0; // rate proportional to rho * Ye
            dEps_dlnT += weight * (dlnk_dlnT + (*dlnf)[r].dlnT);
            dEps_dlnRho += weight * (rhoPower + dlnk_dlnRho + (*dlnf)[r].dlnRho);
        }
        return EnergyDerivatives{dEps_dlnT / (T9 * 1e9), dEps_dlnRho / rho};
    }

    EnergyDerivatives GraphEngine::calculateAutoDiffEpsDerivatives(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        auto* state = scratch::get_state<scratch::GraphEngineScratchPad, true>(ctx);
        const size_t numSpecies = m_networkSpecies.size();
        const size_t numADInputs = numSpecies + 2; // +2 for T9 and rho

        if (comp.getRegisteredSpecies().size() != numSpecies) {
            LOG_ERROR(m_logger, "Input abundance vector size ({}) does not match number of species in the network ({}).",
                      comp.getRegisteredSpecies().size(), numSpecies);
            throw std::invalid_argument("Input abundance vector size does not match number of species in the network.");
        }

        std::vector<double> x(numADInputs);
        const std::vector<double> Y = comp.getMolarAbundanceVector();
        for (size_t i = 0; i < numSpecies; ++i) {
            x[i] = Y[i];
        }

        x[numSpecies] = T9;
        x[numSpecies + 1] = rho;

        // Use reverse mode to get the gradient. W selects which dependent variable we care about, the Eps AD tape only has eps as a dependent variable so we just select set the 0th element to 1.
        std::vector<double> w(numSpecies + 1, 0.0);
        w[numSpecies] = 1.0; // We want the derivative of the energy generation rate

        // Sweep the tape forward to record the function value at x
        assert(state->rhsADFun.has_value() && "AD tape for energy derivatives has not been recorded.");
        state->rhsADFun.value().Forward(0, x);

        // Extract the gradient at the previously evaluated point x using reverse mode
        const std::vector<double> eps_derivatives = state->rhsADFun.value().Reverse(1, w);

        const double dEps_dT9 = eps_derivatives[numSpecies];
        const double dEps_dRho = eps_derivatives[numSpecies + 1];

        // Chain rule to scale from deps/dT9 to deps/dT
        // dT9/dT = 1e-9
        const double dEps_dT = dEps_dT9 * 1e-9;


        return {dEps_dT, dEps_dRho};
    }

    void GraphEngine::generate_jacobian_sparsity_pattern() {
        const size_t outputSize = m_authoritativeADFun.Range();

        // Create a range x range identity pattern
        CppAD::sparse_rc<std::vector<size_t>> patternIn(outputSize, outputSize, outputSize);
        for (size_t i = 0; i < outputSize; ++i) {
            patternIn.set(i, i, i);
        }

        m_authoritativeADFun.rev_jac_sparsity(patternIn, false, false, false, m_full_jacobian_sparsity_pattern);

        m_full_sparsity_set.clear();
        const auto& rows = m_full_jacobian_sparsity_pattern.row();
        const auto& cols = m_full_jacobian_sparsity_pattern.col();
        const size_t nnz = m_full_jacobian_sparsity_pattern.nnz();

        for (size_t k = 0; k < nnz; ++k) {
            if (cols[k] < m_networkSpecies.size() + 1) {
                m_full_sparsity_set.insert(std::make_pair(rows[k], cols[k]));
            }
        }
    }

    void GraphEngine::syncInternalMaps() {
        LOG_INFO(m_logger, "Synchronizing internal maps for REACLIB graph network (serif::network::GraphNetwork)...");
        collectNetworkSpecies();
        populateReactionIDMap();
        populateSpeciesToIndexMap();
        collectAtomicReverseRateAtomicBases();

        recordADTape(); // Record the AD tape for the RHS of the ODE (dY/di and dEps/di) for all independent variables i
        generate_jacobian_sparsity_pattern();

        precomputeNetwork();
        LOG_INFO(m_logger, "Internal maps synchronized. Network contains {} species and {} reactions.",
                 m_networkSpecies.size(), m_reactions.size());
    }

    // --- Network Graph Construction Methods ---
    void GraphEngine::collectNetworkSpecies() {
        m_networkSpecies.clear();
        m_networkSpeciesMap.clear();

        std::set<std::string_view> uniqueSpeciesNames;

        for (const auto& reaction: m_reactions) {
            for (const auto& reactant: reaction->reactants()) {
                uniqueSpeciesNames.insert(reactant.name());
            }
            for (const auto& product: reaction->products()) {
                uniqueSpeciesNames.insert(product.name());
            }
        }

        // Passive species: in the composition (Y_e, screening) but in no reaction, so their abundances never change
        for (const auto& name : m_config->engine.passiveSpecies) {
            if (!fourdst::atomic::species.contains(name)) {
                throw exceptions::BadCollectionError(std::format("Passive species {} is not a known species.", name));
            }
            if (uniqueSpeciesNames.contains(fourdst::atomic::species.at(name).name())) {
                throw exceptions::BadCollectionError(std::format(
                    "Passive species {} takes part in reactions of the network; only species outside it can be passive.", name));
            }
            uniqueSpeciesNames.insert(fourdst::atomic::species.at(name).name());
        }

        for (const auto& name: uniqueSpeciesNames) {
            auto it = fourdst::atomic::species.find(std::string(name));
            if (it != fourdst::atomic::species.end()) {
                m_networkSpecies.push_back(it->second);
                m_networkSpeciesMap.insert({name, it->second});
            } else {
                LOG_ERROR(m_logger, "Species '{}' not found in global atomic species database.", name);
                m_logger->flush_log();
                throw std::runtime_error("Species not found in global atomic species database: " + std::string(name));
            }
        }
        std::ranges::sort(m_networkSpecies, [](const fourdst::atomic::Species& a, const fourdst::atomic::Species& b) -> bool {
            return a.mass() < b.mass(); // Otherwise, sort by mass
        });
    }

    void GraphEngine::populateReactionIDMap() {
        LOG_TRACE_L1(m_logger, "Populating reaction ID map for REACLIB graph network (serif::network::GraphNetwork)...");
        m_reactionIDMap.clear();
        for (auto& reaction: m_reactions) {
            m_reactionIDMap.emplace(reaction->id(), reaction.get());
        }
        LOG_TRACE_L1(m_logger, "Populated {} reactions in the reaction ID map.", m_reactionIDMap.size());
    }

    void GraphEngine::populateSpeciesToIndexMap() {
        m_speciesToIndexMap.clear();
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            m_speciesToIndexMap.insert({m_networkSpecies[i], i});
        }
        m_indexToSpeciesMap.clear();
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            m_indexToSpeciesMap.insert({i, m_networkSpecies[i]});
        }
    }

    // --- Basic Accessors and Queries ---
    const std::vector<fourdst::atomic::Species>& GraphEngine::getNetworkSpecies(scratch::StateBlob &ctx) const {
        return m_networkSpecies;
    }

    const reaction::ReactionSet& GraphEngine::getNetworkReactions(
        scratch::StateBlob& ctx
    ) const {
        return m_reactions;
    }

    bool GraphEngine::involvesSpecies(
        scratch::StateBlob& ctx,
        const fourdst::atomic::Species& species
    ) const {
        const bool found = m_networkSpeciesMap.contains(species.name());
        return found;
    }

    double GraphEngine::compute_reaction_flow(
        scratch::StateBlob& ctx,
        const std::vector<double> &local_abundances,
        const std::vector<double> &screening_factors,
        const std::vector<double> &bare_rates,
        const std::vector<double> &bare_reverse_rates,
        const double rho,
        const size_t reactionCounter,
        const reaction::Reaction &reaction,
        const size_t reactionIndex,
        const PrecomputedReaction &precomputedReaction
    ) const {
        double forwardAbundanceProduct = 1.0;
        for (size_t i = 0; i < precomputedReaction.unique_reactant_indices.size(); ++i) {
            const size_t reactantIndex = precomputedReaction.unique_reactant_indices[i];
            const int power = precomputedReaction.reactant_powers[i];

            const double abundance = local_abundances[reactantIndex];

            double factor;
            if (power == 1) { factor = abundance; }
            else if (power == 2)  { factor = abundance * abundance; }
            else if (power == 3)  { factor = abundance * abundance * abundance; }
            else { factor = std::pow(abundance, static_cast<double>(power)); }

            if (!std::isfinite(factor)) {
                const auto& sp = m_indexToSpeciesMap.at(reactantIndex);
                std::string error_msg = std::format("Non-finite factor encountered in forward abundance product in reaction {} for species {} (Abundance: {}). Check input abundances for validity.", reaction.id(), sp.name(), abundance);
                LOG_CRITICAL(m_logger, "{}", error_msg);
                throw exceptions::BadRHSEngineError(error_msg);
            }

            forwardAbundanceProduct *= factor;
        }

        const double bare_rate = bare_rates.at(reactionCounter);

        const double screeningFactor = screening_factors[reactionCounter];
        const size_t numReactants = m_reactions[reactionIndex].reactants().size();
        const size_t numProducts = m_reactions[reactionIndex].products().size();

        const double forwardMolarReactionFlow = screeningFactor *
                                   bare_rate *
                                   precomputedReaction.symmetry_factor *
                                   forwardAbundanceProduct *
                                   std::pow(rho, numReactants >  1 ? static_cast<double>(numReactants) - 1 : 0.0);
        if (!std::isfinite(forwardMolarReactionFlow)) {
            LOG_CRITICAL(m_logger, "Non-finite forward molar reaction flow computed for reaction '{}'. Check input abundances and rates for validity.", reaction.id());
            throw exceptions::BadRHSEngineError("Non-finite forward molar reaction flow computed.");
        }


        double reverseMolarReactionFlow = 0.0;
        if (precomputedReaction.reverse_symmetry_factor != 0.0 and m_useReverseReactions) {
            const double bare_reverse_rate = bare_reverse_rates.at(reactionCounter);

            double reverseAbundanceProduct = 1.0;
            for (size_t i = 0; i < precomputedReaction.unique_product_indices.size(); ++i) {
                const size_t productIndex = precomputedReaction.unique_product_indices[i];
                reverseAbundanceProduct *= std::pow(local_abundances[productIndex], precomputedReaction.product_powers[i]);
            }

            reverseMolarReactionFlow = screeningFactor *
                                       bare_reverse_rate *
                                       precomputedReaction.reverse_symmetry_factor *
                                       reverseAbundanceProduct *
                                       std::pow(rho, numProducts > 1 ? static_cast<double>(numProducts) - 1 : 0.0);
        }

        return forwardMolarReactionFlow - reverseMolarReactionFlow;
    }

    double GraphEngine::neutrino_energy_loss_per_molar_flow(
        const reaction::Reaction &reaction
    ) const {
        if (!is_reaclib_weak(reaction)) {
            return 0.0;
        }
        // Prefer tabulated average neutrino energies (Bahcall 2002 / weaklib, as used by MESA). The fraction-of-Q
        // heuristic below is only a fallback: it overestimates losses for e.g. p(p,e+)d (0.72 vs 0.267 MeV).
        if (const auto tabulated = reaction::weak::tabulated_neutrino_energy_MeV(reaction)) {
            return *tabulated * m_constants.Na * m_constants.MeV_to_erg;
        }
        const double q_abs = std::abs(reaction.qValue());
        const REACLIB_WEAK_TYPES weakType = get_weak_reaclib_reaction_type(reaction);
        double neutrino_loss_fraction = 0.0;
        switch (weakType) {
            case REACLIB_WEAK_TYPES::BETA_PLUS_DECAY:
                [[fallthrough]];
            case REACLIB_WEAK_TYPES::BETA_MINUS_DECAY:
                neutrino_loss_fraction = 0.5; // Approximate 50% energy loss to neutrinos for beta decays
                break;
            case REACLIB_WEAK_TYPES::ELECTRON_CAPTURE:
                [[fallthrough]];
            case REACLIB_WEAK_TYPES::POSITRON_CAPTURE:
                neutrino_loss_fraction = 1.0;
                break;
            default: ;
        }
        return q_abs * neutrino_loss_fraction * m_constants.Na * m_constants.MeV_to_erg;
    }

    std::pair<double, double> GraphEngine::compute_neutrino_fluxes(
        const double netFlow,
        const PrecomputedReaction &precomputedReaction
    ) const {
        if (!precomputedReaction.emits_neutrino) {
            return {0.0, 0.0};
        }
        return {netFlow * precomputedReaction.neutrino_energy_per_flow, netFlow * m_constants.Na};
    }

    GraphEngine::PrecomputationKernelResults GraphEngine::accumulate_flows_serial(
        scratch::StateBlob& ctx,
        const std::vector<double> &local_abundances,
        const std::vector<double> &screening_factors,
        const std::vector<double> &bare_rates,
        const std::vector<double> &bare_reverse_rates,
        const double rho,
        const reaction::ReactionSet &activeReactions
    ) const {
        PrecomputationKernelResults results;
        results.dydt_vector.resize(m_networkSpecies.size(), 0.0);

        std::vector<double> molarReactionFlows;
        molarReactionFlows.reserve(m_precomputed_reactions.size());

        size_t reactionCounter = 0;
        std::vector<size_t> reactionIndices;
        reactionIndices.reserve(m_precomputed_reactions.size());

        for (const auto& reaction : activeReactions) {
            uint64_t reactionHash = reaction->hash(0);
            const size_t reactionIndex = m_precomputed_reaction_index_map.at(reactionHash);
            reactionIndices.push_back(reactionIndex);
            const PrecomputedReaction& precomputedReaction = m_precomputed_reactions[reactionIndex];

            double netFlow = compute_reaction_flow(
                ctx,
                local_abundances,
                screening_factors,
                bare_rates,
                bare_reverse_rates,
                rho,
                reactionCounter,
                *reaction,
                reactionIndex,
                precomputedReaction);

            molarReactionFlows.push_back(netFlow);

            auto [local_neutrino_loss, local_neutrino_flux] = compute_neutrino_fluxes(netFlow, precomputedReaction);
            results.total_neutrino_energy_loss_rate += local_neutrino_loss;
            results.total_neutrino_flux += local_neutrino_flux;

            reactionCounter++;
        }

        LOG_TRACE_L3(m_logger, "Computed {} molar reaction flows for active reactions. Assembling these into RHS", molarReactionFlows.size());

        reactionCounter = 0;
        for (const auto& [reaction, j]: std::views::zip(activeReactions, reactionIndices)) {
            const auto& precomp = m_precomputed_reactions[j];
            const double R_j = molarReactionFlows[reactionCounter];

            for (size_t i = 0; i < precomp.affected_species_indices.size(); ++i) {
                const size_t speciesIndex = precomp.affected_species_indices[i];

                const int stoichiometricCoefficient = precomp.stoichiometric_coefficients[i];

                const double dydt_increment = static_cast<double>(stoichiometricCoefficient) * R_j;
                results.dydt_vector[speciesIndex] += dydt_increment;
            }
            reactionCounter++;
        }

        results.molar_flows = std::move(molarReactionFlows);
        results.reaction_indices = std::move(reactionIndices);
        return results;
    }

    double GraphEngine::calculateReverseRate(
        const reaction::Reaction &reaction,
        const double T9,
        const double rho,
        const fourdst::composition::CompositionAbstract &comp
    ) const {
        if (!m_useReverseReactions) {
            LOG_TRACE_L3_LIMIT_EVERY_N(std::numeric_limits<int>::max(), m_logger, "Reverse reactions are disabled. Returning 0.0 for reverse rate of reaction '{}'.", reaction.id());
            return 0.0; // If reverse reactions are not used, return 0.0
        }
        const double temp = T9 * 1e9; // Convert T9 to Kelvin

        // Reverse reactions are only relevant for strong reactions (at least during the vast majority of stellar evolution)
        // So here we just let these be dummy values since we know
        // 1. The reaction should always be strong
        // 2. The strong reaction rate is independent of Ye and mue
        //
        // In development builds the assert below will confirm this
        constexpr double Ye = 0.0;
        constexpr double mue = 0.0;

        // It is a logic error to call this function on a weak reaction
        assert(reaction.type() != gridfire::reaction::ReactionType::WEAK);

        // In debug builds we check the units on kB to ensure it is in erg/K. This is removed in release builds to avoid overhead. (Note assert is a no-op in release builds)
        assert(Constants::getInstance().get("kB").unit == "erg / K");

        const double kBMeV = m_constants.kB * 624151; // Convert kB to MeV/K NOTE: This relies on the fact that m_constants.kB is in erg/K!
        const double expFactor = std::exp(-reaction.qValue() / (kBMeV * temp));
        double reverseRate = 0.0;
        // We also let Y be an empy vector since the strong reaction rate is independent of Y
        const double forwardRate = reaction.calculate_rate(T9, rho, Ye, mue, {}, m_indexToSpeciesMap);

        if (reaction.reactants().size() == 2 && reaction.products().size() == 2) {
            reverseRate = calculateReverseRateTwoBody(reaction, T9, forwardRate, expFactor);
        } else {
            LOG_WARNING_LIMIT_EVERY_N(1000000, m_logger, "Reverse rate calculation for reactions with more than two reactants or products is not implemented (reaction id {}).", reaction.id());
        }
        LOG_TRACE_L2_LIMIT_EVERY_N(1000, m_logger, "Calculated reverse rate for reaction '{}': {:.3E} at T9={:.3E}.", reaction.id(), reverseRate, T9);
        return reverseRate;
    }

    double GraphEngine::calculateReverseRateTwoBody(
        const reaction::Reaction &reaction,
        const double T9,
        const double forwardRate,
        const double expFactor
    ) const {
        std::vector<double> reactantPartitionFunctions;
        std::vector<double> productPartitionFunctions;

        reactantPartitionFunctions.reserve(reaction.reactants().size());
        productPartitionFunctions.reserve(reaction.products().size());

        std::unordered_map<fourdst::atomic::Species, int> reactantMultiplicity;
        std::unordered_map<fourdst::atomic::Species, int> productMultiplicity;

        reactantMultiplicity.reserve(reaction.reactants().size());
        productMultiplicity.reserve(reaction.products().size());

        for (const auto& reactant : reaction.reactants()) {
            reactantMultiplicity[reactant] += 1;
        }
        for (const auto& product : reaction.products()) {
            productMultiplicity[product] += 1;
        }
        double reactantSymmetryFactor = 1.0;
        double productSymmetryFactor = 1.0;
        for (const auto& count : reactantMultiplicity | std::views::values) {
            reactantSymmetryFactor *= std::tgamma(count + 1);
        }
        for (const auto& count : productMultiplicity | std::views::values) {
            productSymmetryFactor *= std::tgamma(count + 1);
        }
        const double symmetryFactor = reactantSymmetryFactor / productSymmetryFactor;

        // Accumulate mass terms
        auto mass_op = [](double acc, const auto& species) { return acc * species.a(); };
        const double massNumerator = std::accumulate(
            reaction.reactants().begin(),
            reaction.reactants().end(),
            1.0,
            mass_op
        );
        const double massDenominator = std::accumulate(
            reaction.products().begin(),
            reaction.products().end(),
            1.0,
            mass_op
        );

        // Accumulate partition functions
        auto pf_op = [&](double acc, const auto& species) {
            return acc * m_partitionFunction->evaluate(species.z(), species.a(), T9);
        };
        const double partitionFunctionNumerator = std::accumulate(
            reaction.reactants().begin(),
            reaction.reactants().end(),
            1.0,
            pf_op
        );
        const double partitionFunctionDenominator = std::accumulate(
            reaction.products().begin(),
            reaction.products().end(),
            1.0,
            pf_op
        );

        const double CT = std::pow(massNumerator/massDenominator, 1.5) *
            (partitionFunctionNumerator/partitionFunctionDenominator);

        const double reverseRate = forwardRate * symmetryFactor * CT * expFactor;
        if (!std::isfinite(reverseRate)) {
            return 0.0; // If the reverse rate is not finite, return 0.0
        }
        return reverseRate; // Return the calculated reverse rate

    }

    double GraphEngine::calculateReverseRateTwoBodyDerivative(
        const reaction::Reaction &reaction,
        const double T9,
        const double rho,
        const fourdst::composition::Composition& comp,
        const double reverseRate
    ) const {
        assert(reaction.type() == reaction::ReactionType::LOGICAL_REACLIB || reaction.type() == reaction::ReactionType::REACLIB);

        if (!m_useReverseReactions) {
            LOG_TRACE_L3_LIMIT_EVERY_N(std::numeric_limits<int>::max(), m_logger, "Reverse reactions are disabled. Returning 0.0 for reverse rate of reaction '{}'.", reaction.id());
            return 0.0; // If reverse reactions are not used, return 0.0
        }
        const double d_log_kFwd = reaction.calculate_log_rate_partial_deriv_wrt_T9(T9, rho, {}, {}, {});

        auto log_deriv_pf_op = [&](double acc, const auto& species) {
            const double g = m_partitionFunction->evaluate(species.z(), species.a(), T9);
            const double dg_dT = m_partitionFunction->evaluateDerivative(species.z(), species.a(), T9);
            return (g == 0.0) ? acc : acc + (dg_dT / g);
        };

        const double reactant_log_derivative_sum = std::accumulate(
            reaction.reactants().begin(),
            reaction.reactants().end(),
            0.0,
            log_deriv_pf_op
        );

        const double product_log_derivative_sum = std::accumulate(
            reaction.products().begin(),
            reaction.products().end(),
            0.0,
            log_deriv_pf_op
        );

        const double d_log_C = reactant_log_derivative_sum - product_log_derivative_sum;

        const double d_log_exp = reaction.qValue() / (m_constants.kB * T9 * T9);

        const double log_total_derivative = d_log_kFwd + d_log_C + d_log_exp;

        return reverseRate * log_total_derivative; // Return the derivative of the reverse rate with respect to T9

    }

    bool GraphEngine::isUsingReverseReactions(
        scratch::StateBlob& ctx
    ) const {
        return m_useReverseReactions;
    }

    size_t GraphEngine::getSpeciesIndex(
        scratch::StateBlob& ctx,
        const fourdst::atomic::Species &species
    ) const {
        return m_speciesToIndexMap.at(species); // Returns the index of the species in the stoichiometry matrix
    }

    PrimingReport GraphEngine::primeEngine(
        scratch::StateBlob& ctx,
        const NetIn &netIn
    ) const {

        NetIn fullNetIn;
        fourdst::composition::Composition composition;

        for (const auto& sp : m_networkSpecies) {
            composition.registerSpecies(sp);
            if (netIn.composition.contains(sp)) {
                composition.setMolarAbundance(sp, netIn.composition.getMolarAbundance(sp));
            }
        }

        fullNetIn.composition = composition;
        fullNetIn.temperature = netIn.temperature;
        fullNetIn.density = netIn.density;

        std::optional<std::vector<reaction::ReactionType>> reactionTypesToIgnore = std::nullopt;
        if (!m_useReverseReactions) {
            reactionTypesToIgnore = {reaction::ReactionType::WEAK};
        }

        auto primingReport = primeNetwork(ctx, fullNetIn, *this, reactionTypesToIgnore);

        return primingReport;
    }

    fourdst::composition::Composition GraphEngine::collectComposition(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        double T9,
        double rho
    ) const {
        for (const auto &species: comp.getRegisteredSpecies()) {
            if (!m_networkSpeciesMap.contains(species.name())) {
                throw exceptions::BadCollectionError("Cannot collect composition from GraphEngine as " + std::string(species.name()) + " present in input composition does not exist in the network species map");
            }
        }
        fourdst::composition::Composition result;
        for (const auto& species : m_networkSpecies ) {
            result.registerSpecies(species);
            if (comp.contains(species)) {
                double Y = comp.getMolarAbundance(species);
                if (Y < 0.0 && std::abs(Y) <= 1e-16) {
                    result.setMolarAbundance(species, 0.0);
                } else if (Y < 0.0 && std::abs(Y) >= 1e-16) {
                    throw fourdst::composition::exceptions::InvalidCompositionError(std::format("Molar abundance for species {} is negative (Y = {}). GraphEngine does not support non-physical negative abundances, even if they are very small in magnitude (clamp is 1e-16). Check input composition for validity.", species.name(), Y));
                } else {
                    result.setMolarAbundance(species, Y);
                }
            }
        }
        return result;
    }

    SpeciesStatus GraphEngine::getSpeciesStatus(
        scratch::StateBlob& ctx,
        const fourdst::atomic::Species &species
    ) const {
        if (m_networkSpeciesMap.contains(species.name())) {
            return SpeciesStatus::ACTIVE;
        }
        return SpeciesStatus::NOT_PRESENT;

    }

    std::optional<StepDerivatives<double>> GraphEngine::getMostRecentRHSCalculation(
        scratch::StateBlob& ctx
    ) const {
        const auto *state = scratch::get_state<scratch::GraphEngineScratchPad, true>(ctx);
        if (!state->most_recent_rhs_calculation.has_value()) {
            return std::nullopt;
        }
        return state->most_recent_rhs_calculation.value();
    }

    StepDerivatives<double> GraphEngine::calculateAllDerivativesUsingPrecomputation(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const std::vector<double> &bare_rates,
        const std::vector<double> &bare_reverse_rates,
        const double T9,
        const double rho,
        const reaction::ReactionSet &activeReactions,
        const std::span<const double> abundances
    ) const {
        auto *state = scratch::get_state<scratch::GraphEngineScratchPad, true>(ctx);
        // --- Abundances for the mass-action products and screening (raw integrator iterates when given) ---
        if (abundances.empty()) {
            state->local_abundance_cache = networkAbundances(comp);
        } else {
            state->local_abundance_cache.assign(abundances.begin(), abundances.end());
        }

        LOG_TRACE_L3(m_logger, "Computing screening factors for {} active reactions.", activeReactions.size());
        // --- Calculate screening factors ---
        const std::vector<double> screeningFactors = m_screeningModel->calculateScreeningFactors(
            activeReactions,
            m_networkSpecies,
            state->local_abundance_cache,
            T9,
            rho
        );

        StepDerivatives<double> result;
        std::vector<double> dydt_scratch(m_networkSpecies.size(), 0.0);

        const PrecomputationKernelResults kernel = accumulate_flows_serial(
            ctx,
            state->local_abundance_cache,
            screeningFactors,
            bare_rates,
            bare_reverse_rates,
            rho,
            activeReactions
        );
        dydt_scratch = kernel.dydt_vector;
        result.neutrinoEnergyLossRate = kernel.total_neutrino_energy_loss_rate;
        result.totalNeutrinoFlux = kernel.total_neutrino_flux;

        // load scratch into result.dydt
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            result.dydt[m_networkSpecies[i]] = dydt_scratch[i];
        }

        // --- Calculate the nuclear energy generation rate ---
        double massProductionRate = 0.0; // [mol][s^-1]
        for (const auto & species : m_networkSpecies) {
            massProductionRate += result.dydt[species] * species.mass() * m_constants.u;
        }
        result.nuclearEnergyGenerationRate = -massProductionRate * m_constants.Na * m_constants.c * m_constants.c; // [erg][s^-1][g^-1]
        result.nuclearEnergyGenerationRate -= result.neutrinoEnergyLossRate;
        return result;

    }

    screening::ScreeningType GraphEngine::getScreeningModel(
        scratch::StateBlob& ctx
    ) const {
        return m_screeningType;
    }

    bool GraphEngine::isPrecomputationEnabled(
        scratch::StateBlob& ctx
    ) const {
        return m_usePrecomputation;
    }

    const partition::PartitionFunction & GraphEngine::getPartitionFunction(
        scratch::StateBlob& ctx
    ) const {
        return *m_partitionFunction;
    }

    double GraphEngine::calculateMolarReactionFlow(
        scratch::StateBlob& ctx,
        const reaction::Reaction &reaction,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {

        const double Ye = comp.getElectronAbundance();

        const double bareFlow = calculateMolarReactionFlow<double>(
            reaction,
            comp.getMolarAbundanceVector(),
            T9,
            rho,
            Ye,
            0.0,
            [&comp](const fourdst::atomic::Species& species) -> std::optional<size_t> {
                if (comp.contains(species)) { // Species present in the composition
                    return comp.getSpeciesIndex(species);
                }
                return std::nullopt; // Species not present
            }
        );
        if (bareFlow == 0.0 || m_screeningType == screening::ScreeningType::BARE) {
            return bareFlow;
        }

        // Apply the same screening factor the RHS uses, so this is the forward flow the network actually integrates.
        // Callers ask for every reaction at one composition, so the network-ordered abundances are cached per composition.
        auto *state = scratch::get_state<scratch::GraphEngineScratchPad, true>(ctx);
        const std::size_t compHash = comp.hash();
        if (state->flow_abundances_hash != compHash || state->flow_abundances.size() != m_networkSpecies.size()) {
            state->flow_abundances.assign(m_networkSpecies.size(), 0.0);
            for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
                if (comp.contains(m_networkSpecies[i])) {
                    state->flow_abundances[i] = comp.getMolarAbundance(m_networkSpecies[i]);
                }
            }
            state->flow_abundances_hash = compHash;
        }
        return m_screeningModel->calculateScreeningFactor(reaction, m_networkSpecies, state->flow_abundances, T9, rho) * bareFlow;
    }

    NetworkJacobian GraphEngine::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const std::vector<double> Y = networkAbundances(comp);
        if (m_jacobianMethod == JacobianMethod::ANALYTIC) {
            if (std::optional<NetworkJacobian> jac = analyticJacobian(comp, Y, T9, rho)) {
                return std::move(*jac);
            }
            LOG_TRACE_L1_LIMIT_EVERY_N(1000, m_logger, "Screening model has no analytic derivatives; using the CppAD Jacobian.");
        }
        return autoDiffJacobian(ctx, Y, T9, rho);
    }

    std::optional<NetworkJacobian> GraphEngine::generateAnalyticJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const std::span<const double> abundances
    ) const {
        if (abundances.empty()) {
            return analyticJacobian(comp, networkAbundances(comp), T9, rho);
        }
        return analyticJacobian(comp, abundances, T9, rho);
    }

    NetworkJacobian GraphEngine::generateAutoDiffJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const std::span<const double> abundances
    ) const {
        if (abundances.empty()) {
            return autoDiffJacobian(ctx, networkAbundances(comp), T9, rho);
        }
        return autoDiffJacobian(ctx, abundances, T9, rho);
    }

    std::expected<StepDerivatives<double>, EngineStatus> GraphEngine::calculateRHSAndEnergyFromRawAbundances(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const std::span<const double> rawAbundances,
        const double T9,
        const double rho
    ) const {
        if (rawAbundances.size() != m_networkSpecies.size()) {
            throw exceptions::BadRHSEngineError(std::format(
                "calculateRHSAndEnergyFromRawAbundances: expected {} abundances (one per network species), got {}.",
                m_networkSpecies.size(), rawAbundances.size()));
        }
        // With no negative abundance the iterate is exactly comp, so use the ordinary (cached) evaluation
        if (std::ranges::none_of(rawAbundances, [](const double y) { return y < 0.0; }) || !m_usePrecomputation) {
            return calculateRHSAndEnergy(ctx, comp, T9, rho, false);
        }
        auto* state = scratch::get_state<scratch::GraphEngineScratchPad, true>(ctx);
        std::vector<double> bare_rates;
        std::vector<double> bare_reverse_rates;
        compute_bare_rates(comp, T9, rho, m_reactions, bare_rates, bare_reverse_rates);
        StepDerivatives<double> result = calculateAllDerivativesUsingPrecomputation(
            ctx, comp, bare_rates, bare_reverse_rates, T9, rho, m_reactions, rawAbundances
        );
        state->most_recent_rhs_calculation = result;
        return result;
    }

    NetworkJacobian GraphEngine::generateJacobianMatrixFromRawAbundances(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const std::span<const double> rawAbundances,
        const double T9,
        const double rho
    ) const {
        if (rawAbundances.size() != m_networkSpecies.size()) {
            throw std::invalid_argument(std::format(
                "generateJacobianMatrixFromRawAbundances: expected {} abundances (one per network species), got {}.",
                m_networkSpecies.size(), rawAbundances.size()));
        }
        if (m_jacobianMethod == JacobianMethod::ANALYTIC) {
            if (std::optional<NetworkJacobian> jac = analyticJacobian(comp, rawAbundances, T9, rho)) {
                return std::move(*jac);
            }
        }
        return autoDiffJacobian(ctx, rawAbundances, T9, rho);
    }

    std::vector<double> GraphEngine::networkAbundances(const fourdst::composition::CompositionAbstract &comp) const {
        std::vector<double> Y(m_networkSpecies.size(), 0.0);
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            if (comp.contains(m_networkSpecies[i])) {
                Y[i] = comp.getMolarAbundance(m_networkSpecies[i]);
            }
        }
        return Y;
    }

    std::optional<NetworkJacobian> GraphEngine::analyticJacobian(
        const fourdst::composition::CompositionAbstract &comp,
        const std::span<const double> Y,
        const double T9,
        const double rho
    ) const {
        const size_t numSpecies = m_networkSpecies.size();
        const std::vector<double> abundances(Y.begin(), Y.end());

        // --- Screening factors and their composition derivatives, from the same abundances as the RHS ---
        const std::vector<double> screeningFactors = m_screeningModel->calculateScreeningFactors(
            m_reactions, m_networkSpecies, abundances, T9, rho
        );
        const std::optional<std::vector<double>> dlnf_dY = m_screeningModel->calculateScreeningLogDerivatives(
            m_reactions, m_networkSpecies, abundances, T9, rho
        );
        if (!dlnf_dY) {
            return std::nullopt;
        }

        // --- Rates exactly as the RHS evaluates them (from comp: T9, rho and Ye) ---
        std::vector<double> bare_rates;
        std::vector<double> bare_reverse_rates;
        compute_bare_rates(comp, T9, rho, m_reactions, bare_rates, bare_reverse_rates);
        // Electron-capture rates scale with rho * Ye, and Ye = sum_k Z_k Y_k, so d ln k / dY_k = Z_k / Ye
        const double Ye = comp.getElectronAbundance();

        std::vector<double> J(numSpecies * numSpecies, 0.0); // row-major: J[i * numSpecies + j] = dYdot_i / dY_j
        std::vector<double> partials;
        std::vector<std::pair<size_t, double>> flowDerivative; // (column, dR/dY_column) from the mass-action products

        // y^n and n y^(n-1) as the precomputed RHS forms them (powers are small integers)
        const auto integer_power = [](const double y, const int n) {
            double result = 1.0;
            for (int k = 0; k < n; ++k) { result *= y; }
            return result;
        };
        // prod_u Y[idx_u]^p_u and its partial with respect to each factor, built as a product over the other factors
        // so that it stays exact when an abundance is zero
        const auto product_and_partials = [&](const std::vector<size_t>& indices, const std::vector<int>& powers) {
            partials.assign(indices.size(), 0.0);
            double product = 1.0;
            for (size_t u = 0; u < indices.size(); ++u) {
                product *= integer_power(Y[indices[u]], powers[u]);
            }
            for (size_t m = 0; m < indices.size(); ++m) {
                double partial = powers[m] * integer_power(Y[indices[m]], powers[m] - 1);
                for (size_t u = 0; u < indices.size(); ++u) {
                    if (u != m) { partial *= integer_power(Y[indices[u]], powers[u]); }
                }
                partials[m] = partial;
            }
            return product;
        };

        for (size_t r = 0; r < m_reactions.size(); ++r) {
            const reaction::Reaction& reaction = m_reactions[r];
            const PrecomputedReaction& precomp = m_precomputed_reactions[r];
            const double screening = screeningFactors[r];

            // forward flow F = f k sym rho^(n-1) prod_u Y_u^p_u
            const size_t numReactants = reaction.reactants().size();
            const double forwardScale = screening * bare_rates[r] * precomp.symmetry_factor *
                std::pow(rho, numReactants > 1 ? static_cast<double>(numReactants) - 1 : 0.0);
            const double forwardFlow = forwardScale * product_and_partials(precomp.unique_reactant_indices, precomp.reactant_powers);
            flowDerivative.clear();
            for (size_t m = 0; m < precomp.unique_reactant_indices.size(); ++m) {
                flowDerivative.emplace_back(precomp.unique_reactant_indices[m], forwardScale * partials[m]);
            }

            // reverse flow B = f k_rev sym_rev rho^(m-1) prod_u Y_u^q_u, when the RHS includes it
            double reverseFlow = 0.0;
            if (precomp.reverse_symmetry_factor != 0.0 && m_useReverseReactions) {
                const size_t numProducts = reaction.products().size();
                const double reverseScale = screening * bare_reverse_rates[r] * precomp.reverse_symmetry_factor *
                    std::pow(rho, numProducts > 1 ? static_cast<double>(numProducts) - 1 : 0.0);
                reverseFlow = reverseScale * product_and_partials(precomp.unique_product_indices, precomp.product_powers);
                for (size_t m = 0; m < precomp.unique_product_indices.size(); ++m) {
                    flowDerivative.emplace_back(precomp.unique_product_indices[m], -reverseScale * partials[m]);
                }
            }

            // Terms that touch every column: screening, dR/dY_k += (F - B) dln f/dY_k, and for electron captures
            // dF/dY_k += F Z_k / Ye (Ye is that of the clamped iterate, so only non-negative Y_k contribute)
            const double* dlnf_row = dlnf_dY->data() + r * numSpecies;
            const bool electronCapture = precomp.electron_capture && Ye > 0.0;
            const double netFlow = forwardFlow - reverseFlow;

            for (size_t a = 0; a < precomp.affected_species_indices.size(); ++a) {
                const double nu = precomp.stoichiometric_coefficients[a];
                double* J_row = J.data() + precomp.affected_species_indices[a] * numSpecies;
                for (const auto& [col, dR] : flowDerivative) {
                    J_row[col] += nu * dR;
                }
                if (netFlow != 0.0) {
                    for (size_t col = 0; col < numSpecies; ++col) {
                        double dR = netFlow * dlnf_row[col];
                        if (electronCapture && Y[col] >= 0.0) { dR += forwardFlow * m_networkSpecies[col].z() / Ye; }
                        J_row[col] += nu * dR;
                    }
                }
            }
        }
        return packJacobian(J);
    }

    NetworkJacobian GraphEngine::autoDiffJacobian(
        scratch::StateBlob& ctx,
        const std::span<const double> Y,
        const double T9,
        const double rho
    ) const {
        auto *state = scratch::get_state<scratch::GraphEngineScratchPad, true>(ctx);
        LOG_TRACE_L1_LIMIT_EVERY_N(1000, m_logger, "Generating jacobian matrix for T9={}, rho={}..", T9, rho);
        const size_t numSpecies = m_networkSpecies.size();

        // 1. Pack the input variables into a vector for CppAD
        std::vector<double> adInput(numSpecies + 2, 0.0); // +2 for T9 and rho
        std::ranges::copy(Y, adInput.begin());
        adInput[numSpecies]     = T9;  // T9
        adInput[numSpecies + 1] = rho; // rho

        // 2. Calculate the full jacobian (the tape's outputs are dY/dt for every species, then eps)
        assert(state->rhsADFun.has_value() && "RHS ADFun not recorded before Jacobian generation.");
        const std::vector<double> dotY = state->rhsADFun.value().Jacobian(adInput);

        // 3. Keep the species block
        std::vector<double> J(numSpecies * numSpecies);
        for (size_t i = 0; i < numSpecies; ++i) {
            for (size_t j = 0; j < numSpecies; ++j) {
                J[i * numSpecies + j] = dotY[i * (numSpecies + 2) + j];
            }
        }
        return packJacobian(J);
    }

    NetworkJacobian GraphEngine::packJacobian(const std::span<const double> dense) const {
        const size_t numSpecies = m_networkSpecies.size();
        Eigen::SparseMatrix<double> jacobianMatrix(static_cast<long>(numSpecies), static_cast<long>(numSpecies));
        std::vector<Eigen::Triplet<double>> triplets;
        for (size_t i = 0; i < numSpecies; ++i) {
            for (size_t j = 0; j < numSpecies; ++j) {
                const double value = dense[i * numSpecies + j];
                if (std::abs(value) > MIN_JACOBIAN_THRESHOLD || i == j) { // Always keep diagonal elements to avoid pathological stiffness
                    triplets.emplace_back(i, j, value);
                }
            }
        }
        jacobianMatrix.setFromTriplets(triplets.begin(), triplets.end());
        LOG_TRACE_L1_LIMIT_EVERY_N(1000, m_logger, "Jacobian matrix generated with dimensions: {} rows x {} columns.", jacobianMatrix.rows(), jacobianMatrix.cols());

        auto index_to_species = [this](const size_t index) -> fourdst::atomic::Species {
            if (index < m_networkSpecies.size()) {
                return m_networkSpecies[index];
            }
            throw std::out_of_range("Index out of range in index_to_species mapping.");
        };
        return NetworkJacobian(jacobianMatrix, index_to_species);
    }

    NetworkJacobian GraphEngine::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const std::vector<fourdst::atomic::Species> &activeSpecies
    ) const {
        const size_t k_active = activeSpecies.size();

        // --- 1. Get the list of global indices ---
        std::vector<size_t> active_indices;
        active_indices.reserve(k_active);

        for (const auto& species : activeSpecies) {
            assert(involvesSpecies(ctx, species));
            active_indices.push_back(getSpeciesIndex(ctx, species));
        }

        // --- 2. Build the k x k sparsity pattern ---
        SparsityPattern sparsityPattern;
        sparsityPattern.reserve(k_active * k_active);

        for (const size_t i_global : active_indices) { // k rows
            for (const size_t j_global : active_indices) { // k columns
                sparsityPattern.emplace_back(i_global, j_global);
            }
        }

        // --- 3. Call the sparse reverse-mode implementation ---
        return generateJacobianMatrix(ctx, comp, T9, rho, sparsityPattern);
    }

    NetworkJacobian GraphEngine::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const SparsityPattern &sparsityPattern
    ) const {
        auto *state = scratch::get_state<scratch::GraphEngineScratchPad, true>(ctx);
        // --- Compute the intersection of the requested sparsity pattern with the full sparsity pattern ---
        SparsityPattern intersectionSparsityPattern;
        for (const auto& entry : sparsityPattern) {
            if (m_full_sparsity_set.contains(entry)) {
                intersectionSparsityPattern.push_back(entry);
            }
        }

        // --- Pack the input variables into a vector for CppAD ---
        const size_t numSpecies = m_networkSpecies.size();
        std::vector<double> x(numSpecies + 2, 0.0);
        size_t i = 0;
        for (const auto& species: m_networkSpecies) {
            double Yi = 0.0; // Small floor to avoid issues with zero abundances
            if (comp.contains(species)) {
                Yi = std::max(comp.getMolarAbundance(species), 1e-30);
            }
            x[i] = Yi;
            i++;
        }
        x[numSpecies] = T9;
        x[numSpecies + 1] = rho;

        // --- Convert into CppAD Sparsity pattern ---
        const size_t nnz = intersectionSparsityPattern.size(); // Number of non-zero entries in the sparsity pattern
        std::vector<size_t> row_indices(nnz);
        std::vector<size_t> col_indices(nnz);

        for (size_t k = 0; k < nnz; ++k) {
            row_indices[k] = intersectionSparsityPattern[k].first;
            col_indices[k] = intersectionSparsityPattern[k].second;
        }

        std::vector<double> values(nnz);
        const size_t num_rows_jac = numSpecies + 1; // num species + epsilon
        const size_t num_cols_jac = numSpecies + 2; // +2 for T9 and rho

        CppAD::sparse_rc<std::vector<size_t>> CppAD_sparsity_pattern(num_rows_jac, num_cols_jac, nnz);
        std::size_t sparsity_hash = 0;
        for (size_t k = 0; k < nnz; ++k) {
            size_t local_intersection_hash = utils::hash_combine(intersectionSparsityPattern[k].first, intersectionSparsityPattern[k].second);
            sparsity_hash = utils::hash_combine(sparsity_hash, local_intersection_hash);

            CppAD_sparsity_pattern.set(k, intersectionSparsityPattern[k].first, intersectionSparsityPattern[k].second);
        }

        // --- Check cache for existing subset ---
        if (!state->jacobianSubsetCache.contains(sparsity_hash)) {
            state->jacobianSubsetCache.emplace(sparsity_hash, CppAD_sparsity_pattern);
            state->jac_work.clear();
        } else {
            if (state->jacWorkCache.contains(sparsity_hash)) {
                state->jac_work.clear();
                state->jac_work = state->jacWorkCache.at(sparsity_hash);
            }
        }
        auto& jac_subset = state->jacobianSubsetCache.at(sparsity_hash);

        assert(state->rhsADFun.has_value() && "RHS ADFun not recorded before Jacobian generation.");
        state->rhsADFun.value().sparse_jac_rev(
            x,
            jac_subset, // Sparse Jacobian output
            m_full_jacobian_sparsity_pattern,
            "cppad",
            state->jac_work // Work vector for CppAD
        );

        // --- Stash the now populated work vector in the cache if not already present ---
        if (!state->jacWorkCache.contains(sparsity_hash)) {
            state->jacWorkCache.emplace(sparsity_hash, state->jac_work);
        }

        Eigen::SparseMatrix<double> jacobianMatrix(static_cast<long>(numSpecies), static_cast<long>(numSpecies));
        std::vector<Eigen::Triplet<double> > triplets;
        for (size_t k = 0; k < nnz; ++k) {
            const size_t row = jac_subset.row()[k];
            const size_t col = jac_subset.col()[k];
            const double value = jac_subset.val()[k];

            if (std::abs(value) > MIN_JACOBIAN_THRESHOLD || row == col) { // Always keep diagonal elements to avoid pathological stiffness
                triplets.emplace_back(row, col, value);
            }
        }
        jacobianMatrix.setFromTriplets(triplets.begin(), triplets.end());
        auto index_to_species = [this](const size_t index) -> fourdst::atomic::Species {
            if (index < m_networkSpecies.size()) {
                return m_networkSpecies[index];
            }
            throw std::out_of_range("Index out of range in index_to_species mapping.");
        };
        NetworkJacobian jac(jacobianMatrix, index_to_species);
        return jac;
    }

    void GraphEngine::exportToDot(
        scratch::StateBlob& ctx,
        const std::string &filename
    ) const {
        LOG_TRACE_L1(m_logger, "Exporting network graph to DOT file: {}", filename);

        std::ofstream dotFile(filename);
        if (!dotFile.is_open()) {
            LOG_ERROR(m_logger, "Failed to open file for writing: {}", filename);
            m_logger->flush_log();
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }

        dotFile << "digraph NuclearReactionNetwork {\n";
        dotFile << "    graph [rankdir=LR, splines=true, overlap=false, bgcolor=\"#f0f0f0\"];\n";
        dotFile << "    node [shape=circle, style=filled, fillcolor=\"#a7c7e7\", fontname=\"Helvetica\"];\n";
        dotFile << "    edge [fontname=\"Helvetica\", fontsize=10];\n\n";

        // 1. Define all species as nodes
        dotFile << "    // --- Species Nodes ---\n";
        for (const auto& species : m_networkSpecies) {
            dotFile << "    \"" << species.name() << "\" [label=\"" << species.name() << "\"];\n";
        }
        dotFile << "\n";

        // 2. Define all reactions as intermediate nodes and connect them
        dotFile << "    // --- Reaction Edges ---\n";
        for (const auto& reaction : m_reactions) {
            // Create a unique ID for the reaction node
            std::string reactionNodeId = "reaction_" + std::string(reaction->id());

            // Define the reaction node (small, black dot)
            dotFile << "    \"" << reactionNodeId << "\" [shape=point, fillcolor=black, width=0.1, height=0.1, label=\"\"];\n";

            // Draw edges from reactants to the reaction node
            for (const auto& reactant : reaction->reactants()) {
                dotFile << "    \"" << reactant.name() << "\" -> \"" << reactionNodeId << "\";\n";
            }

            // Draw edges from the reaction node to products
            for (const auto& product : reaction->products()) {
                dotFile << "    \"" << reactionNodeId << "\" -> \"" << product.name() << "\" [label=\"" << reaction->qValue() << " MeV\"];\n";
            }
            dotFile << "\n";
        }

        dotFile << "}\n";
        dotFile.close();
        LOG_TRACE_L1(m_logger, "Successfully exported network to {}", filename);
    }

    void GraphEngine::exportToCSV(
        scratch::StateBlob& ctx,
        const std::string &filename
    ) const {
        LOG_TRACE_L1(m_logger, "Exporting network graph to CSV file: {}", filename);

        std::ofstream csvFile(filename, std::ios::out | std::ios::trunc);
        if (!csvFile.is_open()) {
            LOG_ERROR(m_logger, "Failed to open file for writing: {}", filename);
            m_logger->flush_log();
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
        csvFile << "Reaction;Reactants;Products;Q-value;sources;rates\n";
        for (const auto& reaction : m_reactions) {
            // Dynamic cast to REACLIBReaction to access specific properties
            csvFile << reaction->id() << ";";
            // Reactants
            size_t count = 0;
            for (const auto& reactant : reaction->reactants()) {
                csvFile << reactant.name();
                if (++count < reaction->reactants().size()) {
                    csvFile << ",";
                }
            }
            csvFile << ";";
            count = 0;
            for (const auto& product : reaction->products()) {
                csvFile << product.name();
                if (++count < reaction->products().size()) {
                    csvFile << ",";
                }
            }
            csvFile << ";" << reaction->qValue() << ";";
            // Reaction coefficients
            csvFile << "\n";
        }
        csvFile.close();
        LOG_TRACE_L1(m_logger, "Successfully exported network graph to {}", filename);
    }

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> GraphEngine::getSpeciesTimescales(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        return getSpeciesTimescales(ctx, comp, T9, rho, m_reactions);
    }

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> GraphEngine::getSpeciesTimescales(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const reaction::ReactionSet &activeReactions
    ) const {
        const double Ye = comp.getElectronAbundance();

        auto result = calculateAllDerivatives<double>(
            comp.getMolarAbundanceVector(),
            T9,
            rho,
            Ye,
            0.0,
            [&comp](const fourdst::atomic::Species& species) -> std::optional<size_t> {
                if (comp.contains(species)) { // Species present in the composition
                    return comp.getSpeciesIndex(species);
                }
                return std::nullopt; // Species not present
            },
            [&activeReactions](const reaction::Reaction& reaction) -> bool {
                return activeReactions.contains(reaction);
            }
        );
        const std::map<fourdst::atomic::Species, double>& dydt = result.dydt;

        std::unordered_map<fourdst::atomic::Species, double> speciesTimescales;
        speciesTimescales.reserve(m_networkSpecies.size());
        for (const auto& species : m_networkSpecies) {
            double timescale = std::numeric_limits<double>::infinity();
            if (std::abs(dydt.at(species)) > 0.0) {
                timescale = std::abs(comp.getMolarAbundance(species) / dydt.at(species));
            }
            speciesTimescales.emplace(species, timescale);
        }
        return speciesTimescales;
    }

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> GraphEngine::getSpeciesDestructionTimescales(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        return getSpeciesDestructionTimescales(ctx, comp, T9, rho, m_reactions);
    }

    std::expected<std::unordered_map<fourdst::atomic::Species, double>, EngineStatus> GraphEngine::getSpeciesDestructionTimescales(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const reaction::ReactionSet &activeReactions
    ) const {
        const double Ye = comp.getElectronAbundance();
        const std::vector<double>& Y = comp.getMolarAbundanceVector();

        auto speciesLookup = [&comp](const fourdst::atomic::Species& species) -> std::optional<size_t> {
            if (comp.contains(species)) { // Species present in the composition
                return comp.getSpeciesIndex(species);
            }
            return std::nullopt; // Species not present
        };

        std::unordered_map<fourdst::atomic::Species, double> speciesDestructionTimescales;
        speciesDestructionTimescales.reserve(m_networkSpecies.size());
        for (const auto& species : m_networkSpecies) {
            double netDestructionFlow = 0.0;
            for (const auto& reaction : m_reactions) {
                if (reaction->stoichiometry(species) < 0) {
                    const auto flow = calculateMolarReactionFlow<double>(
                        *reaction,
                        Y,
                        T9,
                        rho,
                        Ye,
                        0.0,
                        speciesLookup
                    );
                    netDestructionFlow += flow;
                }
            }
            double timescale = std::numeric_limits<double>::infinity();
            if (netDestructionFlow != 0.0) {
                timescale = comp.getMolarAbundance(species) / netDestructionFlow;
            }
            speciesDestructionTimescales.emplace(species, timescale);
        }
        return speciesDestructionTimescales;
    }

    fourdst::composition::Composition GraphEngine::project(
        scratch::StateBlob& ctx,
        const NetIn &netIn
    ) const {
        // A species the network cannot evolve would otherwise be silently dropped from the solver state (and from the
        // output composition), losing its mass
        for (const auto& species : netIn.composition.getRegisteredSpecies()) {
            if (!m_networkSpeciesMap.contains(species.name()) && netIn.composition.getMolarAbundance(species) > 0.0) {
                throw exceptions::BadCollectionError(std::format(
                    "Input composition contains {} (molar abundance {:.3e}), which is not in the network. Remove it or "
                    "build a network that includes it.", species.name(), netIn.composition.getMolarAbundance(species)));
            }
        }
        fourdst::composition::Composition baseUpdatedComposition = netIn.composition;
        for (const auto& species : m_networkSpecies) {
            if (!netIn.composition.contains(species)) {
                baseUpdatedComposition.registerSpecies(species);
            }
        }
        return baseUpdatedComposition;
    }

    void GraphEngine::recordADTape() {

        LOG_TRACE_L1(m_logger, "Recording AD tape for the RHS calculation...");

        // Task 1: Set dimensions and initialize the matrix
        const size_t numSpecies = m_networkSpecies.size();
        if (numSpecies == 0) {
            LOG_ERROR(m_logger, "Cannot record AD tape: No species in the network.");
            m_logger->flush_log();
            throw std::runtime_error("Cannot record AD tape: No species in the network.");
        }
        const size_t numADInputs = numSpecies + 2; // Y + T9 + rho

        // --- CppAD Tape Recording ---
        // 1. Declare independent variable (adY)
        //    We also initialize the dummy variable for tape recording (these tell CppAD what the derivative chain looks like).
        //    Their numeric values are irrelevant except for in so far as they avoid numerical instabilities.

        // Distribute total mass fraction uniformly between species in the dummy variable space
        const auto uniformMassFraction = static_cast<CppAD::AD<double>>(1.0 / static_cast<double>(numSpecies));
        std::vector<CppAD::AD<double>> adInput(numADInputs, uniformMassFraction);
        adInput[numSpecies]     = 1.0; // Dummy T9
        adInput[numSpecies + 1] = 1.0; // Dummy rho

        // 3. Declare independent variables (what CppAD will differentiate wrt.)
        //    This also beings the tape recording process.
        CppAD::Independent(adInput);

        const std::vector<CppAD::AD<double>> adY(adInput.begin(), adInput.begin() + static_cast<long>(numSpecies));
        const CppAD::AD<double> adT9  = adInput[numSpecies];
        const CppAD::AD<double> adRho = adInput[numSpecies + 1];

        // Dummy values for Ye and mue to let taping happen
        const CppAD::AD<double> adYe = 1e6;
        const CppAD::AD<double> adMue = 10.0;


        // 5. Call the actual templated function
        // We let T9 and rho be constant, so we pass them as fixed values.
        auto result = calculateAllDerivatives<CppAD::AD<double>>(
            adY,
            adT9,
            adRho,
            adYe,
            adMue,
            [&](const fourdst::atomic::Species& querySpecies) -> size_t {
                return m_speciesToIndexMap.at(querySpecies);
            },
            [](const reaction::Reaction& reaction) -> bool {
                return true; // Use all reactions
            }
        );

        // Extract the raw vector from the associative map
        std::vector<CppAD::AD<double>> dependentVector;
        dependentVector.reserve(result.dydt.size() + 1);
        std::ranges::transform(
            result.dydt,
            std::back_inserter(dependentVector),
            [](const auto& kv) {
                return kv.second;
            }
        );
        dependentVector.push_back(result.nuclearEnergyGenerationRate);

        m_authoritativeADFun.Dependent(adInput, dependentVector);
        m_authoritativeADFun.optimize();

        LOG_TRACE_L1(m_logger, "AD tape recorded successfully for the RHS and Eps calculation. Number of independent variables: {}.", adInput.size());
    }

    void GraphEngine::collectAtomicReverseRateAtomicBases(
    ) {
        m_atomicReverseRates.clear();
        m_atomicReverseRates.reserve(m_reactions.size());

        for (const auto& reaction: m_reactions) {
            if (has_reverse_rate(*reaction) and m_useReverseReactions) {
                m_atomicReverseRates.push_back(std::make_unique<AtomicReverseRate>(*reaction, *this));
            } else {
                m_atomicReverseRates.push_back(nullptr);
            }
        }
    }

    void GraphEngine::precomputeNetwork()  {
        LOG_TRACE_L1(m_logger, "Pre-computing constant components of GraphNetwork state...");

        // --- Reverse map for fast species lookups ---
        std::unordered_map<fourdst::atomic::Species, size_t> speciesIndexMap;
        for (size_t i = 0; i < m_networkSpecies.size(); ++i) {
            speciesIndexMap[m_networkSpecies[i]] = i;
        }

        // --- Energy channel classification (once per network) ---
        std::vector<fourdst::atomic::Species> fuels;
        for (const auto& fuel_name : m_config->engine.energyChannels.fuels) {
            const auto it = fourdst::atomic::species.find(fuel_name);
            if (it == fourdst::atomic::species.end()) {
                throw std::invalid_argument(std::format("Unknown energy channel fuel species '{}'.", fuel_name));
            }
            fuels.push_back(it->second);
        }
        m_energy_channels = classify_energy_channels(m_reactions, fuels, m_config->engine.energyChannels.splitHydrogenPPCNO);

        m_precomputed_reactions.clear();
        m_precomputed_reactions.reserve(m_reactions.size());
        m_precomputed_reaction_index_map.clear();
        m_precomputed_reaction_index_map.reserve(m_reactions.size());

        for (size_t i = 0; i < m_reactions.size(); ++i) {
            const auto& reaction = m_reactions[i];
            PrecomputedReaction precomp;
            precomp.reaction_index = i;
            precomp.reaction_type = reaction.type();
            uint64_t reactionHash = reaction.hash(0);

            precomp.reaction_hash = reactionHash;
            m_precomputed_reaction_index_map[reactionHash] = i;

            // --- Precompute forward reaction information ---
            // Count occurrences for each reactant to determine powers and symmetry
            std::unordered_map<size_t, int> reactantCounts;
            for (const auto& reactant: reaction.reactants()) {
                size_t reactantIndex = speciesIndexMap.at(reactant);
                reactantCounts[reactantIndex]++;
            }

            double symmetryDenominator = 1.0;
            for (const auto& [index, count] : reactantCounts) {
                precomp.unique_reactant_indices.push_back(index);
                precomp.reactant_powers.push_back(count);

                symmetryDenominator *= std::tgamma(count + 1);
            }

            precomp.symmetry_factor = 1.0/symmetryDenominator;

            // --- Precompute reverse reaction information ---
            if (has_reverse_rate(reaction)) {
                std::unordered_map<size_t, int> productCounts;
                for (const auto& product : reaction.products()) {
                    productCounts[speciesIndexMap.at(product)]++;
                }
                double reverseSymmetryDenominator = 1.0;
                for (const auto& [index, count] : productCounts) {
                    precomp.unique_product_indices.push_back(index);
                    precomp.product_powers.push_back(count);
                    reverseSymmetryDenominator *= std::tgamma(count + 1);
                }

                precomp.reverse_symmetry_factor = 1.0/reverseSymmetryDenominator;
            } else {
                precomp.unique_product_indices.clear();
                precomp.product_powers.clear();
                precomp.reverse_symmetry_factor = 0.0; // No reverse reaction for weak reactions
            }

            // --- Precompute stoichiometry information ---
            const auto stoichiometryMap = reaction.stoichiometry();
            precomp.affected_species_indices.reserve(stoichiometryMap.size());
            precomp.stoichiometric_coefficients.reserve(stoichiometryMap.size());

            double massChange = 0.0; // [u]
            for (const auto& [species, coeff] : stoichiometryMap) {
                precomp.affected_species_indices.push_back(speciesIndexMap.at(species));
                precomp.stoichiometric_coefficients.push_back(coeff);
                massChange += static_cast<double>(coeff) * species.mass();
            }

            // --- Precompute energy channel attribution (matches the eps definition in calculateAllDerivativesUsingPrecomputation) ---
            // The neutrino energy lookup matches species names against a table, so it is done once here, not per RHS call
            precomp.emits_neutrino = is_reaclib_weak(reaction);
            precomp.neutrino_energy_per_flow = neutrino_energy_loss_per_molar_flow(reaction);
            precomp.energy_per_flow = -massChange * m_constants.u * m_constants.Na * m_constants.c * m_constants.c
                                      - precomp.neutrino_energy_per_flow;
            precomp.channel_index = m_energy_channels.reaction_channel[i];

            // --- Composition dependence of the rate itself (used by the analytic Jacobian) ---
            const auto* reaclib = dynamic_cast<const reaction::ReaclibReaction*>(&reaction);
            precomp.electron_capture = reaclib != nullptr && reaclib->is_electron_capture();
            precomp.analytic_rate_derivative = reaclib != nullptr;

            m_precomputed_reactions.push_back(std::move(precomp));
        }
        LOG_TRACE_L1(m_logger, "Pre-computation complete. Precomputed data for {} reactions.", m_precomputed_reactions.size());
    }

    bool GraphEngine::AtomicReverseRate::forward(
        const size_t p,
        const size_t q,
        const CppAD::vector<bool> &vx,
        CppAD::vector<bool> &vy,
        const CppAD::vector<double> &tx,
        CppAD::vector<double> &ty
    ) {

        if ( p != 0) { return false; }
        const double T9 = tx[0];


        // We can pass a dummy comp and rho because reverse rates should only be calculated for strong reactions whose
        // rates of progression do not depend on composition or density.
        const fourdst::composition::Composition dummyComp;
        const double reverseRate = m_engine.calculateReverseRate(m_reaction, T9, 0.0, dummyComp);
        ty[0] = reverseRate; // Store the reverse rate in the output vector

        if (vx.size() > 0) {
            vy[0] = vx[0];
        }
        return true;
    }

    bool GraphEngine::AtomicReverseRate::reverse(
        const size_t q,
        const CppAD::vector<double> &tx,
        const CppAD::vector<double> &ty,
        CppAD::vector<double> &px,
        const CppAD::vector<double> &py
    ) {
        const double T9 = tx[0];
        const double reverseRate = ty[0];

        const double derivative = m_engine.calculateReverseRateTwoBodyDerivative(m_reaction, T9, 0, {}, reverseRate);

        px[0] = py[0] * derivative; // Return the derivative of the reverse rate with respect to T9

        return true;
    }

    bool GraphEngine::AtomicReverseRate::for_sparse_jac(
        const size_t q,
        const CppAD::vector<std::set<size_t>> &r,
        CppAD::vector<std::set<size_t>> &s
    ) {
        s[0] = r[0];
        return true;
    }

    bool GraphEngine::AtomicReverseRate::rev_sparse_jac(
        const size_t q,
        const CppAD::vector<std::set<size_t>> &rt,
        CppAD::vector<std::set<size_t>> &st
    ) {
        st[0] = rt[0];
        return true;
    }

    bool GraphEngine::AtomicReverseRate::for_sparse_jac(
        const size_t q,
        const CppAD::vector<bool> &r,
        CppAD::vector<bool> &s,
        const CppAD::vector<double> &x
    ) {
        constexpr size_t n = 1;
        constexpr size_t m = 1;

        CPPAD_ASSERT_KNOWN(r.size() == n * q, "AtomicReverseRate::for_sparse_jac: 'r' size is incorrect.");
        CPPAD_ASSERT_KNOWN(s.size() == m * q, "AtomicReverseRate::for_sparse_jac: 's' size is incorrect.");

        // S = R
        for (size_t j = 0; j < q; j++) {
            // s(0,j) = r(0,j)
            s[j*m] = r[j*n];
        }

        return true;
    }

    bool GraphEngine::AtomicReverseRate::rev_sparse_jac(
        const size_t q,
        const CppAD::vector<bool> &rt,
        CppAD::vector<bool> &st,
        const CppAD::vector<double> &x
    ) {
        constexpr size_t n = 1;
        constexpr size_t m = 1;

        CPPAD_ASSERT_KNOWN(rt.size() == n * q, "AtomicReverseRate::for_sparse_jac: 'r' size is incorrect.");
        CPPAD_ASSERT_KNOWN(st.size() == m * q, "AtomicReverseRate::for_sparse_jac: 's' size is incorrect.");

        // st = rt
        for (size_t j = 0; j < q; j++) {
            // st(j, 0) = rt(j, 0)
            st[j * n] = rt[j * m];
        }

        return true;
    }

}
