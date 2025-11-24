#include "gridfire/reaction/weak/weak_rate_library.h"
#include "gridfire/reaction/weak/weak.h"

#include <array>
#include <ranges>
#include <unordered_map>
#include <expected>
#include <vector>
#include <format>

#include "gridfire/reaction/weak/weak_interpolator.h"

#include "xxhash64.h"
#include "fourdst/atomic/species.h"


namespace {
    std::unordered_map<fourdst::atomic::SpeciesErrorType, std::string> SpeciesErrorTypeMap = {
        {fourdst::atomic::SpeciesErrorType::ELEMENT_SYMBOL_NOT_FOUND, "Element symbol not found (Z out of range)"},
        {fourdst::atomic::SpeciesErrorType::SPECIES_SYMBOL_NOT_FOUND, "Species symbol not found ((A,Z) out of range)"}
    };

    std::string normalize_species_id(const fourdst::atomic::Species& species) {
        auto result = std::string(species.name());
        std::ranges::transform(result, result.begin(), ::tolower);
        std::erase(result, '-');
        return result;
    }

    fourdst::atomic::Species resolve_weak_product(
        const gridfire::rates::weak::WeakReactionType type,
        const fourdst::atomic::Species& reactant
    ) {
        using namespace fourdst::atomic;
        using namespace gridfire::rates::weak;

        int zMod = 0;
        switch (type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                zMod = 1;
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                zMod = -1;
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                zMod = 1;
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                zMod = -1;
                break;
        }
        std::expected<Species, SpeciesErrorType> product = az_to_species(reactant.a(), reactant.z() + zMod);

        if (product.has_value()) {
            return product.value();
        }
        const std::string msg = std::format(
            "Failed to resolve weak reaction product (A: {}, Z: {}) for reactant {} (looked up A: {}, Z: {}): {}",
            reactant.a(),
            reactant.z(),
            reactant.name(),
            reactant.a(),
            reactant.z() + zMod,
            SpeciesErrorTypeMap.at(product.error())
        );
        throw std::runtime_error(msg);
    }

    std::string resolve_weak_id(
        const gridfire::rates::weak::WeakReactionType type,
        const fourdst::atomic::Species& reactant,
        const fourdst::atomic::Species& product
    ) {
        using namespace gridfire::rates::weak;

        std::string id;
        std::string reactant_id = normalize_species_id(reactant);
        std::string product_id = normalize_species_id(product);
        switch (type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                id = std::format("{}(,ν|)e-,{}", reactant_id, product_id);
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                id = std::format("{}(,ν)e+,{}", reactant_id, product_id);
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                id = std::format("{}(e-,ν){}", reactant_id, product_id);
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                id = std::format("{}(e+,ν|){}", reactant_id, product_id);
                break;
        }
        return id;
    }
}

namespace gridfire::rates::weak {
    WeakReactionMap::WeakReactionMap() {
        using namespace fourdst::atomic;


        // ReSharper disable once CppUseStructuredBinding
        for (const auto& weak_reaction_record : UNIFIED_WEAK_DATA) {
            std::expected<Species, SpeciesErrorType> species_result = az_to_species(weak_reaction_record.A, weak_reaction_record.Z);
            if (!species_result.has_value()) {
                const SpeciesErrorType type = species_result.error();
                const std::string msg = std::format(
                    "Failed to load weak reaction data for (A={}, Z={}) with error: {}",
                    weak_reaction_record.A, weak_reaction_record.Z, SpeciesErrorTypeMap.at(type)
                );
                throw std::runtime_error(msg);
            }
            const Species& species = species_result.value();

            if (weak_reaction_record.log_beta_minus > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                m_weak_network[species].push_back(
                    WeakReactionEntry{
                        WeakReactionType::BETA_MINUS_DECAY,
                        weak_reaction_record.t9,
                        weak_reaction_record.log_rhoye,
                        weak_reaction_record.mu_e,
                        weak_reaction_record.log_beta_minus,
                        weak_reaction_record.log_antineutrino_loss_bd
                    }
                );
            }
            if (weak_reaction_record.log_beta_plus > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                m_weak_network[species].push_back(
                    WeakReactionEntry{
                        WeakReactionType::BETA_PLUS_DECAY,
                        weak_reaction_record.t9,
                        weak_reaction_record.log_rhoye,
                        weak_reaction_record.mu_e,
                        weak_reaction_record.log_beta_plus,
                        weak_reaction_record.log_neutrino_loss_ec
                    }
                );
            }
            if (weak_reaction_record.log_electron_capture > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                m_weak_network[species].push_back(
                    WeakReactionEntry{
                        WeakReactionType::ELECTRON_CAPTURE,
                        weak_reaction_record.t9,
                        weak_reaction_record.log_rhoye,
                        weak_reaction_record.mu_e,
                        weak_reaction_record.log_electron_capture,
                        weak_reaction_record.log_neutrino_loss_ec
                    }
                );
            }
            if (weak_reaction_record.log_positron_capture > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                m_weak_network[species].push_back(
                    WeakReactionEntry{
                        WeakReactionType::POSITRON_CAPTURE,
                        weak_reaction_record.t9,
                        weak_reaction_record.log_rhoye,
                        weak_reaction_record.mu_e,
                        weak_reaction_record.log_positron_capture,
                        weak_reaction_record.log_antineutrino_loss_bd
                    }
                );
            }
        }
    }

    std::vector<WeakReactionEntry> WeakReactionMap::get_all_reactions() const {
        std::vector<WeakReactionEntry> reactions;
        for (const auto &species_reactions: m_weak_network | std::views::values) {
            reactions.insert(reactions.end(), species_reactions.begin(), species_reactions.end());
        }
        return reactions;
    }

    std::expected<std::vector<WeakReactionEntry>, WeakMapError> WeakReactionMap::get_species_reactions(
        const fourdst::atomic::Species &species) const {
        if (m_weak_network.contains(species)) {
            return m_weak_network.at(species);
        }
        return std::unexpected(WeakMapError::SPECIES_NOT_FOUND);
    }

    std::expected<std::vector<WeakReactionEntry>, WeakMapError> WeakReactionMap::get_species_reactions(
        const std::string &species_name) const {
        const fourdst::atomic::Species& species = fourdst::atomic::species.at(species_name);
        if (m_weak_network.contains(species)) {
            return m_weak_network.at(species);
        }
        return std::unexpected(WeakMapError::SPECIES_NOT_FOUND);
    }

    WeakReaction::WeakReaction(
        const fourdst::atomic::Species &species,
        const WeakReactionType type,
        const WeakRateInterpolator &interpolator
    ) :
    m_reactant(species),
    m_product(resolve_weak_product(type, species)),
    m_reactants({m_reactant}),
    m_products({m_product}),
    m_reactant_a(species.a()),
    m_reactant_z(species.z()),
    m_product_a(m_product.a()),
    m_product_z(m_product.z()),
    m_id(resolve_weak_id(type, species, m_product)),
    m_type(type),
    m_interpolator(interpolator),
    m_atomic(m_interpolator, m_reactant_a, m_reactant_z, m_type) {}

    double WeakReaction::calculate_rate(
        const double T9,
        const double rho,
        const double Ye,
        const double mue,
        const std::vector<double> &Y,
        const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
    ) const {
        return calculate_rate<double>(T9, rho, Ye, mue, Y, index_to_species_map);
    }

    CppAD::AD<double> WeakReaction::calculate_rate(
        const CppAD::AD<double> T9,
        const CppAD::AD<double> rho,
        const CppAD::AD<double> Ye,
        const CppAD::AD<double> mue,
        const std::vector<CppAD::AD<double>> &Y,
        const std::unordered_map<size_t,fourdst::atomic::Species>& index_to_species_map
    ) const {
        return calculate_rate<CppAD::AD<double>>(T9, rho, Ye, mue, Y, index_to_species_map);
    }

    std::string_view WeakReaction::id() const {
        return m_id;
    }

    const std::vector<fourdst::atomic::Species> & WeakReaction::reactants() const {
        return m_reactants;
    }

    const std::vector<fourdst::atomic::Species> & WeakReaction::products() const {
        return m_products;
    }

    bool WeakReaction::contains(const fourdst::atomic::Species &species) const {
        return contains_reactant(species) || contains_product(species);
    }

    bool WeakReaction::contains_reactant(const fourdst::atomic::Species& species) const {
        if (m_reactant == species) {
            return true;
        }
        return false;
    }

    bool WeakReaction::contains_product(const fourdst::atomic::Species &species) const {
        if (m_product == species) {
            return true;
        }
        return false;
    }

    std::unordered_set<fourdst::atomic::Species> WeakReaction::all_species() const {
        return {m_reactant, m_product};
    }

    std::unordered_set<fourdst::atomic::Species> WeakReaction::reactant_species() const {
        return {m_reactant};
    }

    std::unordered_set<fourdst::atomic::Species> WeakReaction::product_species() const {
        return {m_product};
    }

    size_t WeakReaction::num_species() const {
        return 2;
    }

    int WeakReaction::stoichiometry(const fourdst::atomic::Species &species) const {
        if (species == m_reactant) {
            return -1;
        }
        if (species == m_product) {
            return 1;
        }
        return 0;
    }

    std::unordered_map<fourdst::atomic::Species, int> WeakReaction::stoichiometry() const {
        return {
            {m_reactant, -1},
            {m_product, 1}
        };
    }

    uint64_t WeakReaction::hash(const uint64_t seed) const {
        const std::string reaction_string = std::format(
            "{}:{}({})",
            m_reactant.name(),
            m_product.name(),
            static_cast<int>(m_type)
        );
        return XXHash64::hash(reaction_string.data(), reaction_string.size(), seed);
    }

    double WeakReaction::qValue() const {
        // We ignore neutrino mass as it is negligible compared to other masses here.
        double Q_MeV = 0.0;

        const double parentMass_u = m_reactant.mass();
        const double daughterMass_u = m_product.mass();
        const double electronMass_MeV = m_constants.electronMassMeV;

        const double nuclearMassDiff_MeV = (parentMass_u - daughterMass_u) * m_constants.u_to_MeV;
        switch (m_type) {
            case WeakReactionType::BETA_PLUS_DECAY:
                Q_MeV = nuclearMassDiff_MeV - 2.0 * electronMass_MeV;
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                Q_MeV = nuclearMassDiff_MeV + 2.0 * electronMass_MeV;
                break;
            case WeakReactionType::BETA_MINUS_DECAY: // Same as electron capture so we can simply fall through
            case WeakReactionType::ELECTRON_CAPTURE:
                Q_MeV = nuclearMassDiff_MeV;
                break;
        }
        return Q_MeV;
    }

    double WeakReaction::calculate_energy_generation_rate(
        const double T9,
        const double rho,
        const double Ye,
        const double mue,
        const std::vector<double> &Y,
        const std::unordered_map<size_t, fourdst::atomic::Species> &index_to_species_map
    ) const {
        std::expected<WeakRatePayload, InterpolationError> rates =  m_interpolator.get_rates(
            static_cast<uint16_t>(m_reactant_a),
            static_cast<uint8_t>(m_reactant_z),
            T9,
           std::log10(rho * Ye)
        );

        if (!rates.has_value()) {
            const InterpolationErrorType type = rates.error().type;
            const std::string msg = std::format(
                "Failed to interpolate weak rate for (A={}, Z={}) at T9={}, log10(rho*Ye)={}, mu_e={} with error: {}",
                m_reactant.name(), m_reactant_a, m_reactant_z, T9, std::log10(rho * Ye), mue, InterpolationErrorTypeMap.at(type)
            );
            throw std::runtime_error(msg);
        }

        double logNeutrinoLossRate = 0.0;
        if (m_type == WeakReactionType::BETA_PLUS_DECAY || m_type == WeakReactionType::ELECTRON_CAPTURE) {
            logNeutrinoLossRate = rates->log_antineutrino_loss_bd;
        } else if (m_type == WeakReactionType::BETA_MINUS_DECAY || m_type == WeakReactionType::POSITRON_CAPTURE) {
            logNeutrinoLossRate = rates->log_neutrino_loss_ec;
        }

        const double neutrinoLossRate = std::pow(10, logNeutrinoLossRate);

        const double EDeposited_MeV = qValue() - neutrinoLossRate;

        // We reimplement this logic here instead of calling calculate_rate() to avoid
        // doing the interpolation twice (since the payload has already been interpolated).
        const double logRate = get_log_rate_from_payload(rates.value());
        double lambda = 0.0;
        if (logRate > GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
            lambda = std::pow(10, logRate);
        }
        return lambda * EDeposited_MeV; // returns in MeV / s

    }

    CppAD::AD<double> WeakReaction::calculate_energy_generation_rate(
        const CppAD::AD<double> &T9,
        const CppAD::AD<double> &rho,
        const CppAD::AD<double> &Ye,
        const CppAD::AD<double> &mue,
        const std::vector<CppAD::AD<double>> &Y,
        const std::unordered_map<size_t, fourdst::atomic::Species> &index_to_species_map
    ) const {
        const CppAD::AD<double> log_rhoYe = CppAD::log10(rho * Ye);

        const std::vector<CppAD::AD<double>> ax = {T9, log_rhoYe};
        std::vector<CppAD::AD<double>> ay(2); // 2 outputs are the reaction rate (1/s) and the neutrino loss (MeV)

        m_atomic(ax, ay); // Note: We needed to make m_atomic mutable to allow this call in a const method.

        const CppAD::AD<double> rateConstant = ay[0];
        const CppAD::AD<double> NuLoss = ay[1];

        return rateConstant * (qValue() - NuLoss); // returns in MeV / s
    }

    double WeakReaction::calculate_log_rate_partial_deriv_wrt_T9(
        const double T9,
        const double rho,
        const double Ye,
        const double mue,
        const fourdst::composition::Composition &composition
    ) const {
        const double log_rhoYe = std::log10(rho * Ye);
        std::expected<WeakRateDerivatives, InterpolationError> rates =  m_interpolator.get_rate_derivatives(
            static_cast<uint16_t>(m_reactant_a),
            static_cast<uint8_t>(m_reactant_z),
            T9,
            log_rhoYe
        );
        if (!rates.has_value()) {
            const InterpolationErrorType type = rates.error().type;
            const std::string msg = std::format(
                "Failed to interpolate weak rate for (A={}, Z={}) at T9={}, log10(rho*Ye)={}, mu_e={} with error: {}",
                m_reactant.name(), m_reactant_a, m_reactant_z, T9, log_rhoYe, mue, InterpolationErrorTypeMap.at(type)
            );
            throw std::runtime_error(msg);
        }

        double logRate = 0.0;
        switch (m_type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                logRate = rates->d_log_beta_minus[0];
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                logRate = rates->d_log_beta_plus[0];
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                logRate = rates->d_log_electron_capture[0];
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                logRate = rates->d_log_positron_capture[0];
                break;
        }
        return logRate;
    }

    reaction::ReactionType WeakReaction::type() const {
        return reaction::ReactionType::WEAK;
    }

    std::unique_ptr<reaction::Reaction> WeakReaction::clone() const {
        std::unique_ptr<reaction::Reaction> reaction_ptr = std::make_unique<WeakReaction>(
            m_reactant,
            m_type,
            m_interpolator
        );
        return reaction_ptr;
    }

    bool WeakReaction::is_reverse() const {
        return false;
    }

    const WeakRateInterpolator & WeakReaction::getWeakRateInterpolator() const {
        return m_interpolator;
    }

    double WeakReaction::get_log_rate_from_payload(const WeakRatePayload &payload) const {
        double logRate = 0.0;
        switch (m_type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                logRate = payload.log_beta_minus;
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                logRate = payload.log_beta_plus;
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                logRate = payload.log_electron_capture;
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                logRate = payload.log_positron_capture;
                break;
        }
        return logRate;
    }

    double WeakReaction::get_log_neutrino_loss_from_payload(const WeakRatePayload &payload) const {
        double logNeutrinoLoss = 0.0;
        switch (m_type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                logNeutrinoLoss = payload.log_antineutrino_loss_bd;
                break;
            case WeakReactionType::BETA_PLUS_DECAY: // Same as electron capture so we can simply fall through
            case WeakReactionType::ELECTRON_CAPTURE:
                logNeutrinoLoss = payload.log_neutrino_loss_ec;
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                logNeutrinoLoss = payload.log_antineutrino_loss_bd;
                break;
        }
        return logNeutrinoLoss;
    }

    size_t WeakReaction::countReactantOccurrences(const fourdst::atomic::Species &species) const {
        size_t count = 0;
        for (const auto& reactant : m_reactants) {
            if (reactant == species) {
                count++;
            }
        }
        return count;
    }

    size_t WeakReaction::countProductOccurrences(const fourdst::atomic::Species &species) const {
        size_t count = 0;
        for (const auto& product : m_products) {
            if (product == species) {
                count++;
            }
        }
        return count;
    }

    // Note that the input vector tx is of size 2: [T9, log10(rho*Ye)]
    bool WeakReaction::AtomicWeakRate::forward (
        const size_t p,
        const size_t q,
        const CppAD::vector<bool> &vx,
        CppAD::vector<bool> &vy,
        const CppAD::vector<double> &tx,
        CppAD::vector<double> &ty
    ) {
        // Doing this explicitly (only allowing p == 0) makes forward mode AD impossible for now.
        if (p != 0) {
            return false;
        }
        const double T9 = tx[0];
        const double log10_rhoye = tx[1];

        const std::expected<WeakRatePayload, InterpolationError> result = m_interpolator.get_rates(
            static_cast<uint16_t>(m_a),
            static_cast<uint8_t>(m_z),
            T9,
            log10_rhoye
        );
        if (!result.has_value()) {
            const InterpolationErrorType type = result.error().type;
            const std::string msg = std::format(
                "Failed to interpolate weak rate for (A={}, Z={}) at T9={}, log10(ρ Ye)={}, with error: {}",
                m_a, m_z, T9, log10_rhoye, InterpolationErrorTypeMap.at(type)
            );
            throw std::runtime_error(msg);
        }
        switch (m_type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                ty[0] = std::pow(10, result.value().log_beta_minus);
                ty[1] = std::pow(10, result.value().log_antineutrino_loss_bd);
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                ty[0] = std::pow(10, result.value().log_beta_plus);
                ty[1] = std::pow(10, result.value().log_neutrino_loss_ec);
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                ty[0] = std::pow(10, result.value().log_electron_capture);
                ty[1] = std::pow(10, result.value().log_neutrino_loss_ec);
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                ty[0] = std::pow(10, result.value().log_positron_capture);
                ty[1] = std::pow(10, result.value().log_antineutrino_loss_bd);
                break;
        }

        if (vx.size() > 0) { // Set up the sparsity pattern. This is saying that all input variables affect the output variable.
            const bool any_input_varies = vx[0] || vx[1];
            vy[0] = any_input_varies;
            vy[1] = any_input_varies;
        }
        return true;
    }

    bool WeakReaction::AtomicWeakRate::reverse(
        const size_t q,
        const CppAD::vector<double> &tx,
        const CppAD::vector<double> &ty,
        CppAD::vector<double> &px,
        const CppAD::vector<double> &py
    ) {
        const double T9 = tx[0];
        const double log10_rhoye = tx[1];

        const double forwardPassRate = ty[0]; // This is the rate from the forward pass.
        const double forwardPassNeutrinoLossRate = ty[1]; // This is the neutrino loss rate from the forward pass.

        const std::expected<WeakRateDerivatives, InterpolationError> result = m_interpolator.get_rate_derivatives(
            static_cast<uint16_t>(m_a),
            static_cast<uint8_t>(m_z),
            T9,
            log10_rhoye
        );

        if (!result.has_value()) {
            const InterpolationErrorType type = result.error().type;
            const std::string msg = std::format(
                "Failed to interpolate weak rate derivatives for (A={}, Z={}) at T9={}, log10(ρ Ye)={}, with error: {}",
                m_a, m_z, T9, log10_rhoye, InterpolationErrorTypeMap.at(type)
            );
            throw std::runtime_error(msg);
        }

        // ReSharper disable once CppUseStructuredBinding
        const WeakRateDerivatives derivatives = result.value();

        std::array<double, 2> dLogRate{}; // d(rate)/dT9, d(rate)/dlogRhoYe
        std::array<double, 2> dLogNuLoss{}; // d(nu loss)/dT9, d(nu loss)/dlogRhoYe
        switch (m_type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                dLogRate[0]     = derivatives.d_log_beta_minus[0];
                dLogRate[1]     = derivatives.d_log_beta_minus[1];
                dLogNuLoss[0]   = derivatives.d_log_antineutrino_loss_bd[0];
                dLogNuLoss[1]   = derivatives.d_log_antineutrino_loss_bd[1];
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                dLogRate[0]     = derivatives.d_log_beta_plus[0];
                dLogRate[1]     = derivatives.d_log_beta_plus[1];
                dLogNuLoss[0]   = derivatives.d_log_neutrino_loss_ec[0];
                dLogNuLoss[1]   = derivatives.d_log_neutrino_loss_ec[1];
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                dLogRate[0]     = derivatives.d_log_electron_capture[0];
                dLogRate[1]     = derivatives.d_log_electron_capture[1];
                dLogNuLoss[0]   = derivatives.d_log_neutrino_loss_ec[0];
                dLogNuLoss[1]   = derivatives.d_log_neutrino_loss_ec[1];
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                dLogRate[0]     = derivatives.d_log_positron_capture[0];
                dLogRate[1]     = derivatives.d_log_positron_capture[1];
                dLogNuLoss[0]   = derivatives.d_log_antineutrino_loss_bd[0];
                dLogNuLoss[1]   = derivatives.d_log_antineutrino_loss_bd[1];
                break;
        }

        const double ln10 = std::log(10.0);

        // Contributions from the reaction rate (output 0)
        px[0] = py[0] * forwardPassRate * ln10 * dLogRate[0];
        px[1] = py[0] * forwardPassRate * ln10 * dLogRate[1];

        // Contributions from the neutrino loss rate (output 1)
        px[0] += py[1] * forwardPassNeutrinoLossRate * ln10 * dLogNuLoss[0];
        px[1] += py[1] * forwardPassNeutrinoLossRate * ln10 * dLogNuLoss[1];

        return true;

    }

    bool WeakReaction::AtomicWeakRate::for_sparse_jac(
        const size_t q,
        const CppAD::vector<std::set<size_t> > &r,
        CppAD::vector<std::set<size_t> > &s
    ) {
        std::set<size_t> all_input_deps;
        all_input_deps.insert(r[0].begin(), r[0].end());
        all_input_deps.insert(r[1].begin(), r[1].end());

        // What this is saying is that both output variables depend on all input variables.
        s[0] = all_input_deps;
        s[1] = all_input_deps;

        return true;
    }

    bool WeakReaction::AtomicWeakRate::rev_sparse_jac(
        const size_t q,
        const CppAD::vector<std::set<size_t> > &rt,
        CppAD::vector<std::set<size_t> > &st
    ) {
        // What this is saying is that all input variables may affect both output variables.
        std::set<size_t> all_output_deps;
        all_output_deps.insert(rt[0].begin(), rt[0].end());
        all_output_deps.insert(rt[1].begin(), rt[1].end());

        st[0] = all_output_deps;
        st[1] = all_output_deps;
        return true;

    }

    bool WeakReaction::AtomicWeakRate::for_sparse_jac(
        const size_t q,
        const CppAD::vector<bool> &r,
        CppAD::vector<bool> &s,
        const CppAD::vector<double> &x
    ) {
        constexpr size_t n = 2; // Number of inputs
        constexpr size_t m = 2; // Number of outputs

        CPPAD_ASSERT_KNOWN(r.size() == q * n, "AtomicWeakRate::for_sparse_jac: 'r' size is incorrect!");
        CPPAD_ASSERT_KNOWN(s.size() == q * m, "AtomicWeakRate::for_sparse_jac: 's' size is incorrect!");

        // Both outputs depend on both inputs
        //  s[i + j*m] represents s(i,j) - output i, direction j
        //  r[k + j*n] represents r(k,j) - input k, direction j

        for (size_t j = 0; j < q; j++) {
            // s(0,j) = r(0,j) || r(1,j) --- output 0 depends on both inputs
            s[0 + j*m] = r[0 + j*n] || r[1 + j*n];

            // s(1,j) = r(0,j) || r(1,j) --- output 1 depends on both inputs
            s[1 + j*m] = r[0 + j*n] || r[1 + j*n];
        }

        return true;
    }

    bool WeakReaction::AtomicWeakRate::rev_sparse_jac(
        const size_t q,
        const CppAD::vector<bool> &rt,
        CppAD::vector<bool> &st,
        const CppAD::vector<double> &x
    ) {
        constexpr size_t n = 2; // Number of inputs
        constexpr size_t m = 2; // Number of outputs

        CPPAD_ASSERT_KNOWN(rt.size() == q * m, "AtomicWeakRate::rev_sparse_jac: 'rt' size is incorrect!");
        CPPAD_ASSERT_KNOWN(st.size() == q * n, "AtomicWeakRate::rev_sparse_jac: 'st' size is incorrect!");

        for (size_t j = 0; j < q; j++) {
            //st(0,j) = rt(0,j) || rt(1,j) --- input 0 affects both outputs
            st[0 + j*n] = rt[0 + j*m] || rt[1 + j*m];

            //st(1,j) = rt(0,j) || rt(1,j) --- input 1 affects both outputs
            st[1 + j*n] = rt[0 + j*m] || rt[1 + j*m];
        }

        return true;
    }
}
