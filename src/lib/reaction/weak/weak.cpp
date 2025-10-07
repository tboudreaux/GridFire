#include "gridfire/reaction/weak/weak_rate_library.h"
#include "gridfire/reaction/weak/weak.h"

#include "fourdst/composition/species.h"

#include <array>
#include <ranges>
#include <unordered_map>
#include <expected>
#include <vector>
#include <format>

#include "gridfire/reaction/weak/weak_interpolator.h"

#include "xxhash64.h"


namespace {
    fourdst::atomic::Species resolve_weak_product(
        const gridfire::rates::weak::WeakReactionType type,
        const fourdst::atomic::Species& reactant
    ) {
        using namespace fourdst::atomic;
        using namespace gridfire::rates::weak;

        std::optional<Species> product; // Use optional so that we can start in a valid "null" state
        switch (type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                product = az_to_species(reactant.a(), reactant.z() + 1);
                return product.value();
            case WeakReactionType::BETA_PLUS_DECAY:
                product = az_to_species(reactant.a(), reactant.z() - 1);
                return product.value();
            case WeakReactionType::ELECTRON_CAPTURE:
                product = az_to_species(reactant.a(), reactant.z() - 1);
                return product.value();
            case WeakReactionType::POSITRON_CAPTURE:
                product = az_to_species(reactant.a(), reactant.z() + 1);
                break;
        }
        if (!product.has_value()) {
            throw std::runtime_error("Failed to resolve weak reaction product for reactant: " + std::string(reactant.name()));
        }
        return product.value();
    }

    std::string resolve_weak_id(
        const gridfire::rates::weak::WeakReactionType type,
        const fourdst::atomic::Species& reactant,
        const fourdst::atomic::Species& product
    ) {
        using namespace gridfire::rates::weak;

        std::string id;
        switch (type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                id = std::format("{}(,ν|)e-,{}", reactant.name(), product.name());
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                id = std::format("{}(,ν)e+,{}", reactant.name(), product.name());
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                id = std::format("{}(e-,ν){}", reactant.name(), product.name());
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                id = std::format("{}(e+,ν|){}", reactant.name(), product.name());
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
            Species species = az_to_species(weak_reaction_record.A, weak_reaction_record.Z);

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
        const fourdst::atomic::Species species = fourdst::atomic::species.at(species_name);
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
        CppAD::AD<double> T9,
        CppAD::AD<double> rho,
        CppAD::AD<double> Ye,
        CppAD::AD<double> mue, const std::vector<CppAD::AD<double>> &Y, const std::unordered_map<size_t,fourdst::atomic::Species>& index_to_species_map
    ) const {
        return static_cast<CppAD::AD<double>>(0.0);
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
            case WeakReactionType::BETA_MINUS_DECAY:
                Q_MeV = nuclearMassDiff_MeV;
                break;
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
            std::log10(rho * Ye),
            mue
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
        std::vector<CppAD::AD<double>> ax = {T9, log_rhoYe, mue};
        std::vector<CppAD::AD<double>> ay(1);
        m_atomic(ax, ay); // TODO: Sort out why this isn't working and checkline 222 in weak.h where a similar line is
        //TODO: think about how to get out neutrino loss in a autodiff safe way. This may mean I need to add an extra output to the atomic base
        //      so that I can get out both the rate and the neutrino loss rate. This will also mean that the sparsity pattern will need to
        //      be updated to account for the extra output.
        CppAD::AD<double> rateConstant = ay[0];

    }

    std::unique_ptr<reaction::Reaction> WeakReaction::clone() const {
        std::unique_ptr<reaction::Reaction> reaction_ptr = std::make_unique<WeakReaction>(
            m_reactant,
            m_type,
            m_interpolator
        );
        return reaction_ptr;
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
        const double mu_e = tx[2];

        const std::expected<WeakRatePayload, InterpolationError> result = m_interpolator.get_rates(
            static_cast<uint16_t>(m_a),
            static_cast<uint8_t>(m_z),
            T9,
            log10_rhoye,
            mu_e
        );
        if (!result.has_value()) {
            const InterpolationErrorType type = result.error().type;
            std::string msg = std::format(
                "Failed to interpolate weak rate for (A={}, Z={}) at T9={}, log10(rho*Ye)={}, mu_e={} with error: {}",
                m_a, m_z, T9, log10_rhoye, mu_e, InterpolationErrorTypeMap.at(type)
            );
        }
        switch (m_type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                ty[0] = std::pow(10, result.value().log_beta_minus);
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                ty[0] = std::pow(10, result.value().log_beta_plus);
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                ty[0] = std::pow(10, result.value().log_electron_capture);
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                ty[0] = std::pow(10, result.value().log_positron_capture);
                break;
        }

        if (vx.size() > 0) {
            vy[0] = vx[0] || vx[1] || vx[2]; // Sets the output sparsity pattern
        }
        return true;
    }

    bool WeakReaction::AtomicWeakRate::reverse(
        size_t q,
        const CppAD::vector<double> &tx,
        const CppAD::vector<double> &ty,
        CppAD::vector<double> &px,
        const CppAD::vector<double> &py
    ) {
        const double T9 = tx[0];
        const double log10_rhoye = tx[1];
        const double mu_e = tx[2];

        const std::expected<WeakRateDerivatives, InterpolationError> result = m_interpolator.get_rate_derivatives(
            static_cast<uint16_t>(m_a),
            static_cast<uint8_t>(m_z),
            T9,
            log10_rhoye,
            mu_e
        );

        if (!result.has_value()) {
            const InterpolationErrorType type = result.error().type;
            const std::string msg = std::format(
                "Failed to interpolate weak rate derivatives for (A={}, Z={}) at T9={}, log10(rho*Ye)={}, mu_e={} with error: {}",
                m_a, m_z, T9, log10_rhoye, mu_e, InterpolationErrorTypeMap.at(type)
            );
            throw std::runtime_error(msg);
        }

        WeakRateDerivatives derivatives = result.value();

        double dT9 = 0.0;
        double dRho = 0.0;
        double dMuE = 0.0;
        switch (m_type) {
            case WeakReactionType::BETA_MINUS_DECAY:
                dT9 = py[0] * derivatives.d_log_beta_minus[0];
                dRho = py[0] * derivatives.d_log_beta_minus[1];
                dMuE = py[0] * derivatives.d_log_beta_minus[2];
                break;
            case WeakReactionType::BETA_PLUS_DECAY:
                dT9 = py[0] * derivatives.d_log_beta_plus[0];
                dRho = py[0] * derivatives.d_log_beta_plus[1];
                dMuE = py[0] * derivatives.d_log_beta_plus[2];
                break;
            case WeakReactionType::ELECTRON_CAPTURE:
                dT9 = py[0] * derivatives.d_log_electron_capture[0];
                dRho = py[0] * derivatives.d_log_electron_capture[1];
                dMuE = py[0] * derivatives.d_log_electron_capture[2];
                break;
            case WeakReactionType::POSITRON_CAPTURE:
                dT9 = py[0] * derivatives.d_log_positron_capture[0];
                dRho = py[0] * derivatives.d_log_positron_capture[1];
                dMuE = py[0] * derivatives.d_log_positron_capture[2];
                break;
        }

        px[0] = py[0] * dT9;       // d(rate)/dT9
        px[1] = py[0] * dRho;      // d(rate)/dlogRhoYe
        px[2] = py[0] * dMuE;      // d(rate)/dMuE

        return true;

    }

    bool WeakReaction::AtomicWeakRate::for_sparse_jac(
        size_t q,
        const CppAD::vector<std::set<size_t> > &r,
        CppAD::vector<std::set<size_t> > &s
    ) {
        s[0] = r[0];
        s[0].insert(r[1].begin(), r[1].end());
        s[0].insert(r[2].begin(), r[2].end());

        return true;
    }

    bool WeakReaction::AtomicWeakRate::rev_sparse_jac(
        size_t q,
        const CppAD::vector<std::set<size_t> > &rt,
        CppAD::vector<std::set<size_t> > &st
    ) {
        // What this is saying is that each of the three input variables (T9, rho, Ye)
        // all only affect the output variable (the rate) since there is only
        // one output variable.
        st[0] = rt[0];
        st[1] = rt[0];
        st[2] = rt[0];
        return true;

    }








}
