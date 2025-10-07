#pragma once

#define GRIDFIRE_WEAK_REACTION_LIB_SENTINEL -60.0

#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/weak/weak_types.h"
#include "gridfire/reaction/weak/weak_interpolator.h"

#include "gridfire/engine/engine_abstract.h"

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/constants/const.h"

#include "cppad/cppad.hpp"

#include <memory>
#include <unordered_map>
#include <expected>
#include <vector>
#include <string>
#include <string_view>
#include <unordered_set>



namespace gridfire::rates::weak {
    class WeakReactionMap {
    public:
        WeakReactionMap();
        ~WeakReactionMap() = default;

        [[nodiscard]] std::vector<WeakReactionEntry> get_all_reactions() const;

        [[nodiscard]] std::expected<std::vector<WeakReactionEntry>, WeakMapError> get_species_reactions(
            const fourdst::atomic::Species &species) const;

        [[nodiscard]] std::expected<std::vector<WeakReactionEntry>, WeakMapError> get_species_reactions(
            const std::string &species_name) const;
    private:
        std::unordered_map<fourdst::atomic::Species, std::vector<WeakReactionEntry>> m_weak_network;
    };

    class WeakReaction final : public reaction::Reaction {
    public:
        explicit WeakReaction(
            const fourdst::atomic::Species &species,
            WeakReactionType type,
            const WeakRateInterpolator& interpolator
        );
        [[nodiscard]] double calculate_rate(
            double T9,
            double rho,
            double Ye,
            double mue,
            const std::vector<double> &Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const override;
        [[nodiscard]] CppAD::AD<double> calculate_rate(
            CppAD::AD<double> T9,
            CppAD::AD<double> rho,
            CppAD::AD<double> Ye,
            CppAD::AD<double> mue,
            const std::vector<CppAD::AD<double>> &Y,
            const std::unordered_map<size_t,fourdst::atomic::Species>& index_to_species_map
        ) const override;
        [[nodiscard]] std::string_view id() const override { return m_id; }
        [[nodiscard]] const std::vector<fourdst::atomic::Species> &reactants() const override { return {m_reactant}; }
        [[nodiscard]] const std::vector<fourdst::atomic::Species> &products() const override { return {m_product}; }
        [[nodiscard]] bool contains(const fourdst::atomic::Species &species) const override;
        [[nodiscard]] bool contains_reactant(const fourdst::atomic::Species &species) const override;
        [[nodiscard]] bool contains_product(const fourdst::atomic::Species &species) const override;
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> all_species() const override;
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> reactant_species() const override;
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> product_species() const override;
        [[nodiscard]] size_t num_species() const override { return 2; } // Always 2 for weak reactions
        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, int> stoichiometry() const override;
        [[nodiscard]] int stoichiometry(const fourdst::atomic::Species &species) const override;
        [[nodiscard]] uint64_t hash(uint64_t seed) const override;
        [[nodiscard]] double qValue() const override;
        [[nodiscard]] double calculate_energy_generation_rate(
            double T9,
            double rho,
            double Ye,
            double mue,
            const std::vector<double>& Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const override;
        [[nodiscard]] CppAD::AD<double> calculate_energy_generation_rate(
            const CppAD::AD<double> &T9,
            const CppAD::AD<double> &rho,
            const CppAD::AD<double> &Ye,
            const CppAD::AD<double> &mue,
            const std::vector<CppAD::AD<double>> &Y,
            const std::unordered_map<size_t, fourdst::atomic::Species> &index_to_species_map
        ) const override;
        [[nodiscard]] double calculate_forward_rate_log_derivative(
            double T9,
            double rho,
            double Ye,
            double mue,
            const fourdst::composition::Composition& composition
        ) const override;
        [[nodiscard]] reaction::ReactionType type() const override { return reaction::ReactionType::WEAK; }
        [[nodiscard]] std::unique_ptr<Reaction> clone() const override;
        [[nodiscard]] bool is_reverse() const override { return false; };

    private:
        template<IsArithmeticOrAD T>
        T calculate_rate(
            T T9,
            T rho,
            T Ye,
            T mue,
            const std::vector<T> &Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const;

        double get_log_rate_from_payload(const WeakRatePayload& payload) const;

    private:
        class AtomicWeakRate final : public CppAD::atomic_base<double> {
        public:
            AtomicWeakRate(
                const WeakRateInterpolator& interpolator,
                const size_t a,
                const size_t z,
                const WeakReactionType type
            ) :
            CppAD::atomic_base<double>(std::format("atomic-{}-{}-weak-rate", a, z)),
            m_interpolator(interpolator),
            m_a(a),
            m_z(z) ,
            m_type(type) {}

            bool forward(
                size_t p,
                size_t q,
                const CppAD::vector<bool>& vx,
                CppAD::vector<bool>& vy,
                const CppAD::vector<double>& tx,
                CppAD::vector<double>& ty
            ) override;
            bool reverse(
                size_t q,
                const CppAD::vector<double>& tx,
                const CppAD::vector<double>& ty,
                CppAD::vector<double>& px,
                const CppAD::vector<double>& py
            ) override;
            bool for_sparse_jac(
                size_t q,
                const CppAD::vector<std::set<size_t>>&r,
                CppAD::vector<std::set<size_t>>& s
            ) override;
            bool rev_sparse_jac(
                size_t q,
                const CppAD::vector<std::set<size_t>>&rt,
                CppAD::vector<std::set<size_t>>& st
            ) override;

        private:
            const WeakRateInterpolator& m_interpolator;
            const size_t m_a;
            const size_t m_z;
            const WeakReactionType m_type;
        };

        struct constants {
            const fourdst::constant::Constants& c = fourdst::constant::Constants::getInstance();
            fourdst::constant::Constant neutronMassG = c.get("mN");
            fourdst::constant::Constant protonMassG = c.get("mP");
            fourdst::constant::Constant electronMassG = c.get("mE");
            fourdst::constant::Constant speedOfLight = c.get("c");
            fourdst::constant::Constant eVgRelation = c.get("eV_kg"); // note that despite the symbol this is in g NOT kg
            fourdst::constant::Constant MeV2Erg = c.get("MeV_to_erg");
            fourdst::constant::Constant amu = c.get("u");

            double MeVgRelation = eVgRelation.value * 1.0e6;
            double MeVPerGraph =  1.0/MeVgRelation;

            double neutronMassMeV = neutronMassG.value * MeVgRelation;
            double protonMassMeV = protonMassG.value * MeVgRelation;
            double electronMassMeV = electronMassG.value * MeVgRelation;

            double u_to_MeV = (amu.value * speedOfLight.value * speedOfLight.value)/MeV2Erg.value;
        };

    private:
        const constants m_constants;
        fourdst::atomic::Species m_reactant;
        fourdst::atomic::Species m_product;

        size_t m_reactant_a;
        size_t m_reactant_z;
        size_t m_product_a;
        size_t m_product_z;

        std::string m_id;
        WeakReactionType m_type;

        const WeakRateInterpolator& m_interpolator;

        AtomicWeakRate m_atomic;

    };

    template<IsArithmeticOrAD T>
    T WeakReaction::calculate_rate(
        T T9,
        T rho,
        T Ye,
        T mue,
        const std::vector<T> &Y,
        const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
    ) const {
        T log_rhoYe = CppAD::log10(rho * Ye);

        T rateConstant = static_cast<T>(0.0);
        if constexpr (std::is_same_v<T, CppAD::AD<double>>) { // Case where T is an AD type
            std::vector<T> ax = {T9, log_rhoYe, mue};
            std::vector<T> ay(1);
            m_atomic(ax, ay);
            rateConstant = static_cast<T>(ay[0]);
        } else { // The case where T is of type double
            const std::expected<WeakRatePayload, InterpolationError> result = m_interpolator.get_rates(
                static_cast<uint16_t>(m_reactant_a),
                static_cast<uint8_t>(m_reactant_z),
                T9,
                log_rhoYe,
                mue
            );

            if (!result.has_value()) {
                const InterpolationErrorType type = result.error().type;
                const std::string msg = std::format(
                    "Failed to interpolate weak rate for (A={}, Z={}) at T9={}, log10(rho*Ye)={}, mu_e={} with error: {}",
                    m_reactant.name(), m_reactant_a, m_reactant_z, T9, log_rhoYe, mue, InterpolationErrorTypeMap.at(type)
                );
                throw std::runtime_error(msg);
            }
            const WeakRatePayload payload = result.value();
            const double logRate = get_log_rate_from_payload(payload);
            if (logRate <= GRIDFIRE_WEAK_REACTION_LIB_SENTINEL) {
                rateConstant = static_cast<T>(0.0);
            } else {
                rateConstant = CppAD::pow(10.0, logRate);
            }

        }
        return rateConstant;
    }
}



