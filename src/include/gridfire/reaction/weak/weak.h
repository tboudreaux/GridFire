#pragma once

#define GRIDFIRE_WEAK_REACTION_LIB_SENTINEL (-60.0)

#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/weak/weak_types.h"
#include "gridfire/reaction/weak/weak_interpolator.h"

#include "gridfire/engine/engine_abstract.h"

#include "fourdst/atomic/atomicSpecies.h"
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
    /**
     * @class WeakReactionMap
     * @brief Index of available weak reactions keyed by species.
     *
     * Builds an in-memory map from the compiled weak-rate tables and provides
     * simple query helpers to retrieve all weak reactions or those that involve
     * a particular nuclide.
     *
     * Implementation summary: the constructor iterates over UNIFIED_WEAK_DATA and
     * inserts entries keyed by the parent Species. For each channel (β−, β+, e−-capture,
     * e+-capture), if the tabulated log10(rate) is above the sentinel (-60), a
     * WeakReactionEntry is pushed containing the grids t9, log10(rho*Ye), mu_e, the log10(rate),
     * and the corresponding log10(neutrino loss) column.
     */
    class WeakReactionMap {
    public:
        /**
         * @brief Construct the map by loading all weak reaction entries.
         * @post All valid reactions from the compiled data are available via
         *       get_all_reactions() and get_species_reactions().
         * Implementation: iterates UNIFIED_WEAK_DATA, filters any log(rate) <= -60,
         * and groups entries by parent Species.
         */
        WeakReactionMap();
        ~WeakReactionMap() = default;

        /**
         * @brief Return a flat list of all weak reaction entries.
         * @return Vector of WeakReactionEntry records.
         * @par Example
         * @code
         * WeakReactionMap map;
         * auto all = map.get_all_reactions();
         * // iterate or group as needed
         * @endcode
         */
        [[nodiscard]] std::vector<WeakReactionEntry> get_all_reactions() const;

        /**
         * @brief Get all weak reaction entries for a given species.
         * @param species Nuclide to query (A,Z).
         * @return expected<vector<WeakReactionEntry>, WeakMapError>
         *         containing reactions on success or SPECIES_NOT_FOUND on failure.
         * @par Example
         * @code
         * using fourdst::atomic::Species;
         * WeakReactionMap map;
         * Species fe52 = fourdst::atomic::az_to_species(52, 26);
         * if (auto res = map.get_species_reactions(fe52); res) {
         *     for (const auto& e : *res) { } // use e
         * } else {
         *     // handle WeakMapError::SPECIES_NOT_FOUND
         * }
         * @endcode
         */
        [[nodiscard]] std::expected<std::vector<WeakReactionEntry>, WeakMapError> get_species_reactions(
            const fourdst::atomic::Species &species
        ) const;

        /**
         * @brief Get all weak reaction entries for a given species by name.
         * @param species_name Symbolic name (e.g., "Fe52").
         * @return expected<vector<WeakReactionEntry>, WeakMapError>
         *         containing reactions on success or SPECIES_NOT_FOUND on failure.
         * @par Example
         * @code
         * WeakReactionMap map;
         * if (auto res = map.get_species_reactions("Fe52"); res) {
         *     // use *res
         * }
         * @endcode
         */
        [[nodiscard]] std::expected<std::vector<WeakReactionEntry>, WeakMapError> get_species_reactions(
            const std::string &species_name
        ) const;
    private:
        std::unordered_map<fourdst::atomic::Species, std::vector<WeakReactionEntry>> m_weak_network;
    };

    /**
     * @class WeakReaction
     * @brief Concrete Reaction representing a single weak process (beta±, e−/e+ capture).
     *
     * Wraps interpolation logic for tabulated weak rates and provides both scalar and AD
     * interfaces for rate and energy generation. The reactants/products are the parent/daughter
     * nuclei of the weak process.
     *
     * @details the product nucleus is resolved from (A,Z) and channel via
     * simple charge-changing rules (β−: Z+1; β+: Z−1; e− capture: Z−1; e+ capture: Z+1).
     * The reaction ID is formatted like "Parent(channel)Product" with ν/ν̄ decorations, and
     * an internal CppAD atomic (AtomicWeakRate) is prepared for AD energy calculations.
     */
    class WeakReaction final : public reaction::Reaction {
    public:
        /**
         * @brief Construct a WeakReaction for a specific weak channel and parent species.
         * @param species Parent nuclide undergoing the weak process.
         * @param type The weak reaction channel (beta−, beta+, e− capture, e+ capture).
         * @param interpolator Reference to a WeakRateInterpolator providing tabulated data.
         * @pre The product nuclide must be resolvable for the given (species, type).
         * @post Object is ready to compute rates using the provided interpolator.
         * @throws std::runtime_error If the product species cannot be resolved for the channel
         *         (product resolution uses the charge-changing rules described above).
         */
        explicit WeakReaction(
            const fourdst::atomic::Species &species,
            WeakReactionType type,
            const WeakRateInterpolator& interpolator
        );
        /**
         * @brief Scalar weak reaction rate λ(T9, rho, Ye, μe) in 1/s.
         *
         * @details Performs a single interpolation of the weak-rate tables at (T9, log10(rho*Ye), μe).
         * If the selected log10(rate) is ≤ sentinel (-60), returns 0; otherwise returns 10^{log10(rate)}.
         * On interpolation failure, throws with a message including (A,Z) and the state point.
         *
         * @param T9 Temperature in GK (1e9 K).
         * @param rho Mass density (g cm^-3).
         * @param Ye Electron fraction.
         * @param mue Electron chemical potential (MeV).
         * @param Y Composition vector (unused for weak channels).
         * @param index_to_species_map Index-to-species map (unused for weak channels).
         * @return Reaction rate (1/s).
         * @throws std::runtime_error On interpolation failure.
         * @par Example
         * @code
         * double lambda = rxn.calculate_rate(2.0, 1e8, 0.4, 1.5, {}, {});
         * @endcode
         */
        [[nodiscard]] double calculate_rate(
            double T9,
            double rho,
            double Ye,
            double mue,
            const std::vector<double> &Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const override;
        /**
         * @brief AD-enabled weak reaction rate λ(T9, rho, Ye, μe) in 1/s.
         *
         * @details Current implementation returns 0.0. AD support is provided for the energy-generation
         * overload below using an internal CppAD atomic that evaluates both the rate and neutrino
         * loss consistently. A future implementation may mirror that atomic here and return the AD rate.
         *
         * @param T9 Temperature in GK (AD type).
         * @param rho Mass density (g cm^-3, AD type).
         * @param Ye Electron fraction (AD type).
         * @param mue Electron chemical potential (MeV, AD type).
         * @param Y Composition vector (unused for weak channels).
         * @param index_to_species_map Index-to-species map (unused for weak channels).
         * @return Reaction rate (1/s) as CppAD::AD<double> (currently 0.0).
         * @par Example
         * @code
         * using AD = CppAD::AD<double>;
         * AD lambda_ad = rxn.calculate_rate(AD(3.0), AD(1e7), AD(0.5), AD(2.0), {}, {});
         * @endcode
         */
        [[nodiscard]] CppAD::AD<double> calculate_rate(
            CppAD::AD<double> T9,
            CppAD::AD<double> rho,
            CppAD::AD<double> Ye,
            CppAD::AD<double> mue,
            const std::vector<CppAD::AD<double>> &Y,
            const std::unordered_map<size_t,fourdst::atomic::Species>& index_to_species_map
        ) const override;

        /**
         * @brief Unique identifier string for the weak channel.
         * @return A stable string view (e.g., "Fe52(e-,ν)Mn52").
         */
        [[nodiscard]] std::string_view id() const override;

        /**
         * @brief Reactants list (parent nuclide only).
         * @return Vector with the parent species.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species> &reactants() const override;

        /**
         * @brief Products list (daughter nuclide only).
         * @return Vector with the daughter species.
         */
        [[nodiscard]] const std::vector<fourdst::atomic::Species> &products() const override;

        /**
         * @brief Check if a species participates in this weak reaction.
         */
        [[nodiscard]] bool contains(const fourdst::atomic::Species &species) const override;

        /**
         * @brief Check if a species is the reactant (parent).
         */
        [[nodiscard]] bool contains_reactant(const fourdst::atomic::Species &species) const override;

        /**
         * @brief Check if a species is the product (daughter).
         */
        [[nodiscard]] bool contains_product(const fourdst::atomic::Species &species) const override;

        /**
         * @brief Set of both parent and daughter species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> all_species() const override;

        /**
         * @brief Singleton set containing only the parent species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> reactant_species() const override;

        /**
         * @brief Singleton set containing only the daughter species.
         */
        [[nodiscard]] std::unordered_set<fourdst::atomic::Species> product_species() const override;

        /**
         * @brief Number of unique species involved (always 2 for weak reactions).
         */
        [[nodiscard]] size_t num_species() const override;

        /**
         * @brief Full stoichiometry map: parent -1, daughter +1.
         */
        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, int> stoichiometry() const override;

        /**
         * @brief Stoichiometric coefficient for a species: -1 (parent), +1 (daughter), 0 otherwise.
         */
        [[nodiscard]] int stoichiometry(const fourdst::atomic::Species &species) const override;

        /**
         * @brief Content-based 64-bit hash for this reaction.
         */
        [[nodiscard]] uint64_t hash(uint64_t seed) const override;

        /**
         * @brief Q-value (MeV) based on nuclear mass differences and channel.
         *
         * Computes Q = (M_parent − M_daughter) c^2 converted to MeV. For β+ decay subtract 2 m_e c^2,
         * for e+ capture add 2 m_e c^2; for β− and e− capture it is just the nuclear mass difference.
         * Neutrino rest mass is ignored.
         */
        [[nodiscard]] double qValue() const override;

        /**
         * @brief Net energy generation rate (MeV/s) for this weak reaction.
         *
         * Interpolates once to obtain both the log10(rate) and the appropriate log10(neutrino-loss)
         * for the channel, converts to linear values, computes E_deposited = Q − ν_loss, and returns
         * λ · E_deposited. Throws on interpolation failure.
         *
         * Channel mapping for neutrino-loss column:
         *  - β− decay and e+ capture: use log_antineutrino_loss_bd
         *  - β+ decay and e− capture: use log_neutrino_loss_ec
         *
         * @param T9 Temperature in GK.
         * @param rho Density in g cm^-3.
         * @param Ye Electron fraction.
         * @param mue Electron chemical potential (MeV).
         * @param Y Composition vector (unused for weak channels).
         * @param index_to_species_map Index-to-species map (unused for weak channels).
         * @return Energy generation rate in MeV/s.
         * @throws std::runtime_error On interpolation failure.
         * @par Example
         * @code
         * double eps = rxn.calculate_energy_generation_rate(3.0, 1e7, 0.5, 2.0, {}, {});
         * @endcode
         */
        [[nodiscard]] double calculate_energy_generation_rate(
            double T9,
            double rho,
            double Ye,
            double mue,
            const std::vector<double>& Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const override;

        /**
         * @brief AD-enabled net energy generation rate (MeV/s).
         *
         * Uses an internal CppAD atomic to compute two outputs at once: the rate λ and the neutrino
         * loss ν_loss at (T9, log10(rho*Ye), μe). Returns λ · (Q − ν_loss). The atomic throws on
         * interpolation failure.
         *
         * @throws std::runtime_error If the atomic rate evaluation fails to interpolate.
         */
        [[nodiscard]] CppAD::AD<double> calculate_energy_generation_rate(
            const CppAD::AD<double> &T9,
            const CppAD::AD<double> &rho,
            const CppAD::AD<double> &Ye,
            const CppAD::AD<double> &mue,
            const std::vector<CppAD::AD<double>> &Y,
            const std::unordered_map<size_t, fourdst::atomic::Species> &index_to_species_map
        ) const override;

        /**
         * @brief Logarithmic temperature sensitivity of the rate at the given state.
         *
         * Implementation status: requests derivative tables from the interpolator and throws on
         * failure; otherwise the function is not yet implemented and does not return a value.
         * Avoid calling until implemented.
         *
         * @param T9 Temperature in GK.
         * @param rho Density in g cm^-3.
         * @param Ye Electron fraction.
         * @param mue Electron chemical potential (MeV).
         * @param composition Composition context (not used by weak channels presently).
         * @return d ln λ / d ln T9.
         * @throws std::runtime_error On interpolation failure.
         */
        [[nodiscard]] double calculate_log_rate_partial_deriv_wrt_T9(
            double T9,
            double rho,
            double Ye,
            double mue,
            const fourdst::composition::Composition& composition
        ) const override;

        /**
         * @brief Reaction type tag for runtime dispatch.
         */
        [[nodiscard]] reaction::ReactionType type() const override;

        /**
         * @brief Polymorphic deep copy.
         */
        [[nodiscard]] std::unique_ptr<Reaction> clone() const override;

        /**
         * @brief Weak reactions are parameterized in the forward sense (never reverse).
         */
        [[nodiscard]] bool is_reverse() const override;

        /**
         * @brief Access the underlying rate interpolator used by this reaction.
         */
        [[nodiscard]] const WeakRateInterpolator& getWeakRateInterpolator() const;

        [[nodiscard]] std::optional<std::vector<reaction::RateCoefficientSet>> getRateCoefficients() const override { return std::nullopt; }

    private:
        /**
         * @brief Internal unified implementation for scalar/AD rate evaluation.
         * @tparam T double or CppAD::AD<double>.
         * @param T9, rho, Ye, mue Thermodynamic state.
         * @param Y Composition vector (unused for weak channels).
         * @param index_to_species_map Index-to-species map (unused for weak channels).
         * @return Reaction rate (1/s) as T. For double, performs table interpolation and returns
         *         0 when the tabulated log10(rate) ≤ sentinel; for AD, calls the atomic and returns
         *         the first output.
         * @pre T9 > 0, rho > 0, 0 < Ye <= 1.
         * @post No persistent state is modified.
         * @throws std::runtime_error If interpolation fails (double path) or the atomic fails (AD path).
         */
        template<IsArithmeticOrAD T>
        T calculate_rate(
            T T9,
            T rho,
            T Ye,
            T mue,
            const std::vector<T> &Y,
            const std::unordered_map<size_t, fourdst::atomic::Species>& index_to_species_map
        ) const;

        /**
         * @brief Extract the channel-specific log10(rate) from an interpolated payload.
         * Mapping: β−→log_beta_minus, β+→log_beta_plus, e− capture→log_electron_capture,
         * e+ capture→log_positron_capture.
         */
        double get_log_rate_from_payload(const WeakRatePayload& payload) const;

        /**
         * @brief Extract the channel-specific log10(neutrino loss) from a payload.
         * Mapping: β−/e+ capture→log_antineutrino_loss_bd; β+/e− capture→log_neutrino_loss_ec.
         */
        double get_log_neutrino_loss_from_payload(const WeakRatePayload& payload) const;

    public:
        [[nodiscard]] size_t countReactantOccurrences(const fourdst::atomic::Species &species) const override;

        [[nodiscard]] size_t countProductOccurrences(const fourdst::atomic::Species &species) const override;

    private:
        /**
         * @brief CppAD atomic that wraps weak-rate interpolation for AD evaluation.
         *
         * Forward pass computes two outputs (λ, ν_loss) by interpolating the tables at the
         * provided state; reverse pass uses derivative tables to backpropagate adjoints for
         * all three inputs (T9, log10(rho*Ye), μe). Sparsity routines declare full dependence
         * of both outputs on all inputs.
         */
        class AtomicWeakRate final : public CppAD::atomic_base<double> {
        public:
            /**
             * @brief Construct the atomic operation for a specific (A,Z) and channel.
             * @param interpolator Rate source.
             * @param a Mass number A of the parent.
             * @param z Proton number Z of the parent.
             * @param type Weak channel.
             */
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

            /**
             * @brief Forward pass: compute rate and neutrino-loss values for AD.
             * On failure to interpolate, throws a std::runtime_error with details; sets output
             * sparsity such that both outputs depend on all inputs when any input is variable.
             */
            bool forward(
                size_t p,
                size_t q,
                const CppAD::vector<bool>& vx,
                CppAD::vector<bool>& vy,
                const CppAD::vector<double>& tx,
                CppAD::vector<double>& ty
            ) override;

            /**
             * @brief Reverse pass: propagate adjoints using tabulated derivatives.
             * Uses d log10 columns, converting to linear-scale derivatives via ln(10) scaling and
             * chain rule with the forward-pass outputs.
             */
            bool reverse(
                size_t q,
                const CppAD::vector<double>& tx,
                const CppAD::vector<double>& ty,
                CppAD::vector<double>& px,
                const CppAD::vector<double>& py
            ) override;

            /**
             * @brief Forward-mode sparsity for Jacobian.
             */
            bool for_sparse_jac(
                size_t q,
                const CppAD::vector<std::set<size_t>>&r,
                CppAD::vector<std::set<size_t>>& s
            ) override;

            /**
             * @brief Reverse-mode sparsity for Jacobian.
             */
            bool rev_sparse_jac(
                size_t q,
                const CppAD::vector<std::set<size_t>>&rt,
                CppAD::vector<std::set<size_t>>& st
            ) override;

            bool for_sparse_jac(
                size_t q,
                const CppAD::vector<bool> &r,
                CppAD::vector<bool> &s,
                const CppAD::vector<double> &x
            ) override;

            bool rev_sparse_jac(
                size_t q,
                const CppAD::vector<bool> &rt,
                CppAD::vector<bool> &st,
                const CppAD::vector<double> &x
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

        std::vector<fourdst::atomic::Species> m_reactants;
        std::vector<fourdst::atomic::Species> m_products;

        size_t m_reactant_a;
        size_t m_reactant_z;
        size_t m_product_a;
        size_t m_product_z;

        std::string m_id;
        WeakReactionType m_type;

        const WeakRateInterpolator& m_interpolator;

        mutable AtomicWeakRate m_atomic;

    };

    // template implementation lives at the end of the header for AD instantiation
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
            std::vector<T> ax = {T9, log_rhoYe};
            std::vector<T> ay(2);
            m_atomic(ax, ay);
            rateConstant = static_cast<T>(ay[0]);
        } else { // The case where T is of type double
            std::expected<WeakRatePayload, InterpolationError> result = m_interpolator.get_rates(
                static_cast<uint16_t>(m_reactant_a),
                static_cast<uint8_t>(m_reactant_z),
                T9,
                log_rhoYe
            );

            // TODO: Clean this up. When a bit of code needs this many comments to make it clear it is bad code
            if (!result.has_value()) {
                bool okayToClamp = true;
                const auto&[errorType, boundsErrorInfo] = result.error();

                // The logic here is
                //   1. If there is any bounds error in T9 then we do not allow clamping. T9 should be a large enough grid
                //      that the user should not be asking for values outside the grid.
                //   2. If there is no bounds error in T9, but there is a bounds error in log_rhoYe, then we only allow
                //      clamping if the query value is below the minimum of the grid. If it is above the maximum
                //      of the grid, then we do not allow clamping. The reason for this is that at high density,
                //      screening and other effects can make a significant difference to the rates, and
                //      the user should be aware that they are asking for a value outside the grid.

                // There are a couple of safety asserts in here that are only active in debug builds. These are to
                // ensure that our assumptions about the error information are correct. These should really never
                // be triggered, but if they are, they will help us to identify any issues.
                if (errorType == InterpolationErrorType::BOUNDS_ERROR) {
                    assert(boundsErrorInfo.has_value()); // must be true if type is BOUNDS_ERROR, removed in release builds

                    if (boundsErrorInfo->contains(TableAxes::T9)) {
                        okayToClamp = false;
                    } else {
                        assert(boundsErrorInfo->contains(TableAxes::LOG_RHOYE)); // must be true if T9 is not, removed in release builds
                        const BoundsErrorInfo& boundsError = boundsErrorInfo->at(TableAxes::LOG_RHOYE);

                        if (boundsError.queryValue > boundsError.axisMaxValue) {
                            okayToClamp = false;
                        }

                        assert(boundsError.queryValue < boundsError.axisMinValue); // Given the above logic, this must be true, removed in release builds
                    }
                }

                if (!okayToClamp) {
                    const InterpolationErrorType type = result.error().type;
                    const std::string msg = std::format(
                        "Failed to interpolate weak rate for {} (A={}, Z={}) at T9={}, log10(rho*Ye)={}, with error: {}. Clamping disallowed due to either query value being out of bounds in T9 or being above the maximum in log10(rho*Ye).",
                        m_reactant.name(), m_reactant_a, m_reactant_z, T9, log_rhoYe, InterpolationErrorTypeMap.at(type)
                    );
                    throw std::runtime_error(msg);
                }

                // In the case we get here the error was a bounds error in log_rhoYe and the query value was below the minimum of the grid
                //   so the solution is to clamp the query value to the minimum of the grid and try again.
                result = m_interpolator.get_rates(
                    static_cast<uint16_t>(m_reactant_a),
                    static_cast<uint8_t>(m_reactant_z),
                    T9,
                    boundsErrorInfo->at(TableAxes::LOG_RHOYE).axisMinValue
                );

                // Check the result again. If it fails this time then we have a real problem and we throw.
                if (!result.has_value()) {
                    const InterpolationErrorType type = result.error().type;
                    const std::string msg = std::format(
                        "After clamping, failed to interpolate weak rate for {} (A={}, Z={}) at T9={}, log10(rho*Ye)={}, with error: {}",
                        m_reactant.name(), m_reactant_a, m_reactant_z, T9, log_rhoYe, InterpolationErrorTypeMap.at(type)
                    );
                    throw std::runtime_error(msg);
                }
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

