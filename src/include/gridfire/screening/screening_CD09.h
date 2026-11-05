#pragma once

#include "gridfire/screening/screening_abstract.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/types/types.h"

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"
#include "fourdst/constants/const.h"

#include "cppad/cppad.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <mutex>
#include <numbers>
#include <optional>
#include <shared_mutex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace gridfire::screening {

    /**
     * @class CD09ScreeningModel
     * @brief Global (weak \f$\to\f$ intermediate \f$\to\f$ strong) plasma screening following Chugunov & DeWitt (2009).
     *
     * Implements the interpolation of Chugunov & DeWitt 2009, PRC 80, 014611, eq. (A4). The rate of the reaction
     * \f$i + j\f$ is enhanced by \f$F_{ij} = e^{h_{ij}}\f$ with
     * \f[
     *     h_{ij} = P_{ij} \left[ f_0\!\left(\frac{\Gamma_i}{t_{ij}}\right) + f_0\!\left(\frac{\Gamma_j}{t_{ij}}\right)
     *              - f_0\!\left(\frac{\Gamma_{\rm comp}}{t_{ij}}\right) \right],
     *     \qquad P_{ij} = \frac{C_{ij} + \Gamma_{ij}^2}{1 + \Gamma_{ij}^2},
     * \f]
     * where
     *   - \f$\Gamma_e = e^2 / (a_e k_B T)\f$, \f$a_e = (3 / 4\pi n_e)^{1/3}\f$;
     *   - \f$\Gamma_i = \Gamma_e Z_i^{5/3}\f$, \f$\Gamma_{\rm comp} = \Gamma_e (Z_i + Z_j)^{5/3}\f$;
     *   - \f$\Gamma_{ij} = Z_i Z_j e^2 / (a_{ij} k_B T)\f$ with \f$a_{ij} = (a_i + a_j)/2\f$, \f$a_i = a_e Z_i^{1/3}\f$;
     *   - \f$\tau_{ij} = \left(27\pi^2 \mu_{ij} Z_i^2 Z_j^2 e^4 / 2\hbar^2 k_B T\right)^{1/3}\f$ (CD09 eq. 10;
     *     \f$= 3 E_{\rm Gamow} / k_B T\f$) and \f$\zeta_{ij} = 3 \Gamma_{ij} / \tau_{ij}\f$ (eq. 12);
     *   - \f$t_{ij} = \left(1 + c_1 \zeta_{ij} + c_2 \zeta_{ij}^2 + c_3 \zeta_{ij}^3\right)^{1/3}\f$ (eq. 25), with
     *     \f$c_1 = 0.013\, y_{ij}^2\f$, \f$c_2 = 0.406\, y_{ij}^{0.14}\f$, \f$c_3 = 0.062\, y_{ij}^{0.19} + 1.8/\Gamma_{ij}\f$
     *     and \f$y_{ij} = 4 Z_i Z_j / (Z_i + Z_j)^2\f$;
     *   - \f$f_0\f$ is the Potekhin & Chabrier (2000, PRE 62, 8554) one-component-plasma free energy per ion
     *     (CD09 eq. 24),
     *     \f[
     *         f_0(\Gamma) = A_1\left[\sqrt{\Gamma(A_2 + \Gamma)} - A_2 \ln\left(\sqrt{\Gamma/A_2} + \sqrt{1 + \Gamma/A_2}\right)\right]
     *                     + 2A_3\left[\sqrt{\Gamma} - \arctan\sqrt{\Gamma}\right]
     *                     + B_1\left[\Gamma - B_2 \ln\left(1 + \Gamma/B_2\right)\right]
     *                     + \frac{B_3}{2} \ln\left(1 + \Gamma^2/B_4\right),
     *     \f]
     *     \f$A_1 = -0.907\f$, \f$A_2 = 0.62954\f$, \f$A_3 = -\sqrt{3}/2 - A_1/\sqrt{A_2}\f$, \f$B_1 = 0.00456\f$,
     *     \f$B_2 = 211.6\f$, \f$B_3 = -10^{-4}\f$, \f$B_4 = 0.00462\f$. It tends to the Debye-Hückel limit
     *     \f$-\Gamma^{3/2}/\sqrt{3}\f$ for \f$\Gamma \ll 1\f$ (\f$A_3\f$ is chosen for that) and to the Madelung energy
     *     for \f$\Gamma \gg 1\f$.
     *
     * The prefactor \f$P_{ij}\f$ is as published in eq. (A4), so \f$P_{ij} \to 1\f$ in the strong limit. Physical
     * constants come from libconstants (CODATA 2022), with \f$e^2 = \alpha \hbar c\f$ in Gaussian units.
     *
     * @par Weak-limit matching and electron degeneracy
     * \f$C_{ij} = h_{ij}^{\rm DH} / h_{ij}^{\rm lin}\f$ fixes the weak-coupling limit, with
     * \f$h_{ij}^{\rm lin} = \Gamma_e^{3/2} \left[(Z_i + Z_j)^{5/2} - Z_i^{5/2} - Z_j^{5/2}\right] / \sqrt{3}\f$ (eq. A2).
     * CD09 use the ion-only Debye-Hückel result (eq. A1), appropriate for the strongly degenerate electrons they
     * consider. Here \f$h^{\rm DH}\f$ is the Salpeter (1954) weak-screening exponent, including the electrons'
     * contribution weighted by their degeneracy:
     * \f[
     *     h_{ij}^{\rm DH} = 0.188\, Z_i Z_j \sqrt{\rho \zeta / T_6^3}, \qquad
     *     \zeta = \sum_k Z_k^2 Y_k + \theta_e \sum_k Z_k Y_k, \qquad
     *     \theta_e = \frac{F_{-1/2}(\eta)}{F_{1/2}(\eta)},
     * \f]
     * with \f$F_j\f$ the normalised (non-relativistic) Fermi-Dirac integrals. It therefore reduces to
     * WeakScreeningModel (\f$\theta_e = 1\f$) for non-degenerate electrons and to CD09 eq. (A1) (\f$\theta_e = 0\f$)
     * for degenerate ones. \f$\theta_e\f$ is evaluated from the electron degeneracy parameter
     * \f$\psi = F_{1/2}(\eta) = n_e \lambda_e^3 / 2\f$ (\f$\lambda_e\f$ the electron thermal wavelength) as
     * \f[
     *     \theta_e = \left(1 + a\psi + b\psi^{3/2} + c\psi^2\right)^{-1/3}, \qquad
     *     a = \frac{3}{2\sqrt{2}}, \quad b = -0.1915, \quad c = \frac{\pi}{6},
     * \f]
     * whose \f$a\f$ and \f$c\f$ reproduce the exact non-degenerate (\f$\theta_e \simeq 1 - \psi / 2\sqrt{2}\f$) and
     * degenerate (\f$\theta_e \simeq 3 / 2\eta\f$) limits and whose \f$b\f$ is fitted: it is within 0.7% of the exact
     * ratio for \f$-12 \le \eta \le 80\f$ and smooth, so it serves the CppAD tape and the analytic derivatives alike.
     *
     * @par Triple-alpha
     * Screened as the sum of its two sequential steps, \f$h(^4{\rm He}, ^4{\rm He}) + h(^4{\rm He}, ^8{\rm Be})\f$.
     * Because this model is non-linear in the charges this differs from the \f$3 \times 4\f$ prefactor shortcut of the
     * weak model, but the two agree in the weak limit (\f$2 \cdot 2 + 2 \cdot 4 = 12\f$). Reactions with any other
     * number of reactants are unscreened, as in WeakScreeningModel.
     *
     * @par Analytic derivatives
     * Every quantity depends on the state only through \f$L = \ln\Gamma_e\f$
     * (\f$\partial L/\partial\ln T = -1\f$, \f$\partial L/\partial\ln\rho = 1/3\f$, \f$\partial L/\partial Y_k = Z_k / 3Y_e\f$),
     * \f$\ln\tau_{ij}\f$ (\f$\partial/\partial\ln T = -1/3\f$) and \f$\ln C_{ij}\f$ (through \f$\zeta / Y_e\f$ and
     * \f$\theta_e(\psi)\f$, with \f$\psi \propto \rho Y_e T^{-3/2}\f$), so
     * \f$\partial\ln F_{ij}/\partial Y_k = \alpha_{ij} Z_k + \beta_{ij} Z_k^2\f$, which is cheap to fill.
     *
     * @par Speed
     * The composition- and state-independent constants of each screened pair (powers of the charges, the prefactor
     * of \f$\tau_{ij}\f$, the eq. 25 coefficients) are computed once and cached, and reactions sharing a pair (the same
     * charges and masses, e.g. \f$(p,\gamma)\f$, \f$(p,\alpha)\f$ and \f$(p,n)\f$ on one target) share its evaluation
     * within a call. The cache is keyed by reaction address but checks the reaction's charges and masses on every
     * lookup, so a reaction allocated at a reused address cannot pick up a stale entry; the model is safe to share
     * across threads.
     */
    class CD09ScreeningModel final : public ScreeningModel {
    public:
        [[nodiscard]] std::vector<double> calculateScreeningFactors(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho
        ) const override {
            return calculateFactors_impl<double>(reactions, species, Y, T9, rho);
        }

        /// Only the reactions flagged in @p needed are evaluated; the others get factor 1.
        [[nodiscard]] std::vector<double> calculateScreeningFactors(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho,
            const std::vector<char>& needed
        ) const override {
            return calculateFactors_impl<double>(reactions, species, Y, T9, rho, &needed);
        }

        [[nodiscard]] std::vector<ADDouble> calculateScreeningFactors(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<ADDouble> &Y,
            const ADDouble T9,
            const ADDouble rho
        ) const override {
            return calculateFactors_impl<ADDouble>(reactions, species, Y, T9, rho);
        }

        [[nodiscard]] double calculateScreeningFactor(
            const reaction::Reaction& reaction,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho
        ) const override {
            const PlasmaState<double> state = plasma_impl<double>(species, Y, T9, rho);
            const ResolvedReaction resolved = resolve(reaction);
            double h = 0.0;
            for (std::uint8_t n = 0; n < resolved.n; ++n) {
                h += evaluate_pair<double>(*resolved.pair[n], state).h;
            }
            return capped_factor<double>(h);
        }

        [[nodiscard]] std::optional<std::vector<double>> calculateScreeningLogDerivatives(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho
        ) const override {
            return log_derivatives_impl(reactions, species, Y, T9, rho, nullptr);
        }

        /// Only the reactions flagged in @p needed are evaluated; the others get zero derivatives.
        [[nodiscard]] std::optional<std::vector<double>> calculateScreeningLogDerivatives(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho,
            const std::vector<char>& needed
        ) const override {
            return log_derivatives_impl(reactions, species, Y, T9, rho, &needed);
        }

        [[nodiscard]] std::optional<std::vector<StateLogDerivatives>> calculateScreeningStateLogDerivatives(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho
        ) const override {
            return state_log_derivatives_impl(reactions, species, Y, T9, rho, nullptr);
        }

        /// Only the reactions flagged in @p needed are evaluated; the others get zero derivatives.
        [[nodiscard]] std::optional<std::vector<StateLogDerivatives>> calculateScreeningStateLogDerivatives(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho,
            const std::vector<char>& needed
        ) const override {
            return state_log_derivatives_impl(reactions, species, Y, T9, rho, &needed);
        }

    private:
        [[nodiscard]] std::vector<double> log_derivatives_impl(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho,
            const std::vector<char>* needed
        ) const {
            std::vector<double> dlnf_dY(reactions.size() * species.size(), 0.0);
            const PlasmaState<double> state = plasma_impl<double>(species, Y, T9, rho);
            if (state.valid == 0.0 || !(state.zeta > 0.0)) {
                return dlnf_dY; // every factor is pinned at 1 (matches WeakScreeningModel)
            }
            // d ln C / dY_k = (Z_k^2 + Z_k (theta + dtheta/dlnpsi)) / (2 zeta) - Z_k / (2 Y_e)
            const double dlnC_dZ = 0.5 * (state.theta + state.dtheta_dlnpsi) / state.zeta - 0.5 / state.Ye;
            const double dlnC_dZ2 = 0.5 / state.zeta;

            const std::vector<ResolvedReaction> resolved = resolve(reactions);
            PartialsCache& cache = state_memo(state, T9, rho);
            for (size_t r = 0; r < resolved.size(); ++r) {
                if (needed && !(*needed)[r]) {
                    continue;
                }
                const ReactionPartials d = reaction_partials(resolved[r], state, cache);
                if (!d.active) {
                    continue;
                }
                // dh/dY_k = A_L dL/dY_k + A_C dC/dY_k = alpha Z_k + beta Z_k^2
                const double alpha = d.A_L / (3.0 * state.Ye) + d.A_C_times_C * dlnC_dZ;
                const double beta = d.A_C_times_C * dlnC_dZ2;
                double* row = dlnf_dY.data() + r * species.size();
                for (size_t k = 0; k < species.size(); ++k) {
                    const double Z = species[k].m_z;
                    row[k] = alpha * Z + beta * Z * Z;
                }
            }
            return dlnf_dY;
        }

        [[nodiscard]] std::vector<StateLogDerivatives> state_log_derivatives_impl(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            const double T9,
            const double rho,
            const std::vector<char>* needed
        ) const {
            std::vector<StateLogDerivatives> derivatives(reactions.size());
            const PlasmaState<double> state = plasma_impl<double>(species, Y, T9, rho);
            if (state.valid == 0.0 || !(state.zeta > 0.0)) {
                return derivatives;
            }
            // ln C depends on T and rho only through theta_e(psi), psi proportional to rho T^{-3/2}
            const double dlnC_dlnpsi = 0.5 * state.Ye * state.dtheta_dlnpsi / state.zeta;

            const std::vector<ResolvedReaction> resolved = resolve(reactions);
            PartialsCache& cache = state_memo(state, T9, rho);
            for (size_t r = 0; r < resolved.size(); ++r) {
                if (needed && !(*needed)[r]) {
                    continue;
                }
                if (const ReactionPartials d = reaction_partials(resolved[r], state, cache); d.active) {
                    // dL/dlnT = -1, dlntau/dlnT = -1/3, dlnpsi/dlnT = -3/2; dL/dlnrho = 1/3, dlnpsi/dlnrho = 1
                    derivatives[r] = {
                        .dlnT = -d.A_L - d.A_tau / 3.0 - 1.5 * d.A_C_times_C * dlnC_dlnpsi,
                        .dlnRho = d.A_L / 3.0 + d.A_C_times_C * dlnC_dlnpsi
                    };
                }
            }
            return derivatives;
        }

    private:
        // ---------------------------------------------------------------------------------------------
        // Physical constants (CGS, from libconstants)
        // ---------------------------------------------------------------------------------------------
        struct PhysicalConstants {
            double e2;          ///< \f$e^2\f$ [erg cm] (Gaussian units: \f$e^2 = \alpha\hbar c\f$)
            double kB;          ///< Boltzmann constant [erg/K]
            double hbar;        ///< Reduced Planck constant [erg s]
            double N_A;         ///< Avogadro's number [1/mol]
            double u;           ///< Atomic mass unit [g]
            double psiCoeff;    ///< \f$\psi = {\rm psiCoeff}\,\rho Y_e T^{-3/2}\f$, psiCoeff \f$= (N_A/2)\, h^3 (2\pi m_e k_B)^{-3/2}\f$
        };

        /// Read once from libconstants. Its elementary charge is in coulombs (SI), so the Gaussian e^2 is built from
        /// the fine-structure constant instead; the units of every input are checked.
        [[nodiscard]] static const PhysicalConstants& constants() {
            static const PhysicalConstants k = [] {
                const auto& table = fourdst::constant::Constants::getInstance();
                const auto get = [&table](const std::string& name, const std::string& unit) {
                    const fourdst::constant::Constant constant = table.get(name);
                    if (constant.unit != unit) {
                        throw std::runtime_error("CD09ScreeningModel: libconstants gives " + name + " in '" + constant.unit +
                                                 "', expected '" + unit + "'.");
                    }
                    return constant.value;
                };
                const double hbar = get("hbar", "erg s");
                const double kB = get("kB", "erg / K");
                const double N_A = get("N_a", "1 / mol");
                const double h = get("h", "erg s");
                const double me = get("mE", "g");
                return PhysicalConstants{
                    .e2 = get("alpha", "") * hbar * get("c", "cm / s"),
                    .kB = kB,
                    .hbar = hbar,
                    .N_A = N_A,
                    .u = get("u", "g"),
                    .psiCoeff = 0.5 * N_A * h * h * h / std::pow(2.0 * std::numbers::pi * me * kB, 1.5),
                };
            }();
            return k;
        }

        /// Weak-screening coefficient; must be identical to WeakScreeningModel::prefactor_impl.
        static constexpr double WEAK_COEFF = 0.188;
        /// Low temperature cutoff; must be identical to WeakScreeningModel::prefactor_impl.
        static constexpr double LOW_T9_THRESHOLD = 1e-9;

        /// Potekhin & Chabrier (2000) OCP free-energy fit coefficients, as used by CD09 eq. (24).
        static constexpr double A1 = -0.907;
        static constexpr double A2 = 0.62954;
        static constexpr double B1 = 0.00456;
        static constexpr double B2 = 211.6;
        static constexpr double B3 = -1.0e-4;
        static constexpr double B4 = 0.00462;
        static constexpr double SQRT3 = std::numbers::sqrt3;
        inline static const double SQRT_A2 = std::sqrt(A2);
        inline static const double A3 = -SQRT3 / 2.0 - A1 / std::sqrt(A2);

        /// CD09 eq. (25) coefficients (Gamma-independent parts).
        static constexpr double C1_COEFF = 0.013;
        static constexpr double C2_COEFF = 0.406;
        static constexpr double C3_COEFF = 0.062;
        static constexpr double C3_GAMMA = 1.8; ///< c3 contains + C3_GAMMA / Gamma_ij

        /// Electron degeneracy factor \f$\theta_e = (1 + a\psi + b\psi^{3/2} + c\psi^2)^{-1/3}\f$ (see the class description).
        static constexpr double THETA_A = 1.0606601717798212; ///< \f$a = 3/(2\sqrt{2})\f$: exact non-degenerate slope
        static constexpr double THETA_B = -0.1915;            ///< \f$b\f$, fitted to \f$F_{-1/2}/F_{1/2}\f$ for \f$-12 \le \eta \le 80\f$
        static constexpr double THETA_C = std::numbers::pi / 6.0; ///< \f$c = \Gamma(5/2)^2 / 1.5^3 = \pi/6\f$: exact degenerate limit

        /// Guard on exp(h); far beyond any physical stellar-interior value, keeps the factor finite.
        static constexpr double H_CAP = 200.0;

        // ---------------------------------------------------------------------------------------------
        // State and pair data
        // ---------------------------------------------------------------------------------------------
        template <IsArithmeticOrAD T>
        struct PlasmaState {
            T Ye;             ///< \f$Y_e = \sum_k Z_k Y_k\f$
            T zeta;           ///< \f$\zeta = \sum_k Z_k^2 Y_k + \theta_e \sum_k Z_k Y_k\f$ (Salpeter's \f$\zeta\f$ with electron degeneracy)
            T theta;          ///< electron degeneracy factor \f$\theta_e\f$
            T dtheta_dlnpsi;  ///< \f$d\theta_e / d\ln\psi\f$
            T T_K;            ///< temperature [K]
            T Gamma_e;        ///< electron coupling parameter \f$\Gamma_e\f$
            T cbrtT;          ///< \f$T^{1/3}\f$ (shared by every pair's \f$\tau_{ij}\f$)
            T Gamma_e32;      ///< \f$\Gamma_e^{3/2}\f$ (shared by every pair's \f$h^{\rm lin}_{ij}\f$)
            T weakPrefactor;  ///< \f$0.188 \sqrt{\rho\zeta / T_6^3}\f$
            T valid;          ///< 1 if screening applies, 0 otherwise
        };

        /// Composition- and state-independent constants for one screened pair.
        struct PairSpec {
            std::uint32_t id = 0;       ///< Index of this pair in the cache (for per-call sharing)
            double Zi = 0.0;
            double Zj = 0.0;
            double ZiZj = 0.0;
            double gammaIjFactor = 0.0; ///< \f$\Gamma_{ij}/\Gamma_e = 2 Z_i Z_j / (Z_i^{1/3} + Z_j^{1/3})\f$
            double zi53 = 0.0;          ///< \f$Z_i^{5/3}\f$
            double zj53 = 0.0;          ///< \f$Z_j^{5/3}\f$
            double zc53 = 0.0;          ///< \f$(Z_i + Z_j)^{5/3}\f$
            double eta = 0.0;           ///< \f$(Z_i + Z_j)^{5/2} - Z_i^{5/2} - Z_j^{5/2}\f$
            double tauCoeff = 0.0;      ///< \f$\tau_{ij} T^{1/3}\f$
            double c1 = 0.0;
            double c2 = 0.0;
            double c3a = 0.0;           ///< Gamma-independent part of c3
        };

        /// What identifies a pair: its charges and masses (the masses enter tau through the reduced mass).
        struct PairKey {
            int Zi = 0;
            int Zj = 0;
            double Ai = 0.0;
            double Aj = 0.0;
            bool operator==(const PairKey&) const = default;
        };
        struct PairKeyHash {
            size_t operator()(const PairKey& k) const {
                size_t h = std::hash<int>{}(k.Zi) * 1000003u ^ std::hash<int>{}(k.Zj);
                h = h * 1000003u ^ std::hash<double>{}(k.Ai);
                return h * 1000003u ^ std::hash<double>{}(k.Aj);
            }
        };

        /// The screened pairs of one reaction, as keys (cheap to form; no powers or roots).
        struct PairKeys {
            std::array<PairKey, 2> keys{};
            std::uint8_t n = 0;
        };

        /// A reaction resolved against the cache: pointers to its pairs' constants (stable: kept in a deque).
        struct ResolvedReaction {
            std::array<const PairSpec*, 2> pair{};
            std::uint8_t n = 0;
        };

        /// Cached resolution of one reaction, with the reactants it was made from (checked on every lookup).
        struct CachedReaction {
            ResolvedReaction resolved;
            std::array<int, 3> Z{};
            std::array<double, 3> A{};
            std::uint8_t num_reactants = 0;
        };

        template <IsArithmeticOrAD T>
        struct PairEval {
            T h;
            T S;        ///< f0 bracket
            T P;        ///< (C + Gamma_ij^2) / (1 + Gamma_ij^2)
            T C;        ///< h_weak / h_lin
            T Gamma_ij;
            T tau;
            T zeta;
            T t3;       ///< t_ij^3
            T gi, gj, gc;
        };

        struct ReactionPartials {
            bool active = false;
            double h = 0.0;            ///< exponent (summed over the reaction's pairs)
            double A_L = 0.0;          ///< dh/dln Gamma_e at fixed tau, C
            double A_tau = 0.0;        ///< dh/dln tau at fixed Gamma_e, C
            double A_C_times_C = 0.0;  ///< C * dh/dC (so that dh = ... + A_C_times_C * dln C)
        };

        /**
         * Each pair's exponent and partials at one plasma state, kept per thread: the engine (and BlendedScreeningModel)
         * asks for the factors, composition derivatives and state derivatives at the same state in separate calls, and
         * each pair is shared by several reactions. Keyed by the model instance and everything the exponents depend on
         * (T9, rho, Y_e and zeta), so a different state, or another model instance, starts afresh.
         */
        struct PartialsCache {
            std::uint64_t instance = 0;
            double T9 = std::numeric_limits<double>::quiet_NaN();
            double rho = std::numeric_limits<double>::quiet_NaN();
            double Ye = std::numeric_limits<double>::quiet_NaN();
            double zeta = std::numeric_limits<double>::quiet_NaN();
            std::vector<double> h;                  ///< exponent of each pair
            std::vector<char> haveH;
            std::vector<ReactionPartials> partials; ///< exponent and partials of each pair
            std::vector<char> have;

            void ensure(const size_t size) {
                if (h.size() < size) {
                    h.resize(size);
                    haveH.resize(size, 0);
                    partials.resize(size);
                    have.resize(size, 0);
                }
            }
        };

        /// This thread's memo for the given state of this model (reset if the state or the instance differs).
        [[nodiscard]] PartialsCache& state_memo(const PlasmaState<double>& s, const double T9, const double rho) const {
            thread_local PartialsCache memo;
            if (memo.instance != m_instance || memo.T9 != T9 || memo.rho != rho || memo.Ye != s.Ye || memo.zeta != s.zeta) {
                memo.instance = m_instance;
                memo.T9 = T9;
                memo.rho = rho;
                memo.Ye = s.Ye;
                memo.zeta = s.zeta;
                std::ranges::fill(memo.haveH, 0);
                std::ranges::fill(memo.have, 0);
            }
            return memo;
        }

        inline static std::atomic<std::uint64_t> s_next_instance{1};
        const std::uint64_t m_instance = s_next_instance.fetch_add(1); ///< Distinguishes instances in the per-thread memo

        mutable std::shared_mutex m_cache_mutex;
        mutable std::unordered_map<const reaction::Reaction*, CachedReaction> m_reaction_cache;
        mutable std::unordered_map<PairKey, const PairSpec*, PairKeyHash> m_pair_index;
        mutable std::deque<PairSpec> m_pairs; ///< deque: addresses stay valid as pairs are added

        // ---------------------------------------------------------------------------------------------
        // Helpers
        // ---------------------------------------------------------------------------------------------
        template <IsArithmeticOrAD T>
        [[nodiscard]] static T cbrt_t(const T& x) {
            if constexpr (std::is_same_v<T, double>) {
                return std::cbrt(x);
            } else {
                return CppAD::pow(x, 1.0 / 3.0);
            }
        }

        /**
         * Potekhin & Chabrier (2000) OCP free energy per ion (CD09 eq. 24):
         * \f[
         *     f_0(\Gamma) = A_1\left[\sqrt{\Gamma(A_2 + \Gamma)} - A_2 \ln\left(\sqrt{\Gamma/A_2} + \sqrt{1 + \Gamma/A_2}\right)\right]
         *                 + 2A_3\left[\sqrt{\Gamma} - \arctan\sqrt{\Gamma}\right]
         *                 + B_1\left[\Gamma - B_2 \ln(1 + \Gamma/B_2)\right] + \frac{B_3}{2} \ln(1 + \Gamma^2/B_4)
         * \f]
         */
        template <IsArithmeticOrAD T>
        [[nodiscard]] static T f0(const T& g) {
            const T sg = CppAD::sqrt(g);
            const T asinh_term = CppAD::log(sg / SQRT_A2 + CppAD::sqrt(T(1.0) + g / A2));
            const T term1 = CppAD::sqrt(g * (T(A2) + g)) - A2 * asinh_term;
            const T term2 = sg - CppAD::atan(sg);
            const T term3 = g - B2 * CppAD::log(T(1.0) + g / B2);
            const T term4 = CppAD::log(T(1.0) + g * g / B4);
            return A1 * term1 + 2.0 * A3 * term2 + B1 * term3 + 0.5 * B3 * term4;
        }

        /// \f$ f_0'(\Gamma) = A_1\sqrt{\Gamma/(A_2 + \Gamma)} + A_3\sqrt{\Gamma}/(1 + \Gamma) + B_1\Gamma/(B_2 + \Gamma) + B_3\Gamma/(B_4 + \Gamma^2) \f$
        [[nodiscard]] static double f0_prime(const double g) {
            const double sg = std::sqrt(g);
            return A1 * std::sqrt(g / (A2 + g)) + A3 * sg / (1.0 + g) + B1 * g / (B2 + g) + B3 * g / (B4 + g * g);
        }

        [[nodiscard]] static PairSpec make_pair(const PairKey& key, const std::uint32_t id) {
            const double Zi = key.Zi;
            const double Zj = key.Zj;
            PairSpec p;
            p.id = id;
            p.Zi = Zi;
            p.Zj = Zj;
            p.ZiZj = Zi * Zj;
            p.gammaIjFactor = 2.0 * Zi * Zj / (std::cbrt(Zi) + std::cbrt(Zj));
            p.zi53 = std::pow(Zi, 5.0 / 3.0);
            p.zj53 = std::pow(Zj, 5.0 / 3.0);
            p.zc53 = std::pow(Zi + Zj, 5.0 / 3.0);
            p.eta = std::pow(Zi + Zj, 2.5) - std::pow(Zi, 2.5) - std::pow(Zj, 2.5);

            const PhysicalConstants& k = constants();
            const double mu = key.Ai * key.Aj / (key.Ai + key.Aj) * k.u;
            p.tauCoeff = std::cbrt(
                27.0 * std::numbers::pi * std::numbers::pi * mu * k.e2 * k.e2 * Zi * Zi * Zj * Zj /
                (2.0 * k.hbar * k.hbar * k.kB)
            );

            const double y = 4.0 * Zi * Zj / ((Zi + Zj) * (Zi + Zj));
            p.c1 = C1_COEFF * y * y; // CD09 eq. (25): c1 = 0.013 y^2
            p.c2 = C2_COEFF * std::pow(y, 0.14);
            p.c3a = C3_COEFF * std::pow(y, 0.19);
            return p;
        }

        /// The screened pairs for a reaction (same reaction classification as WeakScreeningModel::H_impl).
        [[nodiscard]] static PairKeys pair_keys(const reaction::Reaction& reaction) {
            PairKeys list;
            const auto& reactants = reaction.reactants();
            if (reactants.size() == 2) {
                const auto& a = reactants[0];
                const auto& b = reactants[1];
                if (a.m_z > 0 && b.m_z > 0) {
                    list.keys[0] = {static_cast<int>(a.m_z), static_cast<int>(b.m_z), a.mass(), b.mass()};
                    list.n = 1;
                }
                return list;
            }
            const bool isTripleAlpha = (
                reactants.size() == 3 &&
                reactants[0].m_z == 2 &&
                reactants[1].m_z == 2 &&
                reactants[2].m_z == 2 &&
                reactants[0] == reactants[1] &&
                reactants[1] == reactants[2]
            );
            if (isTripleAlpha) {
                const double mHe = reactants[0].mass();
                list.keys[0] = {2, 2, mHe, mHe};                              // He4 + He4 -> Be8
                list.keys[1] = {2, 4, mHe, fourdst::atomic::Be_8.mass()};     // He4 + Be8 -> C12
                list.n = 2;
            }
            return list;
        }

        /// Whether a cached entry was made from this reaction's reactants (guards against reused addresses).
        [[nodiscard]] static bool matches(const CachedReaction& cached, const reaction::Reaction& reaction) {
            const auto& reactants = reaction.reactants();
            if (reactants.size() != cached.num_reactants) {
                return false;
            }
            for (size_t i = 0; i < reactants.size() && i < 3; ++i) {
                if (static_cast<int>(reactants[i].m_z) != cached.Z[i] || reactants[i].mass() != cached.A[i]) {
                    return false;
                }
            }
            return true;
        }

        /// Build (under the exclusive lock) the cache entry for a reaction.
        [[nodiscard]] CachedReaction insert_locked(const reaction::Reaction& reaction) const {
            CachedReaction entry;
            const auto& reactants = reaction.reactants();
            entry.num_reactants = static_cast<std::uint8_t>(reactants.size());
            for (size_t i = 0; i < reactants.size() && i < 3; ++i) {
                entry.Z[i] = static_cast<int>(reactants[i].m_z);
                entry.A[i] = reactants[i].mass();
            }
            const PairKeys keys = pair_keys(reaction);
            entry.resolved.n = keys.n;
            for (std::uint8_t n = 0; n < keys.n; ++n) {
                auto it = m_pair_index.find(keys.keys[n]);
                if (it == m_pair_index.end()) {
                    m_pairs.push_back(make_pair(keys.keys[n], static_cast<std::uint32_t>(m_pairs.size())));
                    it = m_pair_index.emplace(keys.keys[n], &m_pairs.back()).first;
                }
                entry.resolved.pair[n] = it->second;
            }
            m_reaction_cache[&reaction] = entry;
            return entry;
        }

        /// Resolve every reaction of a set against the cache (one shared lock; an exclusive one only for new reactions).
        [[nodiscard]] std::vector<ResolvedReaction> resolve(const reaction::ReactionSet& reactions) const {
            std::vector<ResolvedReaction> resolved(reactions.size());
            std::vector<size_t> missing;
            {
                std::shared_lock lock(m_cache_mutex);
                size_t r = 0;
                for (const auto& reaction : reactions) {
                    const auto it = m_reaction_cache.find(reaction.get());
                    if (it != m_reaction_cache.end() && matches(it->second, *reaction)) {
                        resolved[r] = it->second.resolved;
                    } else {
                        missing.push_back(r);
                    }
                    ++r;
                }
            }
            if (!missing.empty()) {
                std::unique_lock lock(m_cache_mutex);
                for (const size_t r : missing) {
                    resolved[r] = insert_locked(reactions[r]).resolved;
                }
            }
            return resolved;
        }

        [[nodiscard]] ResolvedReaction resolve(const reaction::Reaction& reaction) const {
            {
                std::shared_lock lock(m_cache_mutex);
                const auto it = m_reaction_cache.find(&reaction);
                if (it != m_reaction_cache.end() && matches(it->second, reaction)) {
                    return it->second.resolved;
                }
            }
            std::unique_lock lock(m_cache_mutex);
            return insert_locked(reaction).resolved;
        }

        [[nodiscard]] size_t num_pairs() const {
            std::shared_lock lock(m_cache_mutex);
            return m_pairs.size();
        }

        template <IsArithmeticOrAD T>
        [[nodiscard]] static PlasmaState<T> plasma_impl(
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<T>& Y,
            const T T9,
            const T rho
        ) {
            const T zero(0.0);
            const T one(1.0);
            const T tiny(1e-300);
            const T t9_floor(LOW_T9_THRESHOLD);

            PlasmaState<T> s;
            s.Ye = zero;
            T zetaIon = zero;
            for (size_t i = 0; i < species.size(); ++i) {
                const T Z = species[i].m_z;
                s.Ye += Z * Y[i];
                zetaIon += Z * Z * Y[i];
            }

            const PhysicalConstants& k = constants();
            const T T9_safe = CppAD::CondExpLe(T9, t9_floor, t9_floor, T9);
            const T rho_safe = CppAD::CondExpLe(rho, tiny, tiny, rho);
            const T Ye_safe = CppAD::CondExpLe(s.Ye, tiny, tiny, s.Ye);
            s.T_K = T9_safe * 1.0e9;

            // Electron degeneracy: psi = F_{1/2}(eta) = n_e lambda_e^3 / 2, theta_e = F_{-1/2} / F_{1/2}
            const T psi = k.psiCoeff * rho_safe * Ye_safe / (s.T_K * CppAD::sqrt(s.T_K));
            const T sqrtPsi = CppAD::sqrt(psi);
            const T poly = 1.0 + THETA_A * psi + THETA_B * psi * sqrtPsi + THETA_C * psi * psi;
            s.theta = 1.0 / cbrt_t<T>(poly);
            const T theta4 = s.theta * s.theta * s.theta * s.theta;
            s.dtheta_dlnpsi = -(theta4 / 3.0) * (THETA_A * psi + 1.5 * THETA_B * psi * sqrtPsi + 2.0 * THETA_C * psi * psi);
            s.zeta = zetaIon + s.theta * s.Ye;

            // Same low-temperature behaviour as WeakScreeningModel; additionally unscreened for rho, Y_e, zeta <= 0
            const T lowTFlag = CppAD::CondExpLt(T9, t9_floor, zero, one);
            const T rhoFlag = CppAD::CondExpLe(rho, zero, zero, one);
            const T yeFlag = CppAD::CondExpLe(s.Ye, zero, zero, one);
            const T zetaFlag = CppAD::CondExpLe(s.zeta, zero, zero, one);
            s.valid = lowTFlag * rhoFlag * yeFlag * zetaFlag;
            const T zeta_safe = CppAD::CondExpLe(s.zeta, tiny, tiny, s.zeta);

            const T ne = rho_safe * k.N_A * Ye_safe;
            const T a_e = cbrt_t<T>(T(3.0) / (4.0 * std::numbers::pi * ne));
            s.Gamma_e = k.e2 / (a_e * k.kB * s.T_K);
            s.Gamma_e32 = s.Gamma_e * CppAD::sqrt(s.Gamma_e);
            s.cbrtT = cbrt_t<T>(s.T_K);

            const T T6 = T9_safe * 1000.0;
            s.weakPrefactor = WEAK_COEFF * CppAD::sqrt(rho_safe / (T6 * T6 * T6)) * CppAD::sqrt(zeta_safe);
            return s;
        }

        template <IsArithmeticOrAD T>
        [[nodiscard]] static PairEval<T> evaluate_pair(const PairSpec& p, const PlasmaState<T>& s) {
            PairEval<T> e;
            e.Gamma_ij = s.Gamma_e * p.gammaIjFactor;
            e.tau = p.tauCoeff / s.cbrtT;
            e.zeta = 3.0 * e.Gamma_ij / e.tau;

            // c3 zeta^3 = c3a zeta^3 + 1.8 zeta^3 / Gamma_ij = c3a zeta^3 + 5.4 zeta^2 / tau (no 1/Gamma singularity)
            const T z2 = e.zeta * e.zeta;
            e.t3 = 1.0 + p.c1 * e.zeta + p.c2 * z2 + p.c3a * z2 * e.zeta + 3.0 * C3_GAMMA * z2 / e.tau;
            const T t = cbrt_t<T>(e.t3);

            e.gi = s.Gamma_e * p.zi53 / t;
            e.gj = s.Gamma_e * p.zj53 / t;
            e.gc = s.Gamma_e * p.zc53 / t;
            e.S = f0<T>(e.gi) + f0<T>(e.gj) - f0<T>(e.gc);

            // Weak-limit normalisation: C = h_weak / h_lin, h_lin = Gamma_e^{3/2} eta / sqrt(3)
            const T hLin = s.Gamma_e32 * p.eta / SQRT3;
            const T hWeak = p.ZiZj * s.weakPrefactor;
            e.C = hWeak / hLin;

            const T G2 = e.Gamma_ij * e.Gamma_ij;
            e.P = (e.C + G2) / (1.0 + G2);
            e.h = s.valid * e.P * e.S;
            return e;
        }

        template <IsArithmeticOrAD T>
        [[nodiscard]] static T capped_factor(const T& h) {
            const T cap(H_CAP);
            return CppAD::exp(CppAD::CondExpGe(h, cap, cap, h));
        }

        /// Partials of one pair's exponent with respect to ln Gamma_e, ln tau and ln C.
        [[nodiscard]] static ReactionPartials pair_partials(const PairSpec& p, const PlasmaState<double>& s) {
            ReactionPartials d;
            const PairEval<double> e = evaluate_pair<double>(p, s);
            d.h = e.h;

            // dS = G (dL - dln t),  G = sum_k +/- g_k f0'(g_k)
            const double G = e.gi * f0_prime(e.gi) + e.gj * f0_prime(e.gj) - e.gc * f0_prime(e.gc);

            // t^3 = 1 + c1 z + c2 z^2 + c3a z^3 + 5.4 z^2 / tau,  z = 3 Gamma_ij / tau
            //   dt^3/dL      = c1 z + 2 c2 z^2 + 3 c3a z^3 + 10.8 z^2 / tau
            //   dt^3/dlntau  = -(that) - 5.4 z^2 / tau
            const double z = e.zeta;
            const double explicitTau = 3.0 * C3_GAMMA * z * z / e.tau;
            const double q = p.c1 * z + 2.0 * p.c2 * z * z + 3.0 * p.c3a * z * z * z + 2.0 * explicitTau;
            const double dlnt_dL = q / (3.0 * e.t3);
            const double dlnt_dlnTau = -(q + explicitTau) / (3.0 * e.t3);

            const double dS_dL = G * (1.0 - dlnt_dL);
            const double dS_dlnTau = -G * dlnt_dlnTau;

            // P = (C + x) / (1 + x), x = Gamma_ij^2, dx/dL = 2x
            const double x = e.Gamma_ij * e.Gamma_ij;
            const double onePx = 1.0 + x;
            const double dP_dL = 2.0 * x * (1.0 - e.C) / (onePx * onePx);
            const double dP_dC = 1.0 / onePx;

            d.A_L = s.valid * (dP_dL * e.S + e.P * dS_dL);
            d.A_tau = s.valid * (e.P * dS_dlnTau);
            d.A_C_times_C = s.valid * (e.C * dP_dC * e.S);
            return d;
        }

        /// Partial derivatives of a reaction's (summed) exponent, sharing each pair's partials within a call.
        [[nodiscard]] ReactionPartials reaction_partials(
            const ResolvedReaction& reaction,
            const PlasmaState<double>& s,
            PartialsCache& cache
        ) const {
            ReactionPartials d;
            for (std::uint8_t n = 0; n < reaction.n; ++n) {
                const PairSpec& p = *reaction.pair[n];
                if (p.id >= cache.have.size()) {
                    cache.ensure(std::max(num_pairs(), static_cast<size_t>(p.id) + 1));
                }
                if (!cache.have[p.id]) {
                    cache.partials[p.id] = pair_partials(p, s);
                    cache.have[p.id] = 1;
                    cache.h[p.id] = cache.partials[p.id].h;
                    cache.haveH[p.id] = 1;
                }
                const ReactionPartials& pd = cache.partials[p.id];
                d.h += pd.h;
                d.A_L += pd.A_L;
                d.A_tau += pd.A_tau;
                d.A_C_times_C += pd.A_C_times_C;
            }
            // Above the cap the factor is constant, so all derivatives vanish (as in WeakScreeningModel)
            d.active = reaction.n > 0 && s.valid != 0.0 && d.h < H_CAP;
            return d;
        }

        template <IsArithmeticOrAD T>
        [[nodiscard]] std::vector<T> calculateFactors_impl(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<T>& Y,
            const T T9,
            const T rho,
            const std::vector<char>* needed = nullptr
        ) const {
            const PlasmaState<T> state = plasma_impl<T>(species, Y, T9, rho);
            const std::vector<ResolvedReaction> resolved = resolve(reactions);

            if constexpr (std::is_same_v<T, double>) {
                // Each pair is evaluated once per state (shared with the derivative calls at the same state)
                PartialsCache& cache = state_memo(state, T9, rho);
                std::vector<double> factors;
                factors.reserve(resolved.size());
                for (size_t r = 0; r < resolved.size(); ++r) {
                    const ResolvedReaction& reaction = resolved[r];
                    if (needed && !(*needed)[r]) {
                        factors.push_back(1.0);
                        continue;
                    }
                    double h = 0.0;
                    for (std::uint8_t n = 0; n < reaction.n; ++n) {
                        const PairSpec& p = *reaction.pair[n];
                        if (p.id >= cache.haveH.size()) {
                            cache.ensure(std::max(num_pairs(), static_cast<size_t>(p.id) + 1));
                        }
                        if (!cache.haveH[p.id]) {
                            cache.h[p.id] = evaluate_pair<double>(p, state).h;
                            cache.haveH[p.id] = 1;
                        }
                        h += cache.h[p.id];
                    }
                    factors.push_back(capped_factor<double>(h));
                }
                return factors;
            }

            // CppAD path (recording a tape): each pair is evaluated once per call, however many reactions share it
            std::vector<T> pairExponent(num_pairs());
            std::vector<char> have(pairExponent.size(), 0);
            std::vector<T> factors;
            factors.reserve(resolved.size());
            for (const ResolvedReaction& reaction : resolved) {
                T h(0.0);
                for (std::uint8_t n = 0; n < reaction.n; ++n) {
                    const PairSpec& p = *reaction.pair[n];
                    if (p.id >= have.size()) {
                        pairExponent.resize(p.id + 1);
                        have.resize(p.id + 1, 0);
                    }
                    if (!have[p.id]) {
                        pairExponent[p.id] = evaluate_pair<T>(p, state).h;
                        have[p.id] = 1;
                    }
                    h += pairExponent[p.id];
                }
                factors.push_back(capped_factor<T>(h));
            }
            return factors;
        }
    };
}
