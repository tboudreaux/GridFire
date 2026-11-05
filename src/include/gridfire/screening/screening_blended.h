#pragma once

#include "gridfire/screening/screening_abstract.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/reaction/reaction.h"

#include "fourdst/atomic/atomicSpecies.h"

#include <memory>
#include <optional>
#include <vector>

namespace gridfire::screening {

    /**
     * @class BlendedScreeningModel
     * @brief Screening that uses a weak-coupling prescription where it is valid and a strong-coupling one where it is
     *        not, with a smooth handover between them.
     *
     * The model is composed of two screening models, a weak-side one (normally WeakScreeningModel, Salpeter 1954) and
     * a strong-side one (normally CD09ScreeningModel, Chugunov & DeWitt 2009), and blends their exponents
     * \f$h = \ln f\f$ reaction by reaction:
     * \f[
     *     \ln f = (1 - w)\, H + w\, h_{\rm strong}, \qquad H = \ln f_{\rm weak},
     * \f]
     * \f[
     *     w = s^2 (3 - 2s), \qquad
     *     s = \min\left(1, \max\left(0, \frac{\ln(H / H_{\rm weak})}{\ln(H_{\rm strong} / H_{\rm weak})}\right)\right),
     * \f]
     * with \f$[H_{\rm weak}, H_{\rm strong}]\f$ the ScreeningBlendWindow. The handover variable is the weak
     * prescription's own exponent: for Salpeter screening \f$H = Z_1 Z_2 e^2 / \lambda_D k_B T\f$ is exactly the small
     * parameter of the weak-screening expansion, so the weak result is used while it is small and is phased out as it
     * grows. \f$w\f$ and \f$dw/ds\f$ are continuous, so \f$\ln f\f$ and its derivatives are continuous across the
     * window edges.
     *
     * @par Why blend
     * Weak (Salpeter) screening is accurate in weakly coupled plasmas such as the solar core, where it agrees with
     * mean-field calculations to about a percent (Gruzinov & Bahcall 1998), but it grows without bound as coupling
     * increases and is far too large in dense matter. CD09 is accurate at strong coupling but its eq. (A4)
     * interpolation under-screens heavier pairs in weakly coupled plasmas (about 14% below Salpeter for
     * \f$^{14}{\rm N}(p,\gamma)\f$ in the Sun). With the default window (\f$H\f$ from 0.5 to 1.5) every pp-chain and
     * CNO reaction in the Sun uses Salpeter, and dense or degenerate matter (low-mass stellar cores, red-giant and
     * white-dwarf cores) uses CD09.
     *
     * @par Derivatives
     * With \f$\partial_x\f$ any of \f$Y_k\f$, \f$\ln T\f$, \f$\ln\rho\f$,
     * \f[
     *     \partial_x \ln f = (1 - w)\, \partial_x H + w\, \partial_x h_{\rm strong}
     *                       + (h_{\rm strong} - H)\, \frac{6 s (1 - s)}{\ln(H_{\rm strong} / H_{\rm weak})}\,
     *                         \frac{\partial_x H}{H},
     * \f]
     * so analytic derivatives are available whenever both models provide them.
     */
    class BlendedScreeningModel final : public ScreeningModel {
    public:
        /**
         * @param weak Prescription used at weak coupling; its exponent \f$\ln f\f$ is the handover variable.
         * @param strong Prescription used at strong coupling.
         * @param window Where the handover happens (see ScreeningBlendWindow).
         * @throws std::invalid_argument for a missing model, or a window that is not
         *         0 < weakExponent < strongExponent <= 2.
         */
        BlendedScreeningModel(
            std::unique_ptr<ScreeningModel> weak,
            std::unique_ptr<ScreeningModel> strong,
            const ScreeningBlendWindow& window = {}
        );

        [[nodiscard]] std::vector<double> calculateScreeningFactors(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] std::vector<ADDouble> calculateScreeningFactors(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<ADDouble>& Y,
            ADDouble T9,
            ADDouble rho
        ) const override;

        [[nodiscard]] double calculateScreeningFactor(
            const reaction::Reaction& reaction,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] std::optional<std::vector<double>> calculateScreeningLogDerivatives(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] std::optional<std::vector<StateLogDerivatives>> calculateScreeningStateLogDerivatives(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const override;

        /// The handover window.
        [[nodiscard]] const ScreeningBlendWindow& window() const { return m_window; }

        /// Blend weight \f$w\f$ for a weak-side exponent \f$H\f$ (0: weak only, 1: strong only).
        [[nodiscard]] double weight(double H) const;

    private:
        std::unique_ptr<ScreeningModel> m_weak;
        std::unique_ptr<ScreeningModel> m_strong;
        ScreeningBlendWindow m_window;
        double m_lnWeak;   ///< \f$\ln H_{\rm weak}\f$
        double m_lnWidth;  ///< \f$\ln(H_{\rm strong} / H_{\rm weak})\f$

        /// \f$w\f$ and \f$dw/d\ln H\f$ for a weak-side exponent \f$H\f$
        struct Weight {
            double w = 0.0;
            double dw_dlnH = 0.0;
        };
        [[nodiscard]] Weight weight_and_slope(double H) const;

        /// Whether any reaction's weak-side exponent is above the window's lower edge (needs the strong model)
        [[nodiscard]] bool any_blended(const std::vector<double>& weakFactors) const;

        /// Which reactions need the strong prescription (weak-side exponent above the window's lower edge)
        [[nodiscard]] std::vector<char> needs_strong(const std::vector<double>& weakFactors) const;

        template <typename T>
        [[nodiscard]] std::vector<T> blend_factors(const std::vector<T>& weak, const std::vector<T>& strong) const;
    };
}
