#include "gridfire/screening/screening_blended.h"

#include "cppad/cppad.hpp"

#include <cmath>
#include <format>
#include <stdexcept>
#include <utility>

namespace gridfire::screening {
    BlendedScreeningModel::BlendedScreeningModel(
        std::unique_ptr<ScreeningModel> weak,
        std::unique_ptr<ScreeningModel> strong,
        const ScreeningBlendWindow& window
    ) :
    m_weak(std::move(weak)),
    m_strong(std::move(strong)),
    m_window(window) {
        if (!m_weak || !m_strong) {
            throw std::invalid_argument("BlendedScreeningModel needs both a weak and a strong screening model.");
        }
        if (!(window.weakExponent > 0.0) || !(window.strongExponent > window.weakExponent) || !(window.strongExponent <= 2.0)) {
            throw std::invalid_argument(std::format(
                "BlendedScreeningModel: the blend window must satisfy 0 < weakExponent < strongExponent <= 2 "
                "(2 is where weak screening is capped); got [{}, {}].", window.weakExponent, window.strongExponent));
        }
        m_lnWeak = std::log(window.weakExponent);
        m_lnWidth = std::log(window.strongExponent / window.weakExponent);
    }

    BlendedScreeningModel::Weight BlendedScreeningModel::weight_and_slope(const double H) const {
        if (!(H > m_window.weakExponent)) {
            return {0.0, 0.0};
        }
        if (H >= m_window.strongExponent) {
            return {1.0, 0.0};
        }
        const double s = (std::log(H) - m_lnWeak) / m_lnWidth;
        return {s * s * (3.0 - 2.0 * s), 6.0 * s * (1.0 - s) / m_lnWidth};
    }

    double BlendedScreeningModel::weight(const double H) const {
        return weight_and_slope(H).w;
    }

    template <typename T>
    std::vector<T> BlendedScreeningModel::blend_factors(const std::vector<T>& weak, const std::vector<T>& strong) const {
        std::vector<T> factors(weak.size());
        const T one(1.0);
        const T lnWeak(m_lnWeak);
        const T lnWidth(m_lnWidth);
        const T Hweak(m_window.weakExponent);
        for (size_t r = 0; r < weak.size(); ++r) {
            const T H = CppAD::log(weak[r]);
            const T h = CppAD::log(strong[r]);
            // s from ln H, with H floored at the window's lower edge so the logarithm stays finite on the CppAD tape
            const T H_safe = CppAD::CondExpGt(H, Hweak, H, Hweak);
            const T s_raw = (CppAD::log(H_safe) - lnWeak) / lnWidth;
            const T s = CppAD::CondExpGt(s_raw, one, one, s_raw);
            const T w = s * s * (3.0 - 2.0 * s);
            factors[r] = CppAD::exp((one - w) * H + w * h);
        }
        return factors;
    }

    std::vector<double> BlendedScreeningModel::calculateScreeningFactors(
        const reaction::ReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) const {
        std::vector<double> weak = m_weak->calculateScreeningFactors(reactions, species, Y, T9, rho);
        if (!any_blended(weak)) {
            return weak; // every reaction is below the window (the usual case in main-sequence stars)
        }
        return blend_factors<double>(weak, m_strong->calculateScreeningFactors(reactions, species, Y, T9, rho, needs_strong(weak)));
    }

    std::vector<char> BlendedScreeningModel::needs_strong(const std::vector<double>& weakFactors) const {
        const double threshold = std::exp(m_window.weakExponent);
        std::vector<char> needed(weakFactors.size());
        for (size_t r = 0; r < weakFactors.size(); ++r) {
            needed[r] = weakFactors[r] > threshold ? 1 : 0;
        }
        return needed;
    }

    bool BlendedScreeningModel::any_blended(const std::vector<double>& weakFactors) const {
        const double threshold = std::exp(m_window.weakExponent);
        for (const double f : weakFactors) {
            if (f > threshold) {
                return true;
            }
        }
        return false;
    }

    std::vector<BlendedScreeningModel::ADDouble> BlendedScreeningModel::calculateScreeningFactors(
        const reaction::ReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<ADDouble>& Y,
        const ADDouble T9,
        const ADDouble rho
    ) const {
        return blend_factors<ADDouble>(
            m_weak->calculateScreeningFactors(reactions, species, Y, T9, rho),
            m_strong->calculateScreeningFactors(reactions, species, Y, T9, rho)
        );
    }

    double BlendedScreeningModel::calculateScreeningFactor(
        const reaction::Reaction& reaction,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) const {
        const double H = std::log(m_weak->calculateScreeningFactor(reaction, species, Y, T9, rho));
        const double w = weight(H);
        if (w == 0.0) {
            return std::exp(H);
        }
        const double h = std::log(m_strong->calculateScreeningFactor(reaction, species, Y, T9, rho));
        return std::exp((1.0 - w) * H + w * h);
    }

    std::optional<std::vector<double>> BlendedScreeningModel::calculateScreeningLogDerivatives(
        const reaction::ReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) const {
        const auto dWeak = m_weak->calculateScreeningLogDerivatives(reactions, species, Y, T9, rho);
        if (!dWeak) {
            return std::nullopt;
        }
        const std::vector<double> fWeak = m_weak->calculateScreeningFactors(reactions, species, Y, T9, rho);
        if (!any_blended(fWeak)) {
            return dWeak;
        }
        const std::vector<char> needed = needs_strong(fWeak);
        const auto dStrong = m_strong->calculateScreeningLogDerivatives(reactions, species, Y, T9, rho, needed);
        if (!dStrong) {
            return std::nullopt;
        }
        const std::vector<double> fStrong = m_strong->calculateScreeningFactors(reactions, species, Y, T9, rho, needed);

        const size_t n = species.size();
        std::vector<double> dlnf(reactions.size() * n, 0.0);
        for (size_t r = 0; r < reactions.size(); ++r) {
            const double H = std::log(fWeak[r]);
            const double h = std::log(fStrong[r]);
            const auto [w, dw_dlnH] = weight_and_slope(H);
            // d ln f = (1 - w) dH + w dh + (h - H) (dw/dlnH) dH / H
            const double coupling = (dw_dlnH != 0.0 && H > 0.0) ? (h - H) * dw_dlnH / H : 0.0;
            const double* weakRow = dWeak->data() + r * n;
            const double* strongRow = dStrong->data() + r * n;
            double* row = dlnf.data() + r * n;
            for (size_t k = 0; k < n; ++k) {
                row[k] = (1.0 - w + coupling) * weakRow[k] + w * strongRow[k];
            }
        }
        return dlnf;
    }

    std::optional<std::vector<ScreeningModel::StateLogDerivatives>> BlendedScreeningModel::calculateScreeningStateLogDerivatives(
        const reaction::ReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) const {
        const auto dWeak = m_weak->calculateScreeningStateLogDerivatives(reactions, species, Y, T9, rho);
        if (!dWeak) {
            return std::nullopt;
        }
        const std::vector<double> fWeak = m_weak->calculateScreeningFactors(reactions, species, Y, T9, rho);
        if (!any_blended(fWeak)) {
            return dWeak;
        }
        const std::vector<char> needed = needs_strong(fWeak);
        const auto dStrong = m_strong->calculateScreeningStateLogDerivatives(reactions, species, Y, T9, rho, needed);
        if (!dStrong) {
            return std::nullopt;
        }
        const std::vector<double> fStrong = m_strong->calculateScreeningFactors(reactions, species, Y, T9, rho, needed);

        std::vector<StateLogDerivatives> derivatives(reactions.size());
        for (size_t r = 0; r < reactions.size(); ++r) {
            const double H = std::log(fWeak[r]);
            const double h = std::log(fStrong[r]);
            const auto [w, dw_dlnH] = weight_and_slope(H);
            const double coupling = (dw_dlnH != 0.0 && H > 0.0) ? (h - H) * dw_dlnH / H : 0.0;
            const double a = 1.0 - w + coupling;
            derivatives[r] = {
                .dlnT = a * (*dWeak)[r].dlnT + w * (*dStrong)[r].dlnT,
                .dlnRho = a * (*dWeak)[r].dlnRho + w * (*dStrong)[r].dlnRho
            };
        }
        return derivatives;
    }
}
