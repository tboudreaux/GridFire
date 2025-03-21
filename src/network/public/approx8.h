#ifndef APPROX8_H
#define APPROX8_H

#include <array>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

#include "network.h"

/**
 * @file approx8.h
 * @brief Header file for the Approx8 nuclear reaction network.
 *
 * This file contains the definitions and declarations for the Approx8 nuclear reaction network.
 * The network is based on Frank Timmes' "aprox8" and includes 8 isotopes and various nuclear reactions.
 * The rates are evaluated using a fitting function with coefficients from reaclib.jinaweb.org.
 */

/**
 * @typedef vector_type
 * @brief Alias for a vector of doubles using Boost uBLAS.
 */
typedef boost::numeric::ublas::vector< double > vector_type;

/**
 * @typedef matrix_type
 * @brief Alias for a matrix of doubles using Boost uBLAS.
 */
typedef boost::numeric::ublas::matrix< double > matrix_type;

/**
 * @typedef vec7
 * @brief Alias for a std::array of 7 doubles.
 */
typedef std::array<double,7> vec7;

namespace nnApprox8{

    using namespace boost::numeric::odeint;

    /**
     * @struct Net
     * @brief Contains constants and arrays related to the nuclear network.
     */
    struct Net{
        const static int ih1=0;
        const static int ihe3=1;
        const static int ihe4=2;
        const static int ic12=3;
        const static int in14=4;
        const static int io16=5;
        const static int ine20=6;
        const static int img24=7;

        const static int itemp=img24+1;
        const static int iden =itemp+1;
        const static int iener=iden+1;

        const static int niso=img24+1; // number of isotopes
        const static int nvar=iener+1; // number of variables

        static constexpr std::array<int,niso> aion = {
            1,
            3,
            4,
            12,
            14,
            16,
            20,
            24
        };

        static constexpr std::array<double,niso> mion = {
            1.67262164e-24, 
            5.00641157e-24, 
            6.64465545e-24, 
            1.99209977e-23,
            2.32462686e-23,
            2.65528858e-23,
            3.31891077e-23,
            3.98171594e-23
        };

    };

    /**
     * @brief Multiplies two arrays and sums the resulting elements.
     * @param a First array.
     * @param b Second array.
     * @return Sum of the product of the arrays.
     * @example
     * @code
     * vec7 a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
     * vec7 b = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
     * double result = sum_product(a, b);
     * @endcode
     */
    double sum_product( const vec7 &a, const vec7 &b);

    /**
     * @brief Returns an array of T9 terms for the nuclear reaction rate fit.
     * @param T Temperature in GigaKelvin.
     * @return Array of T9 terms.
     * @example
     * @code
     * double T = 1.5;
     * vec7 T9_array = get_T9_array(T);
     * @endcode
     */
    vec7 get_T9_array(const double &T);

    /**
     * @brief Evaluates the nuclear reaction rate given the T9 array and coefficients.
     * @param T9 Array of T9 terms.
     * @param coef Array of coefficients.
     * @return Evaluated rate.
     * @example
     * @code
     * vec7 T9 = get_T9_array(1.5);
     * vec7 coef = {1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
     * double rate = rate_fit(T9, coef);
     * @endcode
     */
    double rate_fit(const vec7 &T9, const vec7 &coef);

    /**
     * @brief Calculates the rate for the reaction p + p -> d.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double pp_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction p + d -> he3.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double dp_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction he3 + he3 -> he4 + 2p.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double he3he3_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction he3(he3,2p)he4.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double he3he4_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction he4 + he4 + he4 -> c12.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double triple_alpha_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction c12 + p -> n13.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double c12p_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction c12 + he4 -> o16.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double c12a_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction n14(p,g)o15 - o15 + p -> c12 + he4.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double n14p_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction n14(a,g)f18 assumed to go on to ne20.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double n14a_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction n15(p,a)c12 (CNO I).
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double n15pa_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction n15(p,g)o16 (CNO II).
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double n15pg_rate(const vec7 &T9);

    /**
     * @brief Calculates the fraction for the reaction n15(p,g)o16.
     * @param T9 Array of T9 terms.
     * @return Fraction of the reaction.
     */
    double n15pg_frac(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction o16(p,g)f17 then f17 -> o17(p,a)n14.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double o16p_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction o16(a,g)ne20.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double o16a_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction ne20(a,g)mg24.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double ne20a_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction c12(c12,a)ne20.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double c12c12_rate(const vec7 &T9);

    /**
     * @brief Calculates the rate for the reaction c12(o16,a)mg24.
     * @param T9 Array of T9 terms.
     * @return Rate of the reaction.
     */
    double c12o16_rate(const vec7 &T9);

    /**
     * @struct Jacobian
     * @brief Functor to calculate the Jacobian matrix for implicit solvers.
     */
    struct Jacobian {
        /**
         * @brief Calculates the Jacobian matrix.
         * @param y State vector.
         * @param J Jacobian matrix.
         * @param t Time.
         * @param dfdt Derivative of the state vector.
         */
        void operator() ( const vector_type &y, matrix_type &J, double /* t */, vector_type &dfdt );
    };

    /**
     * @struct ODE
     * @brief Functor to calculate the derivatives for the ODE solver.
     */
    struct ODE {
        /**
         * @brief Calculates the derivatives of the state vector.
         * @param y State vector.
         * @param dydt Derivative of the state vector.
         * @param t Time.
         */
        void operator() ( const vector_type &y, vector_type &dydt, double /* t */);
    };

    /**
     * @class Approx8Network
     * @brief Class for the Approx8 nuclear reaction network.
     */
    class Approx8Network : public nuclearNetwork::Network {
    public:
        /**
         * @brief Evaluates the nuclear network.
         * @param netIn Input parameters for the network.
         * @return Output results from the network.
         */
        virtual nuclearNetwork::NetOut evaluate(const nuclearNetwork::NetIn &netIn);

        /**
         * @brief Sets whether the solver should use a stiff method.
         * @param stiff Boolean indicating if a stiff method should be used.
         */
        void setStiff(bool stiff);

        /**
         * @brief Checks if the solver is using a stiff method.
         * @return Boolean indicating if a stiff method is being used.
         */
        bool isStiff() { return m_stiff; }
    private:
        vector_type m_y;
        double m_tmax;
        double m_dt0;
        bool m_stiff;

        /**
         * @brief Converts the input parameters to the internal state vector.
         * @param netIn Input parameters for the network.
         * @return Internal state vector.
         */
        vector_type convert_netIn(const nuclearNetwork::NetIn &netIn);
    };

} // namespace nnApprox8

#endif