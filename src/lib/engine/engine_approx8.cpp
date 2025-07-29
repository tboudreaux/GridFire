/* ***********************************************************************
//
//   Copyright (C) 2025 -- The 4D-STAR Collaboration
//   File Author: Emily Boudreaux
//   Last Modified: March 21, 2025
//
//   4DSSE is free software; you can use it and/or modify
//   it under the terms and restrictions the GNU General Library Public
//   License version 3 (GPLv3) as published by the Free Software Foundation.
//
//   4DSSE is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Library General Public License for more details.
//
//   You should have received a copy of the GNU Library General Public License
//   along with this software; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// *********************************************************************** */
#include <cmath>
#include <stdexcept>
#include <array>

#include <boost/numeric/odeint.hpp>

#include "fourdst/constants/const.h"
#include "fourdst/config/config.h"
#include "quill/LogMacros.h"

#include "gridfire/engine/engine_approx8.h"
#include "gridfire/network.h"

/* Nuclear reaction network in cgs units based on Frank Timmes' "approx8".
   At this time it does neither screening nor neutrino losses. It includes
   the following 8 isotopes:

   h1
   he3
   he4
   c12
   n14
   o16
   ne20
   mg24

   and the following nuclear reactions:

   ---pp chain---
   p(p,e+)d
   d(p,g)he3
   he3(he3,2p)he4

   ---CNO cycle---
   c12(p,g)n13 - n13 -> c13 + p -> n14
   n14(p,g)o15 - o15 + p -> c12 + he4
   n14(a,g)f18 - proceeds to ne20
   n15(p,a)c12 - / these two n15 reactions are  \  CNO I
   n15(p,g)o16 - \ used to calculate a fraction /  CNO II
   o16(p,g)f17 - f17 + e -> o17 and then o17 + p -> n14 + he4

   ---alpha captures---
   c12(a,g)o16 
   triple alpha
   o16(a,g)ne20
   ne20(a,g)mg24
   c12(c12,a)ne20
   c12(o16,a)mg24

At present the rates are all evaluated using a fitting function.
The coefficients to the fit are from reaclib.jinaweb.org .

*/

namespace gridfire::approx8{

	// using namespace std;
	using namespace boost::numeric::odeint;

	//helper functions
	// a function to multiply two arrays and then sum the resulting elements: sum(a*b)
	double sum_product( const vec7 &a, const vec7 &b){
		double sum=0;
		for (size_t i=0; i < a.size(); i++) {
			sum += a[i] * b[i];
		}
		return sum;
	}

	// the fit to the nuclear reaction rates is of the form:
	// exp( a0 + a1/T9 + a2/T9^(1/3) + a3*T9^(1/3) + a4*T9 + a5*T9^(5/3) + log(T9) )
	// this function returns an array of the T9 terms in that order, where T9 is the temperatures in GigaKelvin
	vec7 get_T9_array(const double &T) {
		vec7 arr;
		const double T9=1e-9*T;
		const double T913=pow(T9,1./3.);

		arr[0]=1;
		arr[1]=1/T9;
		arr[2]=1/T913;
		arr[3]=T913;
		arr[4]=T9;
		arr[5]=pow(T9,5./3.);
		arr[6]=log(T9);

		return arr;
	}

	// this function uses the two preceding functions to evaluate the rate given the T9 array and the coefficients
	double rate_fit(const vec7 &T9, const vec7 &coef){
		return exp(sum_product(T9,coef));
	}

	// p + p -> d; this, like some of the other rates, this is a composite of multiple fits
	double pp_rate(const vec7 &T9) {
		constexpr vec7 a1 = {-34.78630,              0,-3.511930, 3.100860, -0.1983140,    1.262510e-2, -1.025170};
		constexpr vec7 a2 = { -4.364990e+1,-2.460640e-3,-2.750700,-4.248770e-1,1.598700e-2,-6.908750e-4,-2.076250e-1};
		return rate_fit(T9,a1) + rate_fit(T9,a2);
	}

	// p + d -> he3
	double dp_rate(const vec7 &T9) {
		constexpr vec7 a1 = {7.528980, 0, -3.720800, 0.8717820, 0, 0,-0.6666670};
		constexpr vec7 a2 = {8.935250, 0, -3.720800, 0.1986540, 0, 0, 0.3333330};
		return rate_fit(T9,a1) + rate_fit(T9,a2);
	}

	// he3 + he3 -> he4 + 2p
	double he3he3_rate(const vec7 &T9){
		constexpr vec7 a = {2.477880e+01,0,-12.27700,-0.1036990,-6.499670e-02,1.681910e-02,-6.666670e-01};
		return rate_fit(T9,a);
	}

	// he3(he3,2p)he4 ** (missing both coefficients but have a reaction)
	double he3he4_rate(const vec7 &T9){
		constexpr vec7 a1 = {1.560990e+01,0.000000e+00,-1.282710e+01,-3.082250e-02,-6.546850e-01,8.963310e-02,-6.666670e-01};
		constexpr vec7 a2 = {1.770750e+01,0.000000e+00,-1.282710e+01,-3.812600e+00,9.422850e-02,-3.010180e-03,1.333330e+00};
		return rate_fit(T9,a1) + rate_fit(T9,a2);
	}

	// he4 + he4 + he4 -> c12 ** (missing middle coefficient but have other two)
	double triple_alpha_rate(const vec7 &T9){
		constexpr vec7 a1 = {-9.710520e-01,0.000000e+00,-3.706000e+01,2.934930e+01,-1.155070e+02,-1.000000e+01,-1.333330e+00};
		constexpr vec7 a2 = {-1.178840e+01,-1.024460e+00,-2.357000e+01,2.048860e+01,-1.298820e+01,-2.000000e+01,-2.166670e+00};
		constexpr vec7 a3 = {-2.435050e+01,-4.126560e+00,-1.349000e+01,2.142590e+01,-1.347690e+00,8.798160e-02,-1.316530e+01};
		return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3);
	}

	// c12 + p -> n13
	double c12p_rate(const vec7 &T9){
		constexpr vec7 a1={1.714820e+01,0.000000e+00,-1.369200e+01,-2.308810e-01,4.443620e+00,-3.158980e+00,-6.666670e-01};
		constexpr vec7 a2={1.754280e+01,-3.778490e+00,-5.107350e+00,-2.241110e+00,1.488830e-01,0.000000e+00,-1.500000e+00};
		return rate_fit(T9,a1) + rate_fit(T9,a2);
	}

	// c12 + he4 -> o16 ** (missing first coefficient but have the second)
	double c12a_rate(const vec7 &T9){
		constexpr vec7 a1={6.965260e+01,-1.392540e+00,5.891280e+01,-1.482730e+02,9.083240e+00,-5.410410e-01,7.035540e+01};
		constexpr vec7 a2={2.546340e+02,-1.840970e+00,1.034110e+02,-4.205670e+02,6.408740e+01,-1.246240e+01,1.373030e+02};
		return rate_fit(T9,a1) + rate_fit(T9,a2);
	}

	// n14(p,g)o15 - o15 + p -> c12 + he4
	double n14p_rate(const vec7 &T9){
		constexpr vec7 a1={1.701000e+01,0.000000e+00,-1.519300e+01,-1.619540e-01,-7.521230e+00,-9.875650e-01,-6.666670e-01};
		constexpr vec7 a2={2.011690e+01,0.000000e+00,-1.519300e+01,-4.639750e+00,9.734580e+00,-9.550510e+00,3.333330e-01};
		constexpr vec7 a3={7.654440e+00,-2.998000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
		constexpr vec7 a4={6.735780e+00,-4.891000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.820000e-02};
		return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3) + rate_fit(T9,a4);
	}

	// n14(a,g)f18 assumed to go on to ne20
	double n14a_rate(const vec7 &T9){
		constexpr vec7 a1={2.153390e+01,0.000000e+00,-3.625040e+01,0.000000e+00,0.000000e+00,-5.000000e+00,-6.666670e-01};
		constexpr vec7 a2={1.968380e-01,-5.160340e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
		constexpr vec7 a3={1.389950e+01,-1.096560e+01,-5.622700e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
		return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3);
	}

	// n15(p,a)c12 (CNO I)
	double n15pa_rate(const vec7 &T9){
		constexpr vec7 a1 = {2.747640e+01,0.000000e+00,-1.525300e+01,1.593180e+00,2.447900e+00,-2.197080e+00,-6.666670e-01};
		constexpr vec7 a2 = {-4.873470e+00,-2.021170e+00,0.000000e+00,3.084970e+01,-8.504330e+00,-1.544260e+00,-1.500000e+00};
		constexpr vec7 a3 = {2.089720e+01,-7.406000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
		constexpr vec7 a4 = {-6.575220e+00,-1.163800e+00,0.000000e+00,2.271050e+01,-2.907070e+00,2.057540e-01,-1.500000e+00};
		return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3) + rate_fit(T9,a4);
	}

	// n15(p,g)o16 (CNO II)
	double n15pg_rate(const vec7 &T9){
		constexpr vec7 a1 = {2.001760e+01,0.000000e+00,-1.524000e+01,3.349260e-01,4.590880e+00,-4.784680e+00,-6.666670e-01};
		constexpr vec7 a2 = {6.590560e+00,-2.923150e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
		constexpr vec7 a3 = {1.454440e+01,-1.022950e+01,0.000000e+00,0.000000e+00,4.590370e-02,0.000000e+00,-1.500000e+00};
		return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3);
	}

	double n15pg_frac(const vec7 &T9){
		const double f1=n15pg_rate(T9);
		const double f2=n15pa_rate(T9);
		return f1/(f1+f2);
	}

	// o16(p,g)f17 then f17 -> o17(p,a)n14
	double o16p_rate(const vec7 &T9){
		constexpr vec7 a={1.909040e+01,0.000000e+00,-1.669600e+01,-1.162520e+00,2.677030e-01,-3.384110e-02,-6.666670e-01};
		return rate_fit(T9,a);
	}

		// o16(a,g)ne20
	double o16a_rate(const vec7 &T9){
		constexpr vec7 a1={2.390300e+01,0.000000e+00,-3.972620e+01,-2.107990e-01,4.428790e-01,-7.977530e-02,-6.666670e-01};
		constexpr vec7 a2={3.885710e+00,-1.035850e+01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
		constexpr vec7 a3={9.508480e+00,-1.276430e+01,0.000000e+00,-3.659250e+00,7.142240e-01,-1.075080e-03,-1.500000e+00};
		return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3);
	}

		// ne20(a,g)mg24
	double ne20a_rate(const vec7 &T9){
		constexpr vec7 a1={2.450580e+01,0.000000e+00,-4.625250e+01,5.589010e+00,7.618430e+00,-3.683000e+00,-6.666670e-01};
		constexpr vec7 a2={-3.870550e+01,-2.506050e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
		constexpr vec7 a3={1.983070e+00,-9.220260e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
		constexpr vec7 a4={-8.798270e+00,-1.278090e+01,0.000000e+00,1.692290e+01,-2.573250e+00,2.089970e-01,-1.500000e+00};
		return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3) + rate_fit(T9,a4);
	}

		// c12(c12,a)ne20
	double c12c12_rate(const vec7 &T9){
		constexpr vec7 a={6.128630e+01,0.000000e+00,-8.416500e+01,-1.566270e+00,-7.360840e-02,-7.279700e-02,-6.666670e-01};
		return rate_fit(T9,a);
	}

		// c12(o16,a)mg24
	double c12o16_rate(const vec7 &T9){
		constexpr vec7 a={4.853410e+01,3.720400e-01,-1.334130e+02,5.015720e+01,-3.159870e+00,1.782510e-02,-2.370270e+01};
		return rate_fit(T9,a);
	}


	// for Boost ODE solvers either a struct or a class is required

	// a Jacobian matrix for implicit solvers

    void Jacobian::operator() ( const vector_type &y, matrix_type &J, double /* t */, vector_type &dfdt ) const {
        fourdst::constant::Constants& constants = fourdst::constant::Constants::getInstance();
        const double avo = constants.get("N_a").value;
        const double clight = constants.get("c").value;
        // EOS
        const vec7 T9=get_T9_array(y[Approx8Net::iTemp]);

        // evaluate rates once per call
        const double rpp=pp_rate(T9);
        const double r33=he3he3_rate(T9);
        const double r34=he3he4_rate(T9);
        const double r3a=triple_alpha_rate(T9);
        const double rc12p=c12p_rate(T9);
        const double rc12a=c12a_rate(T9);
        const double rn14p=n14p_rate(T9);
        const double rn14a=n14a_rate(T9);
        const double ro16p=o16p_rate(T9);
        const double ro16a=o16a_rate(T9);
        const double rne20a=ne20a_rate(T9);
        const double r1212=c12c12_rate(T9);
        const double r1216=c12o16_rate(T9);

        const double pFrac=n15pg_frac(T9);
        const double aFrac=1-pFrac;

        const double   yh1 = y[Approx8Net::ih1];
        const double  yhe3 = y[Approx8Net::ihe3];
        const double  yhe4 = y[Approx8Net::ihe4];
        const double  yc12 = y[Approx8Net::ic12];
        const double  yn14 = y[Approx8Net::in14];
        const double  yo16 = y[Approx8Net::io16];
        const double yne20 = y[Approx8Net::ine20];

        // zero all elements to begin
        for (int i=0; i < Approx8Net::nVar; i++) {
            for (int j=0; j < Approx8Net::nVar; j++) {
				J(i,j)=0.0;
			}
		}

        // h1 jacobian elements
        J(Approx8Net::ih1,Approx8Net::ih1)  = -3*yh1*rpp - 2*yc12*rc12p -2*yn14*rn14p -2*yo16*ro16p;
        J(Approx8Net::ih1,Approx8Net::ihe3) =  2*yhe3*r33 - yhe4*r34;
        J(Approx8Net::ih1,Approx8Net::ihe4) = -yhe3*r34;
        J(Approx8Net::ih1,Approx8Net::ic12) = -2*yh1*rc12p;
        J(Approx8Net::ih1,Approx8Net::in14) = -2*yh1*rn14p;
        J(Approx8Net::ih1,Approx8Net::io16) = -2*yh1*ro16p;

        // he3 jacobian elements
        J(Approx8Net::ihe3,Approx8Net::ih1)   =  yh1*rpp;
        J(Approx8Net::ihe3,Approx8Net::ihe3)  = -2*yhe3*r33 - yhe4*r34;
        J(Approx8Net::ihe3,Approx8Net::ihe4)  = -yhe3*r34;

        // he4 jacobian elements
        J(Approx8Net::ihe4,Approx8Net::ih1)   =  yn14*aFrac*rn14p + yo16*ro16p;
        J(Approx8Net::ihe4,Approx8Net::ihe3)  =  yhe3*r33 - yhe4*r34;
        J(Approx8Net::ihe4,Approx8Net::ihe4)  =  yhe3*r34 - 1.5*yhe4*yhe4*r3a - yc12*rc12a - 1.5*yn14*rn14a - yo16*ro16a - yne20*rne20a;
        J(Approx8Net::ihe4,Approx8Net::ic12)  = -yhe4*rc12a + yc12*r1212 + yo16*r1216;
        J(Approx8Net::ihe4,Approx8Net::in14)  =  yh1*aFrac*rn14p - 1.5*yhe4*rn14a;
        J(Approx8Net::ihe4,Approx8Net::io16)  =  yh1*ro16p - yhe4*ro16a + yc12*r1216;
        J(Approx8Net::ihe4,Approx8Net::ine20) = -yhe4*rne20a;

        // c12 jacobian elements
        J(Approx8Net::ic12,Approx8Net::ih1)   = -yc12*rc12p + yn14*aFrac*rn14p;
        J(Approx8Net::ic12,Approx8Net::ihe4)  = 0.5*yhe4*yhe4*r3a - yhe4*rc12a;
        J(Approx8Net::ic12,Approx8Net::ic12)  = -yh1*rc12p - yhe4*rc12a - yo16*r1216 - 2*yc12*r1212;
        J(Approx8Net::ic12,Approx8Net::in14)  =  yh1*yn14*aFrac*rn14p;
        J(Approx8Net::ic12,Approx8Net::io16)  = -yc12*r1216;

        // n14 jacobian elements
        J(Approx8Net::in14,Approx8Net::ih1)   =  yc12*rc12p - yn14*rn14p + yo16*ro16p;
        J(Approx8Net::in14,Approx8Net::ihe4)  = -yn14*rn14a;
        J(Approx8Net::in14,Approx8Net::ic12)  =  yh1*rc12p;
        J(Approx8Net::in14,Approx8Net::in14)  = -yh1*rn14p - yhe4*rn14a;
        J(Approx8Net::in14,Approx8Net::io16)  =  yo16*ro16p;

        // o16 jacobian elements
        J(Approx8Net::io16,Approx8Net::ih1)   = yn14*pFrac*rn14p - yo16*ro16p;
        J(Approx8Net::io16,Approx8Net::ihe4)  = yc12*rc12a - yo16*ro16a;
        J(Approx8Net::io16,Approx8Net::ic12)  = yhe4*rc12a - yo16*r1216;
        J(Approx8Net::io16,Approx8Net::in14)  = yh1*pFrac*rn14p;
        J(Approx8Net::io16,Approx8Net::io16)  = yh1*ro16p - yc12*r1216 -yhe4*ro16a;

        // ne20 jacobian elements
        J(Approx8Net::ine20,Approx8Net::ihe4)  =  yn14*rn14a + yo16*ro16a - yne20*rne20a;
        J(Approx8Net::ine20,Approx8Net::ic12)  =  yc12*r1212;
        J(Approx8Net::ine20,Approx8Net::in14)  =  yhe4*rn14a;
        J(Approx8Net::ine20,Approx8Net::io16)  =  yo16*ro16a;
        J(Approx8Net::ine20,Approx8Net::ine20) = -yhe4*rne20a;

        // mg24 jacobian elements
        J(Approx8Net::img24,Approx8Net::ihe4) = yne20*rne20a;
        J(Approx8Net::img24,Approx8Net::ic12) = yo16*r1216;
        J(Approx8Net::img24,Approx8Net::io16) = yc12*r1216;
        J(Approx8Net::img24,Approx8Net::ine20) = yhe4*rne20a;

        // energy accounting
        for (int j=0; j<Approx8Net::nIso; j++) {
            for (int i=0; i<Approx8Net::nIso; i++) {
				J(Approx8Net::iEnergy,j) += J(i,j)*Approx8Net::mIon[i];
			}
            J(Approx8Net::iEnergy,j) *= -avo*clight*clight;
        }
    }

	void ODE::operator() ( const vector_type &y, vector_type &dydt, double /* t */) const {
		const fourdst::constant::Constants& constants = fourdst::constant::Constants::getInstance();
		const double avo = constants.get("N_a").value;
		const double clight = constants.get("c").value;

		// EOS
		const double T = y[Approx8Net::iTemp];
		const double den = y[Approx8Net::iDensity];
		const vec7 T9=get_T9_array(T);

		// rates
		const double rpp=den*pp_rate(T9);
		const double r33=den*he3he3_rate(T9);
		const double r34=den*he3he4_rate(T9);
		const double r3a=den*den*triple_alpha_rate(T9);
		const double rc12p=den*c12p_rate(T9);
		const double rc12a=den*c12a_rate(T9);
		const double rn14p=den*n14p_rate(T9);
		const double rn14a=n14a_rate(T9);
		const double ro16p=den*o16p_rate(T9);
		const double ro16a=den*o16a_rate(T9);
		const double rne20a=den*ne20a_rate(T9);
		const double r1212=den*c12c12_rate(T9);
		const double r1216=den*c12o16_rate(T9);

		const double pFrac=n15pg_frac(T9);
		const double aFrac=1-pFrac;

		const double   yh1 = y[Approx8Net::  ih1];
		const double  yhe3 = y[Approx8Net:: ihe3];
		const double  yhe4 = y[Approx8Net:: ihe4];
		const double  yc12 = y[Approx8Net:: ic12];
		const double  yn14 = y[Approx8Net:: in14];
		const double  yo16 = y[Approx8Net:: io16];
		const double yne20 = y[Approx8Net::ine20];

		dydt[Approx8Net::ih1]  = -1.5*yh1*yh1*rpp;
		dydt[Approx8Net::ih1] +=  yhe3*yhe3*r33;
		dydt[Approx8Net::ih1] += -yhe3*yhe4*r34;
		dydt[Approx8Net::ih1] += -2*yh1*yc12*rc12p;
		dydt[Approx8Net::ih1] += -2*yh1*yn14*rn14p;
		dydt[Approx8Net::ih1] += -2*yh1*yo16*ro16p;

		dydt[Approx8Net::ihe3]  = 0.5*yh1*yh1*rpp;
		dydt[Approx8Net::ihe3] += -yhe3*yhe3*r33;
		dydt[Approx8Net::ihe3] += -yhe3*yhe4*r34;

		dydt[Approx8Net::ihe4]  = 0.5*yhe3*yhe3*r33;
		dydt[Approx8Net::ihe4] +=  yhe3*yhe4*r34;
		dydt[Approx8Net::ihe4] += -yhe4*yc12*rc12a;
		dydt[Approx8Net::ihe4] +=  yh1*yn14*aFrac*rn14p;
		dydt[Approx8Net::ihe4] +=  yh1*yo16*ro16p;
		dydt[Approx8Net::ihe4] += -0.5*yhe4*yhe4*yhe4*r3a;
		dydt[Approx8Net::ihe4] += -yhe4*yo16*ro16a;
		dydt[Approx8Net::ihe4] += 0.5*yc12*yc12*r1212;
		dydt[Approx8Net::ihe4] +=  yc12*yo16*r1216;
		dydt[Approx8Net::ihe4] += -yhe4*yne20*rne20a;

		dydt[Approx8Net::ic12]  = (1./6.)*yhe4*yhe4*yhe4*r3a;
		dydt[Approx8Net::ic12] += -yhe4*yc12*rc12a;
		dydt[Approx8Net::ic12] += -yh1*yc12*rc12p;
		dydt[Approx8Net::ic12] +=  yh1*yn14*aFrac*rn14p;
		dydt[Approx8Net::ic12] += -yc12*yc12*r1212;
		dydt[Approx8Net::ic12] += -yc12*yo16*r1216;

		dydt[Approx8Net::in14]  =  yh1*yc12*rc12p;
		dydt[Approx8Net::in14] += -yh1*yn14*rn14p;
		dydt[Approx8Net::in14] +=  yh1*yo16*ro16p;
		dydt[Approx8Net::in14] += -yhe4*yn14*rn14a;

		dydt[Approx8Net::io16]  =  yhe4*yc12*rc12a;
		dydt[Approx8Net::io16] +=  yh1*yn14*pFrac*rn14p;
		dydt[Approx8Net::io16] += -yh1*yo16*ro16p;
		dydt[Approx8Net::io16] += -yc12*yo16*r1216;
		dydt[Approx8Net::io16] += -yhe4*yo16*ro16a;

		dydt[Approx8Net::ine20]  = 0.5*yc12*yc12*r1212;
		dydt[Approx8Net::ine20] +=  yhe4*yn14*rn14a;
		dydt[Approx8Net::ine20] +=  yhe4*yo16*ro16a;
		dydt[Approx8Net::ine20] += -yhe4*yne20*rne20a;

		dydt[Approx8Net::img24]  =  yc12*yo16*r1216;
		dydt[Approx8Net::img24] +=  yhe4*yne20*rne20a;

		dydt[Approx8Net::iTemp] = 0.;
		dydt[Approx8Net::iDensity]  = 0.;

		// energy accounting
		double eNuc = 0.;
		for (int i=0; i<Approx8Net::nIso; i++) {
			eNuc += Approx8Net::mIon[i]*dydt[i];
		}
		dydt[Approx8Net::iEnergy] = -eNuc*avo*clight*clight;
	}

	Approx8Network::Approx8Network() : Network(APPROX8) {}

	NetOut Approx8Network::evaluate(const NetIn &netIn) {
		m_y = convert_netIn(netIn);
		m_tMax = netIn.tMax;
		m_dt0 = netIn.dt0;

		const double stiff_abs_tol = m_config.get<double>("Network:Approx8:Stiff:AbsTol", 1.0e-6);
		const double stiff_rel_tol = m_config.get<double>("Network:Approx8:Stiff:RelTol", 1.0e-6);
		const double nonstiff_abs_tol = m_config.get<double>("Network:Approx8:NonStiff:AbsTol", 1.0e-6);
		const double nonstiff_rel_tol = m_config.get<double>("Network:Approx8:NonStiff:RelTol", 1.0e-6);

		int num_steps = -1;

		if (m_stiff) {
			LOG_DEBUG(m_logger, "Using stiff solver for Approx8Network");
			num_steps = integrate_const(
				make_dense_output<rosenbrock4<double>>(stiff_abs_tol, stiff_rel_tol),
				std::make_pair(ODE(), Jacobian()),
				m_y,
				0.0,
				m_tMax,
				m_dt0
			);

		} else {
			LOG_DEBUG(m_logger, "Using non stiff solver for Approx8Network");
			num_steps = integrate_const (
				make_dense_output<runge_kutta_dopri5<vector_type>>(nonstiff_abs_tol, nonstiff_rel_tol),
				ODE(),
				m_y,
				0.0,
				m_tMax,
				m_dt0
			);
		}

		double ySum = 0.0;
		for (int i = 0; i < Approx8Net::nIso; i++) {
			m_y[i] *= Approx8Net::aIon[i];
			ySum += m_y[i];
		}
		for (int i = 0; i < Approx8Net::nIso; i++) {
			m_y[i] /= ySum;
		}

		NetOut netOut;
		std::vector<double> outComposition;
		outComposition.reserve(Approx8Net::nVar);

		for (int i = 0; i < Approx8Net::nIso; i++) {
			outComposition.push_back(m_y[i]);
		}
		netOut.energy = m_y[Approx8Net::iEnergy];
		netOut.num_steps = num_steps;

		const std::vector<std::string> symbols = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
		netOut.composition = fourdst::composition::Composition(symbols, outComposition);

		return netOut;
	}

	void Approx8Network::setStiff(bool stiff) {
		m_stiff = stiff;
	}

	vector_type Approx8Network::convert_netIn(const NetIn &netIn) {
		vector_type y(Approx8Net::nVar, 0.0);
		y[Approx8Net::ih1] = netIn.composition.getNumberFraction("H-1");
		y[Approx8Net::ihe3] = netIn.composition.getNumberFraction("He-3");
		y[Approx8Net::ihe4] = netIn.composition.getNumberFraction("He-4");
		y[Approx8Net::ic12] = netIn.composition.getNumberFraction("C-12");
		y[Approx8Net::in14] = netIn.composition.getNumberFraction("N-14");
		y[Approx8Net::io16] = netIn.composition.getNumberFraction("O-16");
		y[Approx8Net::ine20] = netIn.composition.getNumberFraction("Ne-20");
		y[Approx8Net::img24] = netIn.composition.getNumberFraction("Mg-24");
		y[Approx8Net::iTemp] = netIn.temperature;
		y[Approx8Net::iDensity] = netIn.density;
		y[Approx8Net::iEnergy] = netIn.energy;

		return y;
	}
};


// main program 

