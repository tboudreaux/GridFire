#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <stdexcept>
#include <array>
#include <numbers>
#include <utility>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

#include "const.h"
//#include "probe.h"
//#include "config.h"
//#include "quill/LogMacros.h"

/* Nuclear reaction network in cgs units based on Frank Timmes' "aprox8".
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


// using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

// these types are required by the Rosenbrock implicit solver
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

// this array is used only in the nuclear reaction rate evaluations
typedef std::array<double,7> vec;

// only need a couple of constants
Constants& constants = Constants::getInstance();
const double avo = constants.get("N_a").value;
const double clight = constants.get("c").value;

// identify the isotopes used in the network.
const int ih1=0;
const int ihe3=1;
const int ihe4=2;
const int ic12=3;
const int in14=4;
const int io16=5;
const int ine20=6;
const int img24=7;

// physical variables; this routine currently does not need to call EOS
// since the temperature and density are assumed constant during the burn
const int itemp=img24+1;
const int iden =itemp+1;
const int iener=iden+1;

// number of isotopes and number of variables
const int niso=img24+1; // number of isotopes
const int nvar=iener+1; // number of variables

// atomic stuff
std::array<int,niso> aion = {1,3,4,12,14,16,20,24};
//std::array<int,niso> zion = {1,2,2, 6, 7, 8,10,12};
//std::array<double,niso> bion = {0,7.71819,28.29603,92.16294,104.65998,127.62093,160.64788,198.25790};
//nion = aion - zion #neutrons
//mion = nion*mn + zion*mp - bion*mev2gr #mass
std::array<double,niso> mion = {1.67262164e-24, 5.00641157e-24, 6.64465545e-24, 1.99209977e-23,
				2.32462686e-23, 2.65528858e-23, 3.31891077e-23, 3.98171594e-23};

//helper functions
// a function to multilpy two arrays and then sum the resulting elements: sum(a*b)
double sum_product( const vec &a, const vec &b){
    if (a.size() != b.size()){
	throw std::runtime_error("Error: array size mismatch in sum_product");
    }

    double sum=0;
    for (size_t i=0; i < a.size(); i++) {
	sum += a[i] * b[i];
    }
    return sum;
}

// the fit to the nuclear reaction rates is of the form:
// exp( a0 + a1/T9 + a2/T9^(1/3) + a3*T9^(1/3) + a4*T9 + a5*T9^(5/3) + log(T9) )
// this function returns an array of the T9 terms in that order, where T9 is the temperatures in GigaKelvin
vec get_T9_array(const double &T) {
    vec arr;
    double T9=1e-9*T;
    double T913=pow(T9,1./3.);
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
double rate_fit(const vec &T9, const vec &coef){return exp(sum_product(T9,coef));}

// p + p -> d; this, like some of the other rates, this is a composite of multiple fits
double pp_rate(const vec &T9) {
    vec a1 = {-34.78630,              0,-3.511930, 3.100860, -0.1983140,    1.262510e-2, -1.025170};
    vec a2 = { -4.364990e+1,-2.460640e-3,-2.750700,-4.248770e-1,1.598700e-2,-6.908750e-4,-2.076250e-1};
    return rate_fit(T9,a1) + rate_fit(T9,a2);
}

// p + d -> he3
double dp_rate(const vec &T9) {
    vec a1 = {7.528980, 0, -3.720800, 0.8717820, 0, 0,-0.6666670};
    vec a2 = {8.935250, 0, -3.720800, 0.1986540, 0, 0, 0.3333330};
    return rate_fit(T9,a1) + rate_fit(T9,a2);
}

// he3 + he3 -> he4 + 2p
double he3he3_rate(const vec &T9){
    vec a = {2.477880e+01,0,-12.27700,-0.1036990,-6.499670e-02,1.681910e-02,-6.666670e-01};
    return rate_fit(T9,a);
}

// he3(he3,2p)he4
double he3he4_rate(const vec &T9){
    vec a1 = {1.560990e+01,0.000000e+00,-1.282710e+01,-3.082250e-02,-6.546850e-01,8.963310e-02,-6.666670e-01};
    vec a2 = {1.770750e+01,0.000000e+00,-1.282710e+01,-3.812600e+00,9.422850e-02,-3.010180e-03,1.333330e+00};
    return rate_fit(T9,a1) + rate_fit(T9,a2);
}

// he4 + he4 + he4 -> c12
double triple_alpha_rate(const vec &T9){
    vec a1 = {-9.710520e-01,0.000000e+00,-3.706000e+01,2.934930e+01,-1.155070e+02,-1.000000e+01,-1.333330e+00};
    vec a2 = {-1.178840e+01,-1.024460e+00,-2.357000e+01,2.048860e+01,-1.298820e+01,-2.000000e+01,-2.166670e+00};
    vec a3 = {-2.435050e+01,-4.126560e+00,-1.349000e+01,2.142590e+01,-1.347690e+00,8.798160e-02,-1.316530e+01};
    return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3);
}

// c12 + p -> n13
double c12p_rate(const vec &T9){
    vec a1={1.714820e+01,0.000000e+00,-1.369200e+01,-2.308810e-01,4.443620e+00,-3.158980e+00,-6.666670e-01};
    vec a2={1.754280e+01,-3.778490e+00,-5.107350e+00,-2.241110e+00,1.488830e-01,0.000000e+00,-1.500000e+00};
    return rate_fit(T9,a1) + rate_fit(T9,a2);
}

// c12 + he4 -> o16
double c12a_rate(const vec &T9){
    vec a1={6.965260e+01,-1.392540e+00,5.891280e+01,-1.482730e+02,9.083240e+00,-5.410410e-01,7.035540e+01};
    vec a2={2.546340e+02,-1.840970e+00,1.034110e+02,-4.205670e+02,6.408740e+01,-1.246240e+01,1.373030e+02};
    return rate_fit(T9,a1) + rate_fit(T9,a2);
}

// n14(p,g)o15 - o15 + p -> c12 + he4
double n14p_rate(const vec &T9){
    vec a1={1.701000e+01,0.000000e+00,-1.519300e+01,-1.619540e-01,-7.521230e+00,-9.875650e-01,-6.666670e-01};
    vec a2={2.011690e+01,0.000000e+00,-1.519300e+01,-4.639750e+00,9.734580e+00,-9.550510e+00,3.333330e-01};
    vec a3={7.654440e+00,-2.998000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    vec a4={6.735780e+00,-4.891000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.820000e-02};
    return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3) + rate_fit(T9,a4);
}

// n14(a,g)f18 assumed to go on to ne20
double n14a_rate(const vec &T9){
    vec a1={2.153390e+01,0.000000e+00,-3.625040e+01,0.000000e+00,0.000000e+00,-5.000000e+00,-6.666670e-01};	
    vec a2={1.968380e-01,-5.160340e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};	
    vec a3={1.389950e+01,-1.096560e+01,-5.622700e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3);
}

// n15(p,a)c12 (CNO I)
double n15pa_rate(const vec &T9){
    vec a1 = {2.747640e+01,0.000000e+00,-1.525300e+01,1.593180e+00,2.447900e+00,-2.197080e+00,-6.666670e-01};
    vec a2 = {-4.873470e+00,-2.021170e+00,0.000000e+00,3.084970e+01,-8.504330e+00,-1.544260e+00,-1.500000e+00};
    vec a3 = {2.089720e+01,-7.406000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    vec a4 = {-6.575220e+00,-1.163800e+00,0.000000e+00,2.271050e+01,-2.907070e+00,2.057540e-01,-1.500000e+00};	
    return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3) + rate_fit(T9,a4);
}

// n15(p,g)o16 (CNO II)
double n15pg_rate(const vec &T9){
    vec a1 = {2.001760e+01,0.000000e+00,-1.524000e+01,3.349260e-01,4.590880e+00,-4.784680e+00,-6.666670e-01};
    vec a2 = {6.590560e+00,-2.923150e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    vec a3 = {1.454440e+01,-1.022950e+01,0.000000e+00,0.000000e+00,4.590370e-02,0.000000e+00,-1.500000e+00};
    return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3);
}

double n15pg_frac(const vec &T9){
    double f1=n15pg_rate(T9);
    double f2=n15pa_rate(T9);
    return f1/(f1+f2);
}

// o16(p,g)f17 then f17 -> o17(p,a)n14
double o16p_rate(const vec &T9){
    vec a={1.909040e+01,0.000000e+00,-1.669600e+01,-1.162520e+00,2.677030e-01,-3.384110e-02,-6.666670e-01};
    return rate_fit(T9,a);
}

    // o16(a,g)ne20
double o16a_rate(const vec &T9){
    vec a1={2.390300e+01,0.000000e+00,-3.972620e+01,-2.107990e-01,4.428790e-01,-7.977530e-02,-6.666670e-01};	
    vec a2={3.885710e+00,-1.035850e+01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    vec a3={9.508480e+00,-1.276430e+01,0.000000e+00,-3.659250e+00,7.142240e-01,-1.075080e-03,-1.500000e+00};
    return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3);
}

    // ne20(a,g)mg24
double ne20a_rate(const vec &T9){
    vec a1={2.450580e+01,0.000000e+00,-4.625250e+01,5.589010e+00,7.618430e+00,-3.683000e+00,-6.666670e-01};
    vec a2={-3.870550e+01,-2.506050e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    vec a3={1.983070e+00,-9.220260e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1.500000e+00};
    vec a4={-8.798270e+00,-1.278090e+01,0.000000e+00,1.692290e+01,-2.573250e+00,2.089970e-01,-1.500000e+00};
    return rate_fit(T9,a1) + rate_fit(T9,a2) + rate_fit(T9,a3) + rate_fit(T9,a4);
}

    // c12(c12,a)ne20
double c12c12_rate(const vec &T9){
    vec a={6.128630e+01,0.000000e+00,-8.416500e+01,-1.566270e+00,-7.360840e-02,-7.279700e-02,-6.666670e-01};
    return rate_fit(T9,a);
}

    // c12(o16,a)mg24
double c12o16_rate(const vec &T9){
    vec a={4.853410e+01,3.720400e-01,-1.334130e+02,5.015720e+01,-3.159870e+00,1.782510e-02,-2.370270e+01};
    return rate_fit(T9,a);
}


// for Boost ODE solvers either a struct or a class is required

// a Jacobian matrix for implicit solvers
struct Jacobian {

    void operator() ( const vector_type &y, matrix_type &J, double /* t */, vector_type &dfdt )
    {
	// EOS
	vec T9=get_T9_array(y[itemp]);

	// evaluate rates once per call
	double rpp=pp_rate(T9);
	double r33=he3he3_rate(T9);
	double r34=he3he4_rate(T9);
	double r3a=triple_alpha_rate(T9);
	double rc12p=c12p_rate(T9);
	double rc12a=c12a_rate(T9);
	double rn14p=n14p_rate(T9);
	double rn14a=n14a_rate(T9);
	double ro16p=o16p_rate(T9);
	double ro16a=o16a_rate(T9);
	double rne20a=ne20a_rate(T9);
	double r1212=c12c12_rate(T9);
	double r1216=c12o16_rate(T9);

	double pfrac=n15pg_frac(T9);
	double afrac=1-pfrac;

	double   yh1 = y[  ih1];
	double  yhe3 = y[ ihe3];
	double  yhe4 = y[ ihe4];
	double  yc12 = y[ ic12];
	double  yn14 = y[ in14];
	double  yo16 = y[ io16];
	double yne20 = y[ine20];

	// zero all elements to begin
	for (int i=0; i < nvar; i++) {
	    for (int j=0; j < nvar; j++) { J(i,j)=0.0; }}
	
	// h1 jacobian elements
	J(ih1,ih1)  = -3*yh1*rpp - 2*yc12*rc12p -2*yn14*rn14p -2*yo16*ro16p;
	J(ih1,ihe3) =  2*yhe3*r33 - yhe4*r34;
	J(ih1,ihe4) = -yhe3*r34;
	J(ih1,ic12) = -2*yh1*rc12p;
	J(ih1,in14) = -2*yh1*rn14p;
	J(ih1,io16) = -2*yh1*ro16p;

	// he3 jacobian elements
	J(ihe3,ih1)   =  yh1*rpp;
	J(ihe3,ihe3)  = -2*yhe3*r33 - yhe4*r34;
	J(ihe3,ihe4)  = -yhe3*r34;
	
	// he4 jacobian elements
	J(ihe4,ih1)   =  yn14*afrac*rn14p + yo16*ro16p;
	J(ihe4,ihe3)  =  yhe3*r33 - yhe4*r34;
	J(ihe4,ihe4)  =  yhe3*r34 - 1.5*yhe4*yhe4*r3a - yc12*rc12a - 1.5*yn14*rn14a - yo16*ro16a - yne20*rne20a;
	J(ihe4,ic12)  = -yhe4*rc12a + yc12*r1212 + yo16*r1216;
	J(ihe4,in14)  =  yh1*afrac*rn14p - 1.5*yhe4*rn14a;
	J(ihe4,io16)  =  yh1*ro16p - yhe4*ro16a + yc12*r1216;
	J(ihe4,ine20) = -yhe4*rne20a;

	// c12 jacobian elements
	J(ic12,ih1)   = -yc12*rc12p + yn14*afrac*rn14p;
	J(ic12,ihe4)  = 0.5*yhe4*yhe4*r3a - yhe4*rc12a;
	J(ic12,ic12)  = -yh1*rc12p - yhe4*rc12a - yo16*r1216 - 2*yc12*r1212;
	J(ic12,in14)  =  yh1*yn14*afrac*rn14p;
	J(ic12,io16)  = -yc12*r1216;

	// n14 jacobian elements
	J(in14,ih1)   =  yc12*rc12p - yn14*rn14p + yo16*ro16p;
	J(in14,ihe4)  = -yn14*rn14a;
	J(in14,ic12)  =  yh1*rc12p;
	J(in14,in14)  = -yh1*rn14p - yhe4*rn14a;
	J(in14,io16)  =  yo16*ro16p;
    
	// o16 jacobian elements
	J(io16,ih1)   = yn14*pfrac*rn14p - yo16*ro16p;
	J(io16,ihe4)  = yc12*rc12a - yo16*ro16a;
	J(io16,ic12)  = yhe4*rc12a - yo16*r1216;
	J(io16,in14)  = yh1*pfrac*rn14p;
	J(io16,io16)  = yh1*ro16p - yc12*r1216 -yhe4*ro16a;

	// ne20 jacobian elements
	J(ine20,ihe4)  =  yn14*rn14a + yo16*ro16a - yne20*rne20a;
	J(ine20,ic12)  =  yc12*r1212;
	J(ine20,in14)  =  yhe4*rn14a;
	J(ine20,io16)  =  yo16*ro16a;
	J(ine20,ine20) = -yhe4*rne20a;

	// mg24 jacobian elements
	J(img24,ihe4) = yne20*rne20a;
	J(img24,ic12) = yo16*r1216;
	J(img24,io16) = yc12*r1216;
	J(img24,ine20) = yhe4*rne20a;

	// energy accounting
	for (int j=0; j<niso; j++) {
	    for (int i=0; i<niso; i++) { J(iener,j) += J(i,j)*mion[i]; }
	    J(iener,j) *= avo*clight*clight;
	}
    } 
};

struct ODE {

    void operator() ( const vector_type &y, vector_type &dydt, double /* t */)
    {
	
	// EOS
	double T=y[itemp];
	double den=y[iden];
	vec T9=get_T9_array(T);

	// rates
	double rpp=den*pp_rate(T9);
	double r33=den*he3he3_rate(T9);
	double r34=den*he3he4_rate(T9);
	double r3a=den*den*triple_alpha_rate(T9);
	double rc12p=den*c12p_rate(T9);
	double rc12a=den*c12a_rate(T9);
	double rn14p=den*n14p_rate(T9);
	double rn14a=n14a_rate(T9);
	double ro16p=den*o16p_rate(T9);
	double ro16a=den*o16a_rate(T9);
	double rne20a=den*ne20a_rate(T9);
	double r1212=den*c12c12_rate(T9);
	double r1216=den*c12o16_rate(T9);

	double pfrac=n15pg_frac(T9);
	double afrac=1-pfrac;

	double   yh1 = y[  ih1];
	double  yhe3 = y[ ihe3];
	double  yhe4 = y[ ihe4];
	double  yc12 = y[ ic12];
	double  yn14 = y[ in14];
	double  yo16 = y[ io16];
	double yne20 = y[ine20];
	
	dydt[ih1]  = -1.5*yh1*yh1*rpp;
	dydt[ih1] +=  yhe3*yhe3*r33;
	dydt[ih1] += -yhe3*yhe4*r34;
	dydt[ih1] += -2*yh1*yc12*rc12p;
	dydt[ih1] += -2*yh1*yn14*rn14p;
	dydt[ih1] += -2*yh1*yo16*ro16p;
	
	dydt[ihe3]  = 0.5*yh1*yh1*rpp;
	dydt[ihe3] += -yhe3*yhe3*r33;
	dydt[ihe3] += -yhe3*yhe4*r34;
    
	dydt[ihe4]  = 0.5*yhe3*yhe3*r33;
	dydt[ihe4] +=  yhe3*yhe4*r34;
	dydt[ihe4] += -yhe4*yc12*rc12a;
	dydt[ihe4] +=  yh1*yn14*afrac*rn14p;
	dydt[ihe4] +=  yh1*yo16*ro16p;
	dydt[ihe4] += -0.5*yhe4*yhe4*yhe4*r3a;
	dydt[ihe4] += -yhe4*yo16*ro16a;
	dydt[ihe4] += 0.5*yc12*yc12*r1212;
	dydt[ihe4] +=  yc12*yo16*r1216;
	dydt[ihe4] += -yhe4*yne20*rne20a;
   
	dydt[ic12]  = (1./6.)*yhe4*yhe4*yhe4*r3a;
	dydt[ic12] += -yhe4*yc12*rc12a;
	dydt[ic12] += -yh1*yc12*rc12p;
	dydt[ic12] +=  yh1*yn14*afrac*rn14p;
	dydt[ic12] += -yc12*yc12*r1212;
	dydt[ic12] += -yc12*yo16*r1216;

	dydt[in14]  =  yh1*yc12*rc12p;
	dydt[in14] += -yh1*yn14*rn14p;
	dydt[in14] +=  yh1*yo16*ro16p;
	dydt[in14] += -yhe4*yn14*rn14a;
    
	dydt[io16]  =  yhe4*yc12*rc12a;
	dydt[io16] +=  yh1*yn14*pfrac*rn14p;
	dydt[io16] += -yh1*yo16*ro16p;
	dydt[io16] += -yc12*yo16*r1216;
	dydt[io16] += -yhe4*yo16*ro16a;

	dydt[ine20]  = 0.5*yc12*yc12*r1212;
	dydt[ine20] +=  yhe4*yn14*rn14a;
	dydt[ine20] +=  yhe4*yo16*ro16a;
	dydt[ine20] += -yhe4*yne20*rne20a;

	dydt[img24]  =  yc12*yo16*r1216;
	dydt[img24] +=  yhe4*yne20*rne20a;

	dydt[itemp] = 0.;
	dydt[iden]  = 0.;

	// energy accounting
	double enuc = 0.;
	for (int i=0; i<niso; i++) { enuc += mion[i]*dydt[i]; }
	dydt[iener] = enuc*avo*clight*clight;
    } 
};


// main program 


int main() {
    
    vector_type y(nvar,0.);

    // initialize the system using a solar composition;
    // initial temp is 10^7 K and density is 100 g/cc
    y[ih1]= 0.708;
    y[ihe3]=2.94e-5;
    y[ihe4]=0.276;
    y[ic12]=0.003;
    y[in14]=0.0011;
    y[io16]=9.62e-3;
    y[ine20]=1.62e-3;
    y[img24]=5.16e-4;

    y[itemp]=1e7;
    y[iden]=1e2;
    y[iener]=0.0;

    // set duration and initial step size; solvers will adjust
    double tmax=3.15e17;
    double dt0=1e10;
    
    //converts mass fraction to number fraction
    double ysum=0;
    for (int i=0; i < niso; i++) {
	y[i] /= aion[i];
	ysum+= y[i];
    }
    for (int i=0; i < niso; i++) {y[i] /= ysum;}



    int stiff=0; // 1 means use Rosenbrock implicit method;
                 // 0 means use explicit 5th order Dormand-Prince Runge-Kutta

    size_t num_of_steps=0;
    if ( stiff ) {
	std::cout << " *** Implicit Rosenbrock method *** " << std::endl;
	
	num_of_steps = integrate_const( make_dense_output<rosenbrock4<double>>(1.0e-6,1.0e-6) ,
					std::make_pair( ODE(), Jacobian() ), y, 0.0, tmax, dt0);
	
    } else {
	std::cout << " *** Explicit RK Dormand-Prince *** " << std::endl;
	
	num_of_steps = integrate_const( make_dense_output<runge_kutta_dopri5<vector_type>>(1.0e-6, 1.0e-6),
					ODE(), y, 0.0, tmax, dt0);
	
    }

    //convert number fraction to mass fraction
    ysum=0;
    for (int i=0; i < niso; i++) {
	y[i] *= aion[i];
	ysum+= y[i];
    }
    for (int i=0; i < niso; i++) {y[i] /= ysum;}
    
    std::cout << "    H1: " << y[ih1] << std::endl;
    std::cout << "   He4: " << y[ihe4] << std::endl;
    std::cout << "energy: " << y[iener] << std::endl;
    std::cout << " steps: " << num_of_steps << std::endl;
    return 0;
    
}
