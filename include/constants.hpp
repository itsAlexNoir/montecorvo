/////////////
///
///    Constants.h
///         Constant miscellanea!!
///
////////


#ifndef ____CONSTANTS__
#define ____CONSTANTS__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <array>
#include <complex.h>
#include <chrono>
#include <ctime>
#define ARMA_DONT_USE_WRAPPER
// Take care with this macro: Use wisely!!!
#define ARMA_NO_DEBUG
#include <armadillo>

typedef std::complex<double> dcomplex;
/* typedef std::vector<std::vector<double>> matrix; */
/* typedef std::vector<std::vector<dcomplex>> cmatrix; */

using namespace std;
//using namespace arma;
//cout.precision(15);

/// Complex numbers
const dcomplex Im {0.0,1.0};
const dcomplex One {1.0,0.0};
const dcomplex Zero {0.0,0.0};

/// pi
const double pi = 4.0*atan(1.0);
const double twopi = 2.0 * pi;
const double one_over_pi = 1.0 / pi;

/// Square root two
const double sqrt_two = sqrt(2.0);
const double one_over_sqrt_two = 1.0 / sqrt_two;

/// Speed of light, in a.u. units
const double speed_light = 137.0; 
const double one_over_speed_light = 1.0 / speed_light;
// Speed of light in meters over seconds 
const double speed_light_is = 2.99792458e8;

/// electron charge
const double charge_electron_au = -1.;
/// Proton mass
const double mass_proton_au = 1836.;
/// Bohr radius (IS)
const double bohr_radius  = 5.29177211e-11;
const double aulength_nm  = 5.29177211e-2;
const double hatree_to_eV = 27.211;
// au time
const double autime_fs = 2.418884326505e-2;
const double autime_s  = 2.418884326505e-17;
// au intensity
const double auintensity = 3.5e16;

#endif  // ____constants__ //


