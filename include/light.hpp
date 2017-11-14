/////////////////////
//// light.h
////
//// Declaration of the light class
////
////////////////////

#ifndef ____LIGHT__
#define ____LIGHT__

#include "constants.hpp"
//#include "params.hpp"
#include "halo_array.hpp"
#include "communicator.hpp"
#include "space.hpp"

class light{
  
private:
  communicator* mycomm;
  space* mygrid;
  int Nx;
  int Ny;
  int Nxglobal;
  int Nyglobal;
  double dx;
  double dy;
  arma::vec x_ax;
  arma::vec y_ax;
  
  double wavelength;
  double k0;
  // Coordinate arrays
  arma::cx_mat field;
  
  // Private functions. For internal use
  arma::cx_mat apply_laplacian(halo_cx_mat& halofunc);
  void communicate_halo_points(halo_cx_mat& halofield);
  
public:
  // Constructor
  light(space& grid, communicator& comm,
	double _wavel = 800.0);
  // Destructor
  ~light(){}
  
  // Functions to get / set the norm of the wavefunction
  double get_norm();
  void normalise();
  double get_wavelength() const {return wavelength;}
  void set_wavelength(const double new_wave)
  {wavelength = new_wave; k0 = twopi / wavelength;}
  
  // Apply Helmholtz equation
  void apply_helmholtz();
  
  // Solve Helmholtz equation
  int solve_helmholtz_eigenproblem(int argc,char **argv, int num_eigen_modes);
  
  // // Time evolution
  // void evolve(const laser& mypulse);
  
  // Functions to save the wavefunction to file
  void save_field_intensity(const string& name);
  void save_observable(const arma::vec& obs, const string& name);
  void save_observable(const arma::cx_vec& obs, const string& name);
  void save_observable(const arma::vec& time, const arma::vec& obs, const string& name);
  void save_observable(const arma::vec& time, const arma::cx_vec& obs, const string& name);
  
};

#endif // ____LIGHT_ //
