/*!
  @file light.hpp
  
  @brief Declaration of the light class
  
*/

#ifndef ____LIGHT__
#define ____LIGHT__

#include <slepceps.h>
#include "constants.hpp"
#include "halo_array.hpp"
#include "communicator.hpp"
#include "fiber.hpp"
#include "io.hpp"

class light{
  
private:
  communicator* mycomm;
  fiber* mygrid;
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

  // Helmholtz matrix (for eigencalculation)
  Mat            Hmat;
  Vec            xr,xi;
  // Eigenproblem solver context
  EPS            eps;
  // Coordinate arrays
  arma::cx_mat field;
  
  // Private functions. For internal use
  arma::cx_mat apply_laplacian(halo_cx_mat& halofunc);
  void communicate_halo_points(halo_cx_mat& halofield);
  
public:
  // Constructor
  light(fiber& grid, communicator& comm,
	double _wavel,int argc, char **argv);
  // Destructor
  ~light(){
    PetscErrorCode ierr;
    // Free PETSC's work space
    ierr = EPSDestroy(&eps);
    ierr = MatDestroy(&Hmat);
    ierr = VecDestroy(&xr);
    ierr = VecDestroy(&xi);
    ierr = SlepcFinalize();
  }
  
  // Functions to get / set the norm of the wavefunction
  double get_norm();
  void normalise();
  double get_wavelength() const {cout << k0 << endl; return wavelength;}
  void set_wavelength(const double new_wave)
  {wavelength = new_wave; k0 = twopi / wavelength;}
  
  // Apply Helmholtz equation
  void apply_helmholtz();
  
  // Solve Helmholtz equation
  int solve_helmholtz_eigenproblem(int argc,char **argv, int num_eigen_modes,
				   bool save_to_hdf5 = true,
				   arma::uvec selected_modes = {0});
  
  // // Time evolution
  // void evolve(const laser& mypulse);
  
  // Functions to save the wavefunction to file
  void save_field_intensity(const string& name);  
};

#endif // ____LIGHT_ //
