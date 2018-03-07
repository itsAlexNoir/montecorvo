/*!
  @file space.hpp

  @brief Declaration of the fiber class.
*/

#ifndef ____FIBER__
#define ____FIBER__

#include "constants.hpp"
#include "halo_array.hpp"
#include "communicator.hpp"

struct shots_strc{
  int Nx;
  int Ny;
  double dx;
  double dy;
  int Nc;
  double** coords;
};
  
class fiber{
  
private:
  // Communicator for parallel stuff
  const communicator *mycomm;
  
  // Number of grid points
  int Nx;
  int Ny;
  int Nxglobal;
  int Nyglobal;
  
  // Finite-difference
  int xrulepts, yrulepts;
  int xfdpts, yfdpts;
  
  // Grid spacings
  double dx;
  double dy;
  
  // Coordinate arrays
  arma::vec x_ax;
  arma::vec y_ax;
  halo_vec xhalo_ax;
  halo_vec yhalo_ax;
  arma::mat nfield;
  arma::mat imag_nfield;
  shots_strc shots;
  
  halo1D_mat xfdcoeffs;
  halo1D_mat yfdcoeffs;
  bool absorption_on;
  
  double xmax, ymax;
  double xminlocal, xmaxlocal;
  double yminlocal, ymaxlocal;
  double xhalomin, xhalomax;
  double yhalomin, yhalomax;
  
  void read_shot_coordinates(string filename, int nx, int ny,
			     double dx, double dy);
  
public:
  // Constructor
  fiber(int nx, int ny, double dx, double dy,
	int xrulepts, int yrulepts,
	bool abs_on, const communicator& comm);
  // Destructor
  ~fiber(){
    for(int ic=0; ic<shots.Nc; ic++)
      delete[] shots.coords[ic];
    
    delete[] shots.coords;
  }
  
  // Functions that returns private values
  int get_Nx() const;
  int get_Ny() const;
  int get_Nxglobal() const;
  int get_Nyglobal() const;
  double get_dx() const;
  double get_dy() const;
  double get_xmax() const;
  double get_ymax() const;
  void absorption_switch(bool state) {absorption_on = state;}
  bool get_absorption_state() {return absorption_on;}
  
  const arma::vec get_x_ax() const;
  const arma::vec get_y_ax() const;
  
  int get_xrulepts() const;
  int get_yrulepts() const;
  
  double get_xfdcoeffs(const int i, const int j);
  double get_yfdcoeffs(const int i, const int j);
  double get_nfield(const int i, const int j);
  double get_imag_nfield(const int i, const int j);
  const arma::mat& get_nfield() const;
  const arma::mat& get_imag_nfield() const;
  
  // Function that print class members
  void print_grid_parameters();

  // Set fiber structure
  double get_refractive_index(string material, double wavelength);
  void set_step_index_fiber(double r0, double n1, double n2);
  void set_circular_honeycomb_fiber(double r0, int no_holes,
				    double n0, double dn,
				    double ddx, double ddy,
				    int exponent=8);
  
  void set_honeycomb_fiber(string filename,
			   double n0, double dn,
			   int Nshx, int Nshy,
			   double deltashx, double deltashy,
			   double ddx, double ddy,
			   int exponent);  
  
  void set_fiber_cladding(double n1, double rclad);
  void set_fiber_absorption(double n1, double rclad, int exponent);
};


inline int fiber::get_Nx() const {return Nx;}
inline int fiber::get_Ny() const {return Ny;}

inline int fiber::get_Nxglobal() const {return Nxglobal;}
inline int fiber::get_Nyglobal() const {return Nyglobal;}

inline double fiber::get_dx() const {return dx;}
inline double fiber::get_dy() const {return dy;}

inline const arma::vec fiber::get_x_ax() const {return x_ax;}
inline const arma::vec fiber::get_y_ax() const {return y_ax;}

inline double fiber::get_xmax() const {return xmax;}
inline double fiber::get_ymax() const {return ymax;}

inline int fiber::get_xrulepts() const {return xrulepts;}
inline int fiber::get_yrulepts() const {return yrulepts;}

inline double fiber::get_xfdcoeffs(const int i, const int j)
{ return xfdcoeffs(i, j);}

inline double fiber::get_yfdcoeffs(const int i, const int j)
{ return yfdcoeffs(i, j);}

inline double fiber::get_nfield(const int i, const int j)
{ return nfield(i, j);}
inline double fiber::get_imag_nfield(const int i, const int j)
{ return imag_nfield(i, j);}

inline const arma::mat& fiber::get_nfield() const { return nfield;}
inline const arma::mat& fiber::get_imag_nfield() const { return imag_nfield;}


#endif // ____FIBER__ //
