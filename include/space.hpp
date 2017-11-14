/////////////////////
//// space.h
////
//// Declaration of the kspace class
////
////////////////////

#ifndef ____SPACE__
#define ____SPACE__

#include "constants.hpp"
//#include "params.hpp"
#include "halo_array.hpp"
#include "communicator.hpp"

class space{
  
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

  halo1D_mat xfdcoeffs;
  halo1D_mat yfdcoeffs;

  double xmax, ymax;
  double xminlocal, xmaxlocal;
  double yminlocal, ymaxlocal;
  double xhalomin, xhalomax;
  double yhalomin, yhalomax;
  
 public:
  // Constructor
  space(int nx, int ny, double dx, double dy,
	int xrulepts, int yrulepts,
	const communicator& comm);
  // Destructor
  ~space(){}
  
  // Functions that returns private values
  int get_Nx() const;
  int get_Ny() const;
  double get_dx() const;
  double get_dy() const;
  double get_xmax() const;
  double get_ymax() const;
  
  const arma::vec get_x_ax() const;
  const arma::vec get_y_ax() const;
  
  int get_xrulepts() const;
  int get_yrulepts() const;
  
  double get_xfdcoeffs(const int i, const int j);
  double get_yfdcoeffs(const int i, const int j);
  double get_nfield(const int i, const int j);
  const arma::mat& get_nfield() const;
  
  // Function that print class members
  void print_grid_parameters();

  // Set refractive index for the space
  void set_refractive_index();
  
};


inline int space::get_Nx() const {return Nx;}
inline int space::get_Ny() const {return Ny;}

inline double space::get_dx() const {return dx;}
inline double space::get_dy() const {return dy;}

inline const arma::vec space::get_x_ax() const {return x_ax;}
inline const arma::vec space::get_y_ax() const {return y_ax;}

inline double space::get_xmax() const {return xmax;}
inline double space::get_ymax() const {return ymax;}

inline int space::get_xrulepts() const {return xrulepts;}
inline int space::get_yrulepts() const {return yrulepts;}

inline double space::get_xfdcoeffs(const int i, const int j)
{ return xfdcoeffs(i, j);}

inline double space::get_yfdcoeffs(const int i, const int j)
{ return yfdcoeffs(i, j);}

inline double space::get_nfield(const int i, const int j)
{ return nfield(i, j);}

inline const arma::mat& space::get_nfield() const { return nfield;}


#endif // ____SPACE__ //
