/////////////////////
//// halo_array.h
////
//// Declaration of the halo_array class
////
////////////////////

#ifndef ____HALO_ARRAY__
#define ____HALO_ARRAY__
#include "constants.hpp"
#include <armadillo>

class halo_vec{

  arma::vec array;
  int Nx;
  int rulepts;
  int fdpts;
  
public:
  /// Constructor
  halo_vec(const int _Nx, const int _rulepts):
    Nx{_Nx},rulepts{_rulepts},
    array(_Nx+_rulepts-1,arma::fill::zeros){fdpts = (rulepts-1)*0.5;}
  /// Destructor
  ~halo_vec(){}
  
  // Overloaded operators
  inline double operator()(const int i) const {return array(i+fdpts);}
  inline halo_vec& operator=(const arma::vec& eqarr)
  {array = eqarr; return *this;}
  inline void zeros() {array.zeros();}
};

class halo1D_mat{

  arma::mat array;
  int Nx, Ny;
  int rulepts;
  int fdpts;
  
public:
  /// Constructor
  halo1D_mat(const int _Nx, const int _Ny, const int _rulepts):
    Nx{_Nx},Ny{_Ny},rulepts{_rulepts},
    array(_Nx+_rulepts-1,_Ny,arma::fill::zeros){fdpts = (rulepts-1)*0.5;}
  /// Destructor
  ~halo1D_mat(){}
  
  // Overloaded operators
  inline double operator()(const int i, const int j) const
  {return array(i+fdpts,j);}
  inline halo1D_mat& operator=(const arma::mat& eqarr)
  {array = eqarr; return *this;}
  inline void zeros() {array.zeros();}
};


class halo_mat{

  arma::mat array;
  int Nx, Ny;
  int xrulepts, yrulepts;
  int xfdpts, yfdpts;
  
public:
  /// Constructor
  halo_mat(const int _Nx, const int _Ny,
	   const int _xrulepts, const int _yrulepts):
    Nx{_Nx}, Ny{_Ny},
    xrulepts{_xrulepts}, yrulepts{_yrulepts},
    array(_Nx+_xrulepts-1,_Ny+_yrulepts-1,arma::fill::zeros)
  {xfdpts = (xrulepts-1)*0.5; yfdpts = (yrulepts-1)*0.5;}
  /// Destructor
  ~halo_mat(){}
  
  // Overloaded operators
  inline const double operator()(const int i, const int j) const
  {return array(i+xfdpts,j+yfdpts);}
  inline halo_mat& operator=(const arma::mat& eqmat)
  {array = eqmat; return *this;}
  inline void zeros() {array.zeros();}
};

class halo_cx_mat{

  arma::cx_mat array;
  int Nx, Ny;
  int xrulepts, yrulepts;
  int xfdpts, yfdpts;
  
public:
  /// Constructor
  halo_cx_mat(const int _Nx, const int _Ny,
	      const int _xrulepts, const int _yrulepts):
    Nx{_Nx}, Ny{_Ny},
    xrulepts{_xrulepts}, yrulepts{_yrulepts},
    array(_Nx+_xrulepts-1,_Ny+_yrulepts-1,arma::fill::zeros)
  {xfdpts = (xrulepts-1)*0.5; yfdpts = (yrulepts-1)*0.5;}
  /// Destructor
  ~halo_cx_mat(){}
  
  // Overloaded operators
  inline dcomplex operator()(const int i, const int j) const
  {return array(i+xfdpts,j+yfdpts);}
  inline halo_cx_mat& operator=(const arma::cx_mat& eqcmat)
  {array = eqcmat; return *this;}

  inline void zeros() {array.zeros();}
  inline void set_value(const int i, const int j, const dcomplex value)
  {array(i+xfdpts, j+yfdpts) = value;}
};

#endif // ____HALO_ARRAY__ //
