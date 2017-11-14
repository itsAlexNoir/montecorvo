////////////
////  space.cpp
////
////  Contains the class functions for the kspace class
////////////

#include <iostream>
#include <math.h>
#include "space.hpp"

void fdweights(double xi, arma::vec& x, int n, int m, arma::mat& c);

space::space(int _Nx, int _Ny, double _dx, double _dy,
	     int _xrulepts, int _yrulepts, const communicator& comm):
  Nx{_Nx}, Ny{_Ny}, dx{_dx}, dy{_dy},
  xrulepts{_xrulepts}, yrulepts{_yrulepts},
  mycomm{&comm},
  x_ax(_Nx,arma::fill::zeros), y_ax(Ny,arma::fill::zeros),
  xhalo_ax(_Nx,_xrulepts), yhalo_ax(_Ny,_yrulepts),
  xfdcoeffs(1,3,_xrulepts), yfdcoeffs(1,3,_yrulepts),
  nfield(_Nx,_Ny,arma::fill::zeros)
{

  // Number of half points of the FD rule
  xfdpts = double(xrulepts - 1) * 0.5;
  yfdpts = double(yrulepts - 1) * 0.5;
  
  // Set number of global points
  Nxglobal = Nx * mycomm->get_numproc1dx();
  Nyglobal = Ny * mycomm->get_numproc1dy();
  
  // Set boundary values
  xmax = double(Nxglobal - 1) * 0.5 * dx;
  ymax = double(Nyglobal - 1) * 0.5 * dy;
  
  xminlocal = -xmax + double(Nx * mycomm->get_ipx()) * dx;
  yminlocal = -ymax + double(Ny * mycomm->get_ipy()) * dy;

  xmaxlocal = xminlocal + double(Nx - 1) * dx;
  ymaxlocal = yminlocal + double(Ny - 1) * dy;
  
  // Creates coordinates axes 
  x_ax = arma::linspace<arma::vec>(xminlocal,xmaxlocal,Nx);
  y_ax = arma::linspace<arma::vec>(yminlocal,ymaxlocal,Ny);

  // Create halo axes
  xhalomin = xminlocal - xfdpts * dx;
  xhalomax = xmaxlocal + xfdpts * dx;
  yhalomin = yminlocal - yfdpts * dy;
  yhalomax = ymaxlocal + yfdpts * dy;
  
  xhalo_ax = arma::linspace<arma::vec>(xhalomin,xhalomax,Nx+xrulepts-1);
  yhalo_ax = arma::linspace<arma::vec>(yhalomin,yhalomax,Ny+yrulepts-1);
  
  // Set FD coefficients  
  arma::vec xpoints = arma::linspace<arma::vec>(-xfdpts,xfdpts,xrulepts) * dx;
  arma::vec ypoints = arma::linspace<arma::vec>(-yfdpts,yfdpts,yrulepts) * dy;
  
  arma::mat xcoeff(xrulepts,3,arma::fill::zeros);
  arma::mat ycoeff(yrulepts,3,arma::fill::zeros);
  
  fdweights(0.0, xpoints, xrulepts, 2, xcoeff);
  fdweights(0.0, ypoints, yrulepts, 2, ycoeff);
  
  xfdcoeffs = xcoeff;
  yfdcoeffs = ycoeff;

  // We set the refractive index of the medium
  set_refractive_index();
  
}

///////////////////////////////////////////////////////////////////

void space::print_grid_parameters()
{
  cout << endl;
  cout << "----- Grid parameters -------------------------------------------------" << endl;
  cout << " Number of grid points in x per proc:  " << Nx << endl;
  cout << " Number of grid points in y per proc:  " << Ny << endl;
  cout << " Number of total grid points in x:     " << Nx * mycomm->get_numproc1dx() << endl;
  cout << " Number of total grid points in y:     " << Ny * mycomm->get_numproc1dy() << endl;
  cout << " Grid spacing in x (au):               " << dx << endl;
  cout << " Grid spacing in y (au):               " << dy << endl;
  cout << " Finite Difference rule in x:          " << xrulepts << endl;
  cout << " Finite Difference rule in y:          " << yrulepts << endl;
  cout << " Extent in x (au):                     " << get_xmax() << endl;
  cout << " Extent in y (au):                     " << get_ymax() << endl;  
  cout << "------------------------------------------------------------------------" << endl;
  cout << endl << endl; 
 
}

///////////////////////////////////////////////////////////////////

void space::set_refractive_index()
{
  
  double rad {0.0};
  double r0 {4.0};
  double n1 {1.445};
  double n2 {1.4378};
  
  double ome {0.1};
  
  for(int ix=0; ix<Nx; ix++)
    for(int iy=0; iy<Ny; iy++)
      {
  	rad = sqrt( x_ax(ix) * x_ax(ix) +
  		    y_ax(iy) * y_ax(iy));
  	if(rad <= r0)
  	  nfield(ix,iy) = n1;
  	else
  	  nfield(ix,iy) = n2;    
      }
  
  // for(int ix=0; ix<Nx; ix++)
  //   for(int iy=0; iy<Ny; iy++)
  //     {
  // 	rad = sqrt( x_ax(ix) * x_ax(ix) +
  // 		    y_ax(iy) * y_ax(iy));
  // 	if(abs(x_ax(ix)) <= r0  && abs(y_ax(iy))<=r0)
  // 	  nfield(ix,iy) = n1;
  // 	else
  // 	  nfield(ix,iy) = n2;    
  //     }
  
  // for(int ix=0; ix<Nx; ix++)
  //   for(int iy=0; iy<Ny; iy++)
  //     {
  // 	nfield(ix, iy) = 0.5 * ome * ome *
  // 	  ( x_ax(ix) * x_ax(ix) + y_ax(iy) * y_ax(iy) );
  //     }	
}

///////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------
//
//  SUBROUTINE fdweights
//
///> \brief Create finite-difference weights.
///> \details Create coefficients to be used for
///> finite-difference rules.
///> Reference:
///> Generation of Finite Difference Formulas on Arbitrarily
///> Spaced Grids, Bengt Fornberg,
///> Mathematics of compuation, 51, 184, 1988, 699-706
///
///> \param[in] xi Location at which approximations are to be accurate (Central grid point).
///> \param[in] x Grid point locations x(0:n).
///> \param[in] n Number of grid points. Also, order of the finite-difference rule.
///> \param[in] m Order of the highest derivative (it includes 0th-order derivative).
///> \param[out] Finite-difference coefficients.
///
///----------------------------------------------------------------------------

void fdweights(double xi, arma::vec& x, int n, int m, arma::mat& c)
{  
  
  double c1, c4, c5;
  c1 = 1.0;
  c4 = x(0) - xi;
  
  for(int i=0; i<n; i++)
    for(int k=0; k<m; k++)
      c(i,k) = 0.0;
  
  c(0,0) = 1.0;
  
  for (int i=1; i < n; ++i){
    int mn = std::min(i, m);
    double c2 = 1.0;
    c5 = c4;
    c4 = x(i) - xi;
    for (int j=0; j<i; ++j){
      double c3 = x(i) - x(j);
      c2 = c2*c3;
      if (j == i-1){
	for (int k=mn; k>=1; --k)
	  c(i, k) = c1 * (k*c(i-1,k-1) - c5*c(i-1,k)) /c2;
	
	c(i, 0) = -c1*c5*c(i-1, 0)/c2;
      }
      for (int k=mn; k>=1; --k){
	c(j, k) = (c4*c(j, k) - k*c(j, k-1)) /c3;
      }
      c(j, 0) = c4*c(j, 0)/c3;
    }
    c1 = c2;
  }
}  

///////////////////////////////////////////////////////////////////
