/*! \file fiber.cpp
  \brief Definitions for the class fiber.

    File containing the definitions of the class fiber.
*/

#include <iostream>
#include <math.h>
#include "fiber.hpp"

void fdweights(double xi, arma::vec& x, int n, int m, arma::mat& c);

/*! \class fiber
    \brief Fiber class.

    A more detailed class description.
*/
fiber::fiber(int _Nx, int _Ny, double _dx, double _dy,
	     int _xrulepts, int _yrulepts, bool abs_on,
	     const communicator& comm):
  Nx{_Nx}, Ny{_Ny}, dx{_dx}, dy{_dy},
  xrulepts{_xrulepts}, yrulepts{_yrulepts},
  absorption_on{abs_on},
  mycomm{&comm},
  x_ax(_Nx,arma::fill::zeros), y_ax(Ny,arma::fill::zeros),
  xhalo_ax(_Nx,_xrulepts), yhalo_ax(_Ny,_yrulepts),
  xfdcoeffs(1,3,_xrulepts), yfdcoeffs(1,3,_yrulepts),
  nfield(_Nx,_Ny,arma::fill::zeros),
  imag_nfield(_Nx,_Ny,arma::fill::zeros)
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
  
  // Allocate this variable for the first time.
  // This operation is needed for memory allocation.
  shots.Nc = 2;
  shots.coords = new double*[shots.Nc];
  for(int ic=0; ic<shots.Nc; ic++)
    shots.coords[ic] = new double[2];
  
  // // Save gridpoints to file
  // ofstream xptsfile;
  // ofstream yptsfile;
  
  // string filename;
  // int iproc = mycomm->get_iprocessor();
  // filename = "./gridpoints/xpts." + to_string(iproc) + ".dat";
  // xptsfile.open(filename);
  // filename = "./gridpoints/ypts." + to_string(iproc) + ".dat";
  // yptsfile.open(filename);

  // for(int ix=-xfdpts; ix<Nx+xfdpts; ix++)
  //   xptsfile << ix << " " << xhalo_ax(ix) << endl;

  // for(int iy=-yfdpts; iy<Ny+yfdpts; iy++)
  //   yptsfile << iy << " " << yhalo_ax(iy) << endl;

  // xptsfile.close();
  // yptsfile.close();

}

///////////////////////////////////////////////////////////////////

void fiber::print_grid_parameters()
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

double fiber::get_refractive_index(string material, double wavelength)
{
  double wavesq = wavelength * wavelength;
  double nsq {0.0};
  
  if(material=="YAG")
    nsq = 1.0 + 2.28200 * wavesq / (wavesq - 0.01185) +
      3.27644 * wavesq / (wavesq - 282.734);
  else
    cout << "No Sellmeier formula for this material" << endl;
  
  return sqrt(nsq);
}

///////////////////////////////////////////////////////////////////

void fiber::set_step_index_fiber(double r0, double n1, double n2)
{
  
  double rad {0.0};
  
  for(int iy=0; iy<Ny; iy++)
    for(int ix=0; ix<Nx; ix++)
      {
  	rad = sqrt( x_ax(ix) * x_ax(ix) +
  		    y_ax(iy) * y_ax(iy));
  	if(rad <= r0)
  	  nfield(ix,iy) = n1;
  	else
  	  nfield(ix,iy) = n2;    
      }  
}

///////////////////////////////////////////////////////////////////

void fiber::set_circular_honeycomb_fiber(double r0, int no_sides,
					 double xfactor, double yfactor,
					 double n0, double dn,
					 double ddx, double ddy,
					 int exponent)
{
  
  double argx, argy;  
  int no_holes = (int) round(2.0 * pi * r0 / double(no_holes));
  double dtheta {twopi / no_holes};
  arma::vec angholes = arma::linspace(0.0, twopi-dtheta,no_holes);
  
  arma::vec x0 = xfactor * r0 * arma::cos(angholes);
  arma::vec y0 = yfactor * r0 * arma::sin(angholes);

  
  for(int ih=0; ih<no_holes; ih++)
    for(int iy=0; iy<Ny; iy++)
      for(int ix=0; ix<Nx; ix++)
	{
	  argx = ( x_ax(ix) - x0(ih) ) / ddx;
	  argy = ( y_ax(iy) - y0(ih) ) / ddy; 
	  argx = pow(argx,exponent);
	  argy = pow(argy,exponent);
	  nfield(ix, iy) += n0 + dn * exp(-argx) * exp(-argy);
	}
  
}

///////////////////////////////////////////////////////////////////

void fiber::read_shot_coordinates(string filename, int nx, int ny,
				  double dx, double dy)
{
  
  shots.Nx  = nx;
  shots.Ny  = ny;
  shots.dx  = dx;
  shots.dy  = dy;

  ifstream shotsfile(filename,ios_base::in);
  arma::Mat<int> binshots(shots.Nx, shots.Ny,arma::fill::zeros); 
  
  // Read positions from file
  if(shotsfile.is_open())
    {
      for(int iy=0; iy<shots.Ny; iy++)
	for(int ix=0; ix<shots.Nx; ix++)
	  shotsfile >> binshots(ix,iy);
      
      shotsfile.close();
    }
  else
    cout << "Unable to open file";
    
  // First delete previous allocation for this variable
  for(int ic=0; ic<shots.Nc; ic++)
    delete[] shots.coords[ic];
  delete[] shots.coords;
  
  shots.Nc = arma::accu(binshots);
  shots.coords = new double*[shots.Nc];
  for(int ic=0; ic<shots.Nc; ic++)
    shots.coords[ic] = new double[2];
  
  // Write shot positions
  ofstream shotcoordsfile("shotsfile.dat",ofstream::out);
  double incx = 0.0;
  int ic = -1;
  for(int iy=0; iy<shots.Ny; iy++)
    {
      for(int ix=0; ix<shots.Nx; ix++)
	if(binshots(ix, iy) == 1)
	  {
	    ic += 1;
	    shots.coords[ic][0] = (- double(shots.Nx) * 0.5 + double(ix) )
	      * shots.dx + incx;
	    shots.coords[ic][1] = (- double(shots.Ny - 1) * 0.5 + double(iy) )
	      * shots.dy;
	    
	    shotcoordsfile << shots.coords[ic][0] << " " << shots.coords[ic][1] << endl;
	  }
      incx += shots.dx * 0.5;
      if(incx==shots.dx) incx=0.0;    
    }
  
  shotcoordsfile.close();
}

///////////////////////////////////////////////////////////////////

void fiber::set_honeycomb_fiber(string filename,
				double n0, double dn,
				int Nshx, int Nshy,
				double deltashx, double deltashy,
				double ddx, double ddy,
				int exponent)
{

  double argx, argy;  
  // Get shot's coordinates
  read_shot_coordinates(filename, Nshx, Nshy,
			deltashx, deltashy);
  
  for(int iy=0; iy<Ny; iy++)
    for(int ix=0; ix<Nx; ix++)
      {	
	nfield(ix, iy) = n0;
	
	for(int ih=0; ih<shots.Nc; ih++)
	  {
	    
	    argx = ( x_ax(ix) - shots.coords[ih][0] ) / ddx;
	    argy = ( y_ax(iy) - shots.coords[ih][1] ) / ddy; 
	    argx = pow(argx,exponent);
	    argy = pow(argy,exponent);
	    nfield(ix, iy) += dn * exp(-argx) * exp(-argy);
	    
	  }
      }
    
}

///////////////////////////////////////////////////////////////////

void fiber::set_fiber_cladding(double n1, double rclad)
{
  
  double rad {0.0};
  for(int iy=0; iy<Ny; iy++)
    for(int ix=0; ix<Nx; ix++)
      {	  
	rad = sqrt(x_ax(ix) * x_ax(ix)
		   + y_ax(iy) * y_ax(iy) );
	
	if(rad > rclad)
	  nfield(ix, iy) = n1;
      }
 
}

void fiber::set_fiber_absorption(double n1, double rclad, int exponent)
{
  
  double rad {0.0};
  double arg_abs;
  double width_abs = xmax - rclad;
  for(int iy=0; iy<Ny; iy++)
    for(int ix=0; ix<Nx; ix++)
      {
	rad = sqrt(x_ax(ix) * x_ax(ix)
		   + y_ax(iy) * y_ax(iy) );
	arg_abs = (rad - rclad) / width_abs;
	arg_abs = pow(arg_abs,exponent);
	
	if(rad > rclad)
	  imag_nfield(ix, iy) = n1 * arg_abs;
	
      }
}

///////////////////////////////////////////////////////////////////

/*!----------------------------------------------------------------------------
  
  @function fdweights
  
  @brief Create finite-difference weights.
  @details Create coefficients to be used for
  finite-difference rules.
  Reference:
  Generation of Finite Difference Formulas on Arbitrarily
  Spaced Grids, Bengt Fornberg,
  Mathematics of compuation, 51, 184, 1988, 699-706
  
  @param[in] xi Location at which approximations are to be accurate (Central grid point).
  @param[in] x Grid point locations x(0:n).
  @param[in] n Number of grid points. Also, order of the finite-difference rule.
  @param[in] m Order of the highest derivative (it includes 0th-order derivative).
  @param[out] Finite-difference coefficients.
  
  ----------------------------------------------------------------------------
*/
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
