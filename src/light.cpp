////////////
////  light.cpp
////
////  Contains the class functions for the light class
////////////
#include <iostream>
#include <math.h>
#include <slepceps.h>
#include "light.hpp"

light::light(const space& grid, const communicator& comm,
	     double _wavel):
  mygrid{&grid}, mycomm{&comm}, wavelength{_wavel},
  field(grid.get_Nx(),grid.get_Ny(),arma::fill::zeros)
{
  // Asign value to private members
  Nx            = mygrid->get_Nx();
  Ny            = mygrid->get_Ny();
  dx            = mygrid->get_dx();
  dy            = mygrid->get_dy();  
  x_ax          = mygrid->get_x_ax();
  y_ax          = mygrid->get_y_ax();
  
  k0 = twopi / wavelength;
  
}

//////////////////////////////////////////////////

double light::get_norm()
{
  double norm = 0.0;
  
  for(int ix=0; ix<Nx; ++ix)
    for(int iy=0; iy<Ny; ++iy)
      norm += real(conj(field(ix, iy))
		   * field(ix, iy));
  
  return norm *= dx * dy;
}

void light::normalise()
{
  double norm = get_norm();
  norm = 1.0 / sqrt(norm);
  field *= norm;
  
}

//////////////////////////////////////////////////

//////////////////////////////////////////////////

// Save wavefunction to file
void light::save_field_intensity(const string& name)
{ 
  ofstream myfile;
  
  string filename = name + ".dat";
  myfile.open(filename);
  
  double prob = 0.0;
  
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      {
	prob = real(conj(field(ix,iy)) * field(ix,iy));
	myfile << prob << endl;
      }
  myfile.close();
}

// For complex observables
void light::save_observable(const arma::vec& obs, const string& name)
{
  ofstream file;
  
  string filename = name + ".dat";
  file.open(filename);
  
  for(auto i=0; i<obs.size(); ++i)
    file << obs(i) << endl;
  
  file.close();
  
}
//////////////////////////////////
// For complex observables
void light::save_observable(const arma::cx_vec& obs, const string& name)
{
  ofstream file;
  
  string filename = name + ".dat";
  file.open(filename);
  
  for(auto i=0; i<obs.size(); ++i)
    file << real(obs(i)) << " " << imag(obs(i)) << endl;
  
  file.close();
  
}

// For complex observables
void light::save_observable(const arma::vec& time, const arma::vec& obs,
			   const string& name)
{
  ofstream file;
  
  string filename = name + ".dat";
  file.open(filename);
  
  for(auto i=0; i<obs.size(); ++i)
    file << time(i) << " " << obs(i) << endl;
  
  file.close();
  
}
//////////////////////////////////

void light::save_observable(const arma::vec& time, const arma::cx_vec& obs,
			   const string& name)
{
  ofstream file;
  
  string filename = name + ".dat";
  file.open(filename);
  
  for(auto i=0; i<obs.size(); ++i)
    file << time(i) << " " << real(obs(i)) << " "
	 << imag(obs(i)) << endl;
  
  file.close();
  
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

arma::cx_mat light::apply_laplacian(const halo_cx_mat& halofunc)
{
  int xfdpts = (mygrid->get_xrulepts() - 1) * 0.5;
  int yfdpts = (mygrid->get_yrulepts() - 1) * 0.5;
  arma::cx_mat lapfield(Nx,Ny,arma::fill::zeros);
  
  // Apply 2nd derivative along the x axis
  for(int ix=0; ix<Nx; ix++)
    for(int ifd=-xfdpts; ifd<xfdpts+1; ifd++)
      for(int iy=0; iy<Ny; iy++)
  	lapfield(ix,iy) += mygrid->get_xfdcoeffs(ifd,2) * halofunc(ix+ifd,iy);
  
  // Apply 2nd derivative along the y axis
  for(int ix=0; ix<Nx; ix++)
    for(int iy=0; iy<Ny; iy++)
      for(int ifd=-yfdpts; ifd<yfdpts+1; ifd++)
  	lapfield(ix,iy) += mygrid->get_yfdcoeffs(ifd,2) * halofunc(ix, iy+ifd);
  
  return lapfield;
}

///////////////////////////////////////////////////////////////////

void light::apply_helmholtz()
{
  // Create field matrix with halo points
  halo_cx_mat halofield(Nx,Ny,
			mygrid->get_xrulepts(),mygrid->get_yrulepts());
  
  // Copy data from field to halo halofield matrices
  arma::cx_mat lapfield(Nx,Ny,arma::fill::zeros);
  
  // for(int ix=0; ix<Nx; ix++)
  //   for(int iy=0; iy<Ny; iy++)
  //     halofield.set_value(ix, iy, field(ix, iy) );

  communicate_halo_points(halofield);
  
  // Apply Laplacian operator
  lapfield = apply_laplacian(halofield);
  
  // Apply potential
  for(int ix=0; ix<Nx; ix++)
    for(int iy=0; iy<Ny; iy++)
      lapfield(ix, iy) += mygrid->get_nfield(ix, iy) * mygrid->get_nfield(ix, iy) *
  	k0 * k0 * field(ix, iy);
  
  field = lapfield;
}

///////////////////////////////////////////////////////////////////

int light::solve_helmholtz_eigenproblem(int argc,char **argv)
{
  
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  PetscScalar    matelem;
  Vec            xr,xi;
  PetscInt       n = Nx * Ny;
  PetscInt       k, l;
  PetscInt       i,Istart,Iend,nev,maxit,its,nconv;
  PetscErrorCode ierr;
  int            xrulepts = mygrid->get_xrulepts();
  int            yrulepts = mygrid->get_yrulepts();
  
  static char help[] = "Solving Helmholtz equation.";
  
  SlepcInitialize(&argc,&argv,(char*)0,help);
  
#if !defined(PETSC_USE_COMPLEX)
  SETERRQ(PETSC_COMM_WORLD,1,"This example requires complex numbers");
#endif
  
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%D\n\n",n);CHKERRQ(ierr);
  
  // Create the matrix that defines the eigensystem
  
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  
  k = -1;
  for(int ix=0; ix<Nx; ix++)
    for(int iy=0; iy<Ny; iy++)
      {
	k += 1;
	field.zeros();
	field(ix, iy) = One;
	
	apply_helmholtz();
	
	l = -1;
	for(int ixx=0; ixx<Nx; ixx++)
	  for(int iyy=0; iyy<Ny; iyy++)
	    {
	      l += 1;
	      if( (abs(k-l)<xrulepts-1) || (abs(k-l)<yrulepts-1))
		{
		  matelem = field(ixx, iyy);
		  ierr = MatSetValue(A, l, k, matelem, INSERT_VALUES);
		  //cout << "l " << l << " k " << k << " mat: " << matelem << endl;
		  CHKERRQ(ierr);
		}
	    }
      } 
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);
  
  // Create the eigensolver and set various options
  /*
    Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  
  // Set operators. In this case, it is a standard eigenvalue problem
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  //ierr = EPSSetType(eps,EPSTRLAN);
  //EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);
  ierr = EPSSetTolerances(eps,1e-8,PETSC_DEFAULT);
  
  // Set solver parameters at runtime
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  
  ///////////////////////////////
  //  Solve the eigensystem
  ///////////////////////////////
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  
  // Optional: Get some information from the solver and display it
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  // Display solution and clean up
   /*
     Get number of converged approximate eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");CHKERRQ(ierr);

    for (i=0;i<nconv;i++) {
      //  Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
      // ki (imaginary part)
      
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
      
      // Compute the relative error associated to each eigenpairs
      ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);
      
#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if (im!=0.0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)re,(double)error);CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }
  
  // for(int ix=0; ix<Nx; ixx)    
  //   cout << PetscRealPart(xr[ix]) << PetscRealPart(xi[ix]) << endl;
  
  // Free work space
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&xr);CHKERRQ(ierr);
  ierr = VecDestroy(&xi);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
  
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

void light::communicate_halo_points(halo_cx_mat& halofield)
{
  
  int numproc1dx = mycomm->get_numproc1dx();
  int numproc1dy = mycomm->get_numproc1dy();
  int maxproc1dx = mycomm->get_maxproc1dx();
  int maxproc1dy = mycomm->get_maxproc1dy();
  int ipx        = mycomm->get_ipx();
  int ipy        = mycomm->get_ipy();
  
  int xrulepts = mygrid->get_xrulepts();
  int yrulepts = mygrid->get_yrulepts();
  int xfdpts   = (xrulepts - 1) * 0.5; 
  int yfdpts   = (yrulepts - 1) * 0.5; 
  int k        = 0;
  
  if(Nx < xrulepts)
    {
      string message = "Insufficient number of points in x for "
	+ to_string(xrulepts) + " finite difference rule." ;
      mycomm->parallel_stop(message);
    }
  
  if(Ny < yrulepts)
    {
      string message = "Insufficient number of points in y for "
	+ to_string(yrulepts) + " finite difference rule.";
      mycomm->parallel_stop(message);
    }
  
  // Set halo matrix to zero
  halofield.zeros();
  
  // Copy data to the halo matrix
  for(int ix=0; ix<Nx; ix++)
    for(int iy=0; iy<Ny; iy++)
      halofield.set_value(ix, iy, field(ix, iy) );
  
  // Do communication in x
  
  if(numproc1dx > 1)
    {
      int length = xfdpts * Ny;
      
      dcomplex *funcputl = new dcomplex[length];
      dcomplex *funcputr = new dcomplex[length];
      dcomplex *funcgetl = new dcomplex[length];
      dcomplex *funcgetr = new dcomplex[length];
      
      // Set up the send arrays
      k = -1;
      for(int ix=0; ix<xfdpts; ix++)
	for(int iy=0; iy<Ny; iy++)
	  {
	    k += 1;
	    funcputl[k] = field(ix, iy);
	    funcputr[k] = field(Nx-xfdpts+ix, iy);
	  }

      k = -1;
      for(int ix=0; ix<xfdpts; ix++)
	for(int iy=0; iy<Ny; iy++)
	  {
	    k += 1;
	    if(ipx == 0)
	      funcgetl[k] = Zero;
	    if(ipx == maxproc1dx)
	      funcgetr[k] = Zero;
	  }
      
      // Send to the right
      if(ipx != maxproc1dx)
	{
	  int iproto = mycomm->get_iparray()[ipx+1][ipy];
	  mycomm->send_data(funcputr, length, iproto);
	}
      
      // Get from left
      if(ipx != 0)
      	{
      	  int iprofrom = mycomm->get_iparray()[ipx-1][ipy];	  
	  mycomm->get_data(funcgetl, length, iprofrom);
      	}
      
      // Send to left
      if(ipx != 0)
      	{
      	  int iprofrom = mycomm->get_iparray()[ipx-1][ipy];
      	  mycomm->send_data(funcputl, length, iprofrom);
      	}
      
      // Get from right.
      
      if(ipx != maxproc1dx)
      	{
      	  int iprofrom = mycomm->get_iparray()[ipx+1][ipy];
	  
      	  mycomm->get_data(funcgetr, length, iprofrom);
      	}
      
      k = -1;
      for(int ix=0; ix<xfdpts; ix++)
      	for(int iy=0; iy<Ny; iy++)
      	  {
      	    k += 1;
      	    halofield(-xfdpts+ix, iy) = funcgetl[k];
      	    halofield(Nx+ix, iy) = funcgetl[k];
      	  } 
      
      delete[] funcgetl;
      delete[] funcgetr;
      delete[] funcputl;
      delete[] funcputr;
      
    }

  //
  // Do communication in y
  
  if(numproc1dy > 1)
    {
      int length = yfdpts * Nx;
      
      dcomplex* funcputl = new dcomplex[length];
      dcomplex* funcputr = new dcomplex[length];
      dcomplex* funcgetl = new dcomplex[length];
      dcomplex* funcgetr = new dcomplex[length];
      
      // Set up the send arrays
      k = -1;
      for(int ix=0; ix<Nx; ix++)
  	for(int iy=0; iy<yfdpts; iy++)
  	  {
  	    k += 1;
  	    funcputl[k] = field(ix, iy);
  	    funcputr[k] = field(ix, Ny-yfdpts+iy);
  	  }
      
      k = -1;
      for(int ix=0; ix<Nx; ix++)
  	for(int iy=0; iy<yfdpts; iy++)
  	  {
  	    k += 1;
  	    if(ipy == 0)
  	      funcgetl[k] = Zero;
  	    if(ipy == maxproc1dy)
  	      funcgetr[k] = Zero;
  	  }
      
      // Send to the right
      if(ipy != maxproc1dy)
  	{
  	  int iproto = mycomm->get_iparray()[ipx][ipy+1]; 
  	  mycomm->send_data(funcputr, length, iproto);
  	}
      
      // Get form the left
      if(ipy != 0)
  	{
  	  int iproto = mycomm->get_iparray()[ipx][ipy-1];
	  
  	  mycomm->get_data(funcgetl, length, iproto);
  	}
      
      // Send to the left
      if(ipy != 0)
  	{
  	  int iprofrom = mycomm->get_iparray()[ipx][ipy-1];
  	  mycomm->send_data(funcputl, length, iprofrom);
  	}
      
      // Get it from the right.
      
      if(ipy != maxproc1dy)
  	{
  	  int iprofrom = mycomm->get_iparray()[ipx][ipy+1];
	  
  	  mycomm->get_data(funcgetr, length, iprofrom);
  	}
      
      k = -1;
      for(int ix=0; ix<Nx; ix++)
  	for(int iy=0; iy<yfdpts; iy++)
  	  {
  	    k += 1;
  	    halofield(ix, -yfdpts+iy) = funcgetl[k];
  	    halofield(ix, Ny+iy) = funcgetl[k];
  	  } 
      
      delete[] funcgetl;
      delete[] funcgetr;
      delete[] funcputl;
      delete[] funcputr;
      
    }
  
}
