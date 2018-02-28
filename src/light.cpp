////////////
////  light.cpp
////
////  Contains the class functions for the light class
////////////
#include "light.hpp"

light::light(fiber& grid, communicator& comm,
	     double _wavel,int argc,char **argv):
  mygrid{&grid}, mycomm{&comm}, wavelength{_wavel},
  field(grid.get_Nx(),grid.get_Ny(),arma::fill::zeros)
{
  // Asign value to private members
  Nx            = mygrid->get_Nx();
  Ny            = mygrid->get_Ny();
  Nxglobal      = Nx * mycomm->get_numproc1dx();
  Nyglobal      = Ny * mycomm->get_numproc1dy();
  dx            = mygrid->get_dx();
  dy            = mygrid->get_dy();  
  x_ax          = mygrid->get_x_ax();
  y_ax          = mygrid->get_y_ax();
  
  k0 = twopi / wavelength;
  
  //
  // Set up PETSC and SLEPC stuff for eigencalculation
  //
  PetscInt       Nlocal  = Nx * Ny;
  PetscInt       Nglobal = Nxglobal * Nyglobal;
  PetscInt       xrulepts = mygrid->get_xrulepts();
  PetscInt       yrulepts = mygrid->get_yrulepts();
  PetscInt       xfdpts = (xrulepts - 1)*0.5;
  PetscInt       yfdpts = (yrulepts - 1)*0.5;
  PetscErrorCode ierr;
  
  SlepcInitialize(&argc,&argv,(char*)0,NULL);
  
#if !defined(PETSC_USE_COMPLEX)
  SETERRQ(PETSC_COMM_WORLD,1,"This example requires complex numbers");
#endif

  //
  // -- Create PETSC matrix ---
  
  // Create the matrix that defines the eigensystem
  ierr = MatCreate(PETSC_COMM_WORLD,&Hmat);
  //ierr = MatSetSizes(Hmat,PETSC_DECIDE,PETSC_DECIDE,Nglobal,Nglobal);
  ierr = MatSetSizes(Hmat,Nlocal,Nlocal,Nglobal,Nglobal);
  ierr = MatSetType(Hmat, MATAIJ);
  
  ierr = MatSeqAIJSetPreallocation(Hmat, xrulepts+yrulepts, NULL);
  ierr = MatMPIAIJSetPreallocation(Hmat, xrulepts + yrulepts, NULL, yrulepts, NULL);
  ierr = MatSetFromOptions(Hmat);
  //ierr = MatSetOption(Hmat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(Hmat);
  
  // Create eigenvectors
  ierr = MatCreateVecs(Hmat,NULL,&xr);
  ierr = MatCreateVecs(Hmat,NULL,&xi);

  // // Create Eigensolver object
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);  
}

//////////////////////////////////////////////////

double light::get_norm()
{
  double norm {0.0};
  double normproc {0.0};
  
  for(int iy=0; iy<Ny; ++iy)
    for(int ix=0; ix<Nx; ++ix)
      normproc += real(conj(field(ix, iy))
		       * field(ix, iy));
  
  norm = mycomm->sumelements(normproc);
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

arma::cx_mat light::apply_laplacian(halo_cx_mat& halofunc)
{
  int xfdpts = (mygrid->get_xrulepts() - 1) * 0.5;
  int yfdpts = (mygrid->get_yrulepts() - 1) * 0.5;
  arma::cx_mat lapfield(Nx,Ny,arma::fill::zeros);
  
  // Apply 2nd derivative along the x axis
  for(int ifd=-xfdpts; ifd<xfdpts+1; ifd++)
    for(int iy=0; iy<Ny; iy++)
      for(int ix=0; ix<Nx; ix++)
  	lapfield(ix,iy) += mygrid->get_xfdcoeffs(ifd,2) * halofunc(ix+ifd,iy);
  
  // Apply 2nd derivative along the y axis
  for(int iy=0; iy<Ny; iy++)
    for(int ifd=-yfdpts; ifd<yfdpts+1; ifd++)
      for(int ix=0; ix<Nx; ix++)
	lapfield(ix,iy) += mygrid->get_yfdcoeffs(ifd,2) * halofunc(ix,iy+ifd);
  
  return lapfield;
}

///////////////////////////////////////////////////////////////////

void light::apply_helmholtz()
{
  // Create field matrix with halo points
  halo_cx_mat halofield(Nx,Ny,
			mygrid->get_xrulepts(),mygrid->get_yrulepts());
  // and set it to zero
  halofield.zeros();
  
  arma::mat n0 = mygrid->get_nfield();
  arma::mat n1 = mygrid->get_imag_nfield();
  //arma::mat n0(mygrid->get_nfield().memptr(), Nx, Ny, false, true);
  //arma::mat n1(mygrid->get_imag_nfield().memptr(), Nx, Ny, false, true); 
  
  // Copy data from field to halo halofield matrices
  arma::cx_mat lapfield(Nx,Ny,arma::fill::zeros);
  
  communicate_halo_points(halofield);
  
  // Apply Laplacian operator
  lapfield = apply_laplacian(halofield);
  
  // Apply potential
  for(int iy=0; iy<Ny; iy++)
    for(int ix=0; ix<Nx; ix++)
      lapfield(ix, iy) += n0(ix, iy) * n0(ix, iy)
	* k0 * k0 * field(ix, iy);
  
  if(mygrid->get_absorption_state())
    for(int iy=0; iy<Ny; iy++)
      for(int ix=0; ix<Nx; ix++)
	lapfield(ix, iy) -= Im * n1(ix, iy) * n1(ix, iy)
	  * k0 * k0 * field(ix, iy);
  
  field = lapfield;
  
}

///////////////////////////////////////////////////////////////////

int light::solve_helmholtz_eigenproblem(int argc,char **argv, int num_eigen_modes)
{
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki, *eigvec;
  PetscScalar    matelem;
  PetscInt       Nlocal  = Nx * Ny;
  PetscInt       Nglobal = Nxglobal * Nyglobal;
  PetscInt       Istart,Iend,nev,maxit,its,nconv;
  PetscInt       II, JJ, nrow, ncol;
  PetscInt       LocalVecSize;
  PetscInt       xrulepts = mygrid->get_xrulepts();
  PetscInt       yrulepts = mygrid->get_yrulepts();
  PetscInt       xfdpts = (xrulepts - 1)*0.5;
  PetscInt       yfdpts = (yrulepts - 1)*0.5;
  int            maxproc1dx = mycomm->get_maxproc1dx();
  int            maxproc1dy = mycomm->get_maxproc1dy();
  int            iproc      = mycomm->get_iprocessor();
  int            ipx        = mycomm->get_ipx();
  int            ipy        = mycomm->get_ipy();
  PetscErrorCode ierr;
  
  //-----------------------------------------------------
  
  ierr = MatGetOwnershipRange(Hmat,&Istart,&Iend);CHKERRQ(ierr);
  ierr = MatZeroEntries(Hmat);
  ierr = VecZeroEntries(xr);
  ierr = VecZeroEntries(xi);
  
  // Create matrix column vector,
  arma::cx_mat colvec(Nx,Ny,arma::fill::zeros);
  // and vectors containing indexes
  arma::uvec indxpair(2);
  arma::uvec indxvec(1);
  
  // Create matrix with halo points
  halo_cx_mat halovec(Nx,Ny,
		      mygrid->get_xrulepts(),mygrid->get_yrulepts());
  
  for(int iyy=-yfdpts; iyy<Ny+yfdpts; iyy++)
    for(int ixx=-xfdpts; ixx<Nx+xfdpts; ixx++)
      {
  	halovec.zeros();
  	if(ipx==0 && ixx<0)
  	  halovec(ixx, iyy) = Zero;
  	else if(ipy==0 && iyy<0)
  	  halovec(ixx, iyy) = Zero;
  	else if(ipx==maxproc1dx && ixx>=Nx)
  	  halovec(ixx, iyy) = Zero;
  	else if(ipy==maxproc1dy && iyy>=Ny)
  	  halovec(ixx, iyy) = Zero;
  	else
  	  halovec(ixx, iyy) = One;
	
  	ncol = Nx * iyy + ixx;
  	colvec.zeros();
	
  	// // Apply 2nd derivative along the x axis
  	for(int iy=0; iy<Ny; iy++)
  	  for(int ix=0; ix<Nx; ix++)
  	    for(int ifd=-xfdpts; ifd<xfdpts+1; ifd++)
  	      colvec(ix,iy) += mygrid->get_xfdcoeffs(ifd,2) * halovec(ix+ifd,iy);
	
  	// Apply 2nd derivative along the y axis
  	for(int iy=0; iy<Ny; iy++)
  	  for(int ifd=-yfdpts; ifd<yfdpts+1; ifd++)
  	    for(int ix=0; ix<Nx; ix++)
  	      colvec(ix,iy) += mygrid->get_yfdcoeffs(ifd,2) * halovec(ix, iy+ifd);
	
  	// Apply potential
  	for(int iy=0; iy<Ny; iy++)
  	  for(int ix=0; ix<Nx; ix++)
  	    colvec(ix, iy) += mygrid->get_nfield(ix, iy) * mygrid->get_nfield(ix, iy) *
  	      k0 * k0 * halovec(ix, iy);
	
  	//colvec = - colvec;
  	arma::uvec indxvec = arma::find(colvec);
  	indxvec.set_size(size(arma::find(colvec)));
  	indxvec = arma::find(colvec);
	
  	for(int k=0; k<indxvec.size(); k++)
  	  {
  	    nrow  = indxvec(k);
  	    indxpair = arma::ind2sub( arma::size(colvec), nrow );
  	    matelem = colvec(indxpair(0),indxpair(1));
  	    II = iproc * Nlocal + nrow;
  	    JJ = iproc * Nlocal + ncol;
	    ierr = MatSetValues(Hmat, 1, &II, 1, &JJ, &matelem, INSERT_VALUES);
  	  }
	
      }
  
  ierr = MatAssemblyBegin(Hmat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Hmat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //ierr = MatSetOption(Hmat,MAT_SYMMETRIC, PETSC_TRUE);CHKERRQ(ierr);
  
  // Create the eigensolver and set various options
  /*
    Create eigensolver context
  */
  // Set operators. In this case, it is a standard eigenvalue problem
  ierr = EPSSetOperators(eps,Hmat,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  //ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);
  //ierr = EPSSetType(eps,EPSARNOLDI);
  ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);
  ierr = EPSSetTolerances(eps,1e-8,PETSC_DEFAULT);
  ierr = EPSSetDimensions(eps,num_eigen_modes,PETSC_DEFAULT,PETSC_DEFAULT);
  
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
  
  // Open to file to save eigenenergies of the modes
  ofstream energyfile;
  string filename = "energy_modes_lambda_" + to_string(wavelength) + ".dat";
  if (mycomm->get_iprocessor()==0)
    energyfile.open(filename);
  
  if (nconv>0) {
    /*
      Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
		       "        beta^2              beta         ||Ax-kx||/||kx||\n"
		       "   ----------------- ------------------ ------------------\n");CHKERRQ(ierr);
    
    for (int i=0;i<nconv;i++) {
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
        ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %9f%+9fi %12g\n",(double(re)),(double(im)),sqrt(double(re)),sqrt(double(im)),
			   (double)error);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12f       %12g\n",(double(re)),sqrt(double(re)),(double)error);CHKERRQ(ierr);
      }
      // Save eigen energies to file
      if (mycomm->get_iprocessor()==0)
      	energyfile << i << " " << setw(8) << double(re) << " " << double(im)
		   << " " << sqrt(double(re)) << " " << sqrt(double(im)) << endl; 
    }
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }
  
  // Save eigenvectors to file
  PetscInt selected_mode = 0;
  ierr = EPSGetEigenpair(eps,selected_mode,&kr,&ki,xr,xi);CHKERRQ(ierr);
  ofstream modefile;
  filename = "field_dist/mode_no." + to_string(selected_mode)
    + ".wavelength." + to_string(wavelength) + "."
    + to_string(mycomm->get_iprocessor()) + ".dat";
  modefile.open(filename);
  
  // Get eigenvector's values from Petsc Vec type
  ierr = VecGetLocalSize(xr, &LocalVecSize);
  if(LocalVecSize!=Nlocal) 
    mycomm->parallel_stop("Local size of eigen vec not equal as local grid size!");
  
  ierr = VecGetArray(xr, &eigvec);
  for(int ii=0; ii<Nlocal; ii++)    
    modefile << PetscRealPart(eigvec[ii]) << " "
	     << PetscImaginaryPart(eigvec[ii]) << endl;
  
  // Close file
  energyfile.close();
  modefile.close();
  
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
      for(int iy=0; iy<Ny; iy++)
	for(int ix=0; ix<xfdpts; ix++)
	  {
	    k += 1;
	    funcputl[k] = field(ix, iy);
	    funcputr[k] = field(Nx-xfdpts+ix, iy);
	  }
      
      k = -1;
      for(int iy=0; iy<Ny; iy++)
	for(int ix=0; ix<xfdpts; ix++)
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
      for(int iy=0; iy<Ny; iy++)
	for(int ix=0; ix<xfdpts; ix++)
      	  {
      	    k += 1;
	    //halofield(-xfdpts+ix, iy) = funcgetl[k];
      	    //halofield(Nx+ix, iy) = funcgetr[k];
	    halofield.set_value(-xfdpts+ix,iy,funcgetl[k]);
	    halofield.set_value(Nx+ix,iy,funcgetr[k]);
	    
	    // if(ipx=3)
	    //   {
	    // 	cout << funcgetl[k] << endl;
	    // 	cout << "halo: " << halofield(-xfdpts+ix,iy) << endl;
	    //   }
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
      for(int iy=0; iy<yfdpts; iy++)
	for(int ix=0; ix<Nx; ix++)
  	  {
  	    k += 1;
  	    funcputl[k] = field(ix, iy);
  	    funcputr[k] = field(ix, Ny-yfdpts+iy);
  	  }
      
      k = -1;
      for(int iy=0; iy<yfdpts; iy++)
	for(int ix=0; ix<Nx; ix++)
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
      for(int iy=0; iy<yfdpts; iy++)
	for(int ix=0; ix<Nx; ix++)
  	  {
  	    k += 1;
  	    //halofield(ix, -yfdpts+iy) = funcgetl[k];
  	    //halofield(ix, Ny+iy) = funcgetr[k];
  	    halofield.set_value(ix, -yfdpts+iy, funcgetl[k]);
  	    halofield.set_value(ix, Ny+iy, funcgetr[k]);
	    
	    // if(ipx=3)
	    //   {
	    // 	cout << funcgetl[k] << endl;
	    // 	cout << "halo: " << halofield(ix,-yfdpts+iy) << endl;
	    //   }
	    
	  } 
      
      delete[] funcgetl;
      delete[] funcgetr;
      delete[] funcputl;
      delete[] funcputr;
      
    }
  
}
