/*!
  @file spectralscan.cpp
  \brief Program Performing an scan of the 
  propagation constant over a spectral range.
*/
#include "constants.hpp"
#include "communicator.hpp"
#include "io.hpp"
#include "params.hpp"
#include "fiber.hpp"
#include "light.hpp"

int main(int argc,char **argv){

  /// Read parameters for this calculation
  // if(argc<1)
  //   {
  //     cout << "You must provide an input file to the program." << endl;
  //     cout << "Aborting program..." << endl;
  //   } 
  //input myinput(argv[1]);
  input myinput("./input.inp");
  
  int nprocsx = myinput.read_integer("nprocsx", 1);
  int nprocsy = myinput.read_integer("nprocsy", 1);
  
  ///////////////////////////////////////////////////////////////////////
  // Initialise communications between processors
  ///////////////////////////////////////////////////////////////////////
  
  communicator comm(argc, argv, nprocsx, nprocsy);
  
  // Print starting date
  auto point_time = chrono::system_clock::now();
  auto start_date = chrono::system_clock::to_time_t(point_time);
  if(comm.get_iprocessor()==0)
    cout << endl << "Program starts at " << ctime(&start_date) << endl;
  
  ///////////////////////////////////////////////////////////////////////
  // Set main calculation parameters
  ///////////////////////////////////////////////////////////////////////
  
  int Nx                  = myinput.read_integer("Nx", 100);
  int Ny                  = myinput.read_integer("Ny", 100);
  int rulepts             = myinput.read_integer("rulepts", 5);
  
  double dx               = myinput.read_double("dx", 0.2);
  double dy               = myinput.read_double("dy", 0.2);
  
  double wavelength       {0.0};
  double frequency        {0.0};
  double nw               {0.0};
  double frequency0       {0.0};
  bool frequency_scan     = myinput.read_boolean("frequency_scan", false);
  double wavelength0      = myinput.read_double("wavelength0", 0.8);
  double wavelength1      = myinput.read_double("wavelength1", 0.85);
  double no_frequency     = myinput.read_double("no_frequency", 1);
  double no_wavelength    = myinput.read_double("no_wavelength", 1);
  int num_modes           = myinput.read_integer("num_modes", 1);
  int num_saved_modes     = myinput.read_integer("num_saved_modes", 1);
  bool cladding_on        = myinput.read_boolean("cladding_on", false);
  double rclad            = myinput.read_double("cladding_radius", 10.0);
  bool abs_switch         = myinput.read_boolean("abs_state", false);
  int abs_mask_exponent   = myinput.read_integer("abs_mask_exponent",2);
  double deltan_abs       = myinput.read_double("dn_abs", -0.002);

  string fiber_structure  = myinput.read_string("fiber_structure","step_index");
  string shots_filename   = myinput.read_string("shots_filename", "shots");
  double n0, n1           {0.0};
  double deltan           = myinput.read_double("deltan", 1e-4);
  int no_sides            = myinput.read_integer("no_sides",1);
  double hole_radius      = myinput.read_double("hole_radius",1.0);
  int Nshx                = myinput.read_integer("Nshx", 11);
  int Nshy                = myinput.read_integer("Nshx", 11);
  double delta_shx        = myinput.read_double("delta_shx", 3.0);
  double delta_shy        = myinput.read_double("delta_shy", 6.0);
  double sigma_shx        = myinput.read_double("sigma_shx", 0.5);
  double sigma_shy        = myinput.read_double("sigma_shy", 2.5);
  
  bool save_to_hdf5;
  arma::uvec selected_modes(num_saved_modes);
  if(num_saved_modes>0)
    {
      save_to_hdf5 = true;
      selected_modes = arma::regspace<arma::uvec>(0,num_saved_modes-1);
    }
  else
    save_to_hdf5 = false;
  
  ///////////////////////////////////////////////////////////////////////
  // Set fiber and pulse light objects
  ///////////////////////////////////////////////////////////////////////

  // Print out communication
  if(comm.get_iprocessor()==0)  
    comm.print_communicator_parameters();
  
  // Initialise space
  fiber grid(Nx, Ny, dx, dy, rulepts, rulepts, abs_switch, comm);
  if(comm.get_iprocessor()==0)
    grid.print_grid_parameters();

  // Initialize electromagnetic field
  light pulse(grid, comm, wavelength0, argc, argv);
  
  ///////////////////////////////////////////////////////////////////////
  // Start calculation propagation constant (beta)
  ///////////////////////////////////////////////////////////////////////
  
  // Prepare loop for spectral scan
  double delta_frequency {0.0};
  double delta_wavelength {0.0};
  if(frequency_scan)
    {
      double lower_frequency = pulse.get_frequency(wavelength0);
      double upper_frequency = pulse.get_frequency(wavelength1);
      delta_frequency = fabs(upper_frequency - lower_frequency) / no_frequency;
      if(lower_frequency<=upper_frequency)
	frequency0 = lower_frequency;
      else
	frequency0 = upper_frequency;
      nw = no_frequency;
    }
  else
    {
      delta_wavelength = fabs(wavelength1 - wavelength0) / no_wavelength;
      nw = no_wavelength;
    }
  
  // Print starting date
  point_time = chrono::system_clock::now();
  start_date = chrono::system_clock::to_time_t(point_time);
  if(comm.get_iprocessor()==0)
    cout << endl << "Eigenproblem calculation starts at " << ctime(&start_date) << endl;
  
  auto t0 = chrono::high_resolution_clock::now();
  
  // Start loop over wavelengths
  for(int iw=0; iw<nw; iw++)
    {
      if(frequency_scan)
	{
	  frequency = frequency0 + double(iw) * delta_frequency;
	  wavelength = pulse.get_wavelength(frequency);
	}
      else
	wavelength = wavelength0 + double(iw) * delta_wavelength; 
      
      pulse.set_wavelength(wavelength);
      // Set fiber's structure
      n0 = grid.get_refractive_index("YAG", wavelength);
      n1 = n0 + deltan_abs;
      
      if(fiber_structure=="regular_honeycomb")
	grid.set_honeycomb_fiber(shots_filename, n0, deltan, Nshx, Nshy,
				 delta_shx, delta_shy, sigma_shx, sigma_shy, 8);
      else if(fiber_structure=="circular_honeycomb")
	grid.set_circular_honeycomb_fiber(hole_radius, no_sides,
					  1.0, 1.6,
					  n0, deltan,
					  sigma_shx,
					  sigma_shy, 8);
      else if(fiber_structure=="step_index")
	grid.set_step_index_fiber(rclad, n0, n1);
      else
	cout << "You must choose an structure for the fiber." << endl;
      
      if(cladding_on)
	grid.set_fiber_cladding(n1, rclad);
      
      if(abs_switch)
	grid.set_fiber_absorption(n1, rclad, abs_mask_exponent);
      
      // Save fiber profile to HDF5 file
      string filename = "fiber_profile.n0." + to_string(n0);
      save2DMatrix_hdf5(filename, grid.get_nfield(), &grid, &comm );
      
      // Solve!
      pulse.solve_helmholtz_eigenproblem(argc, argv, num_modes,
       					 save_to_hdf5, selected_modes);
    }
  
  auto t1 = chrono::high_resolution_clock::now();
  
  // Measuring duration of atom creation
  chrono::milliseconds ms   = chrono::duration_cast<chrono::milliseconds>(t1-t0);
  chrono::minutes min       = chrono::duration_cast<chrono::minutes>(ms);
  chrono::seconds sec       = chrono::duration_cast<chrono::seconds>(ms % chrono::minutes(1));
  chrono::milliseconds msec = chrono::duration_cast<chrono::milliseconds>(ms % chrono::seconds(1) );
  chrono::nanoseconds nsec  = chrono::duration_cast<chrono::nanoseconds>(ms % chrono::milliseconds(1) );
  
  // chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds>(t1-t0);
  // chrono::seconds sec       = chrono::duration_cast<chrono::seconds>(ms);
  // chrono::milliseconds msec = chrono::duration_cast<chrono::milliseconds>(ms % chrono::seconds(1) );
  // chrono::nanoseconds nsec = chrono::duration_cast<chrono::nanoseconds>(ms % chrono::milliseconds(1) );
  
  if(comm.get_iprocessor()==0)
    cout << "Initial EM field set in: " << setfill('0') << setw(2) << min.count() << " minutes, "
	 << setw(2) << sec.count() << " seconds, "
	 << setw(3) << msec.count() << " milliseconds, "
	 << setw(3) << nsec.count() << " nanoseconds." << endl;
  
  // Print ending date
  point_time = chrono::system_clock::now();
  auto end_date = chrono::system_clock::to_time_t(point_time);
  if(comm.get_iprocessor()==0)
    cout << endl << "Eigenproblem ends at " << ctime(&end_date) << endl;
  
  ///////////////////////////////////////////////////////////////////////
  //  End of the program 
  ///////////////////////////////////////////////////////////////////////

}
    
