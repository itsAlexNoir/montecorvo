////////////////////////////////////////////
///////////      Helmholtz.cpp       ///////
//////////      Main program         ///////
////////////////////////////////////////////

#include "constants.hpp"
#include "communicator.hpp"
//#include "params.hpp"
#include "space.hpp"
#include "light.hpp"

int main(int argc,char **argv){

  ///////////////////////////////////////////////////////////////////////
  // Initialise communications between processors
  ///////////////////////////////////////////////////////////////////////
  
  communicator comm(argc, argv, 1, 15);
  
  // Print starting date
  auto point_time = chrono::system_clock::now();
  auto start_date = chrono::system_clock::to_time_t(point_time);
  if(comm.get_iprocessor()==0)
    cout << endl << "Program starts at " << ctime(&start_date) << endl;
  
  ///////////////////////////////////////////////////////////////////////
  // Set main calculation parameters
  ///////////////////////////////////////////////////////////////////////

  //int Nx{401};
  // int Nx{35};
  // int Ny {501};
  int Ny{35};
  int Nx {501};
  
  //int Nx{27};
  //int Ny{301};
  
  int rulepts{5};
  
  double dx{0.2};
  double dy{0.2};

  double wavelength0      {0.8};
  double wavelength       {0.0};
  int no_wavelength       {1};
  double delta_wavelength {1e-3};
  int num_modes           {4};
  
  string shots_filename   {"5capas-core1-codigo.txt"};
  double n0               {0.0};
  double deltan           {-0.004};
  int Nshx                {11};
  int Nshy                {11};
  double delta_shx        {6.0};
  double delta_shy        {5.0};
  double sigma_shx        {0.5};
  double sigma_shy        {2.5};
  
  // int nholes              {15};
  // double holes_radius     {15.0};
  
  ///////////////////////////////////////////////////////////////////////
  // Set fiber and pulse light objects
  ///////////////////////////////////////////////////////////////////////

  // Print out communication
  if(comm.get_iprocessor()==0)  
    comm.print_communicator_parameters();
  
  // Initialise space
  space grid(Nx, Ny, dx, dy, rulepts, rulepts, comm);
  if(comm.get_iprocessor()==0)
    grid.print_grid_parameters();

  // Initialize electromagnetic field
  light pulse(grid, comm, wavelength0, argc, argv);
  
  ///////////////////////////////////////////////////////////////////////
  // Start calculation propagation constant (beta)
  ///////////////////////////////////////////////////////////////////////
  
  // Print starting date
  point_time = chrono::system_clock::now();
  start_date = chrono::system_clock::to_time_t(point_time);
  if(comm.get_iprocessor()==0)
    cout << endl << "Eigenproblem calculation starts at " << ctime(&start_date) << endl;
  
  auto t0 = chrono::high_resolution_clock::now();
  
  // Start loop over wavelengths
  for(int iw=0; iw<no_wavelength; iw++)
    {
      wavelength = wavelength0 + double(iw) * delta_wavelength;
      pulse.set_wavelength(wavelength);
      // Set fiber's structure
      n0 = grid.get_refractive_index("YAG", wavelength);
      
      //grid.set_step_index_fiber(4.0, 1.445, 1.4378);
      
      // grid.set_circular_honeycomb_fiber(holes_radius, nholes,
      // 				    n0, deltan,
      // 				    delta_holex,
      // 				    delta_holey);
      
      grid.set_honeycomb_fiber(shots_filename, n0, deltan, Nshx, Nshy,
			  delta_shx, delta_shy, sigma_shx, sigma_shy, 8);
      
      pulse.solve_helmholtz_eigenproblem(argc, argv, num_modes);
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
    
