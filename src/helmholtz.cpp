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
  
  communicator comm(argc, argv, 3, 1);
  
  // Print starting date
  auto point_time = chrono::system_clock::now();
  auto start_date = chrono::system_clock::to_time_t(point_time);
  if(comm.get_iprocessor()==0)
    cout << endl << "Program starts at " << ctime(&start_date) << endl;
  
  int Nx{113};
  int Ny{335};
  
  int rulepts{5};
  
  double dx{0.3};
  double dy{0.3};
  
  double wavelength{0.8};
  
  // Initialise space
  space grid(Nx, Ny, dx, dy, rulepts, rulepts, comm);
  if(comm.get_iprocessor()==0)
    grid.print_grid_parameters();
  
  // Initialize electromagnetic field
  light pulse(grid, comm, wavelength);
  
  // Print starting date
  point_time = chrono::system_clock::now();
  start_date = chrono::system_clock::to_time_t(point_time);
  if(comm.get_iprocessor()==0)
    cout << endl << "Eigenproblem calculation starts at " << ctime(&start_date) << endl;
  
  // Create wavefunction and set it to the ground state
  auto t0 = chrono::high_resolution_clock::now();
  
  pulse.solve_helmholtz_eigenproblem(argc, argv);
  
  auto t1 = chrono::high_resolution_clock::now();
  
  // Measuring duration of atom creation
  chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds>(t1-t0);
  chrono::seconds sec       = chrono::duration_cast<chrono::seconds>(ms);
  chrono::milliseconds msec = chrono::duration_cast<chrono::milliseconds>(ms % chrono::seconds(1) );
  chrono::nanoseconds nsec = chrono::duration_cast<chrono::nanoseconds>(ms % chrono::milliseconds(1) );
  
  if(comm.get_iprocessor()==0)
    cout << "Initial EM field set in: " << setfill('0') << setw(2) << sec.count() << " seconds, "
	 << setw(3) << msec.count() << " milliseconds, "
	 << setw(3) << nsec.count() << " nanoseconds." << endl;
  
  // Print ending date
  point_time = chrono::system_clock::now();
  auto end_date = chrono::system_clock::to_time_t(point_time);
  if(comm.get_iprocessor()==0)
    cout << endl << "Eigenproblem ends at " << ctime(&end_date) << endl;
  
  // Print starting date
  // point_time = chrono::system_clock::now();
  // start_date = chrono::system_clock::to_time_t(point_time);
  // cout << endl << "Evolution starts at " << ctime(&start_date) << endl;

  // t0 = chrono::high_resolution_clock::now();
  
  // // Time evolve the atom
  // myhydrogen.evolve(pulse);

  // t1 = chrono::high_resolution_clock::now();
  
  // // Measuring Time
  // ms = chrono::duration_cast<chrono::milliseconds>(t1-t0);
  // chrono::hours hh    = chrono::duration_cast<chrono::hours>(ms);
  // chrono::minutes min = chrono::duration_cast<chrono::minutes>(ms % chrono::hours(1));
  // sec                 = chrono::duration_cast<chrono::seconds>(ms % chrono::minutes(1) );
  // msec                = chrono::duration_cast<chrono::milliseconds>(ms % chrono::seconds(1) );
  
  // cout << "Execution took: " << setfill('0') << setw(2) << hh.count() << " hours, "
  //      << setw(2) << min.count() << " minutes, "
  //      << setw(2) << sec.count() << " seconds, " << setw(3) << msec.count() << " milliseconds."
  //      << endl;
								       
  // Print ending date
  // point_time = chrono::system_clock::now();
  // end_date = chrono::system_clock::to_time_t(point_time);
  // cout << endl << "Program ends at " << ctime(&end_date) << endl;
  
}
    
