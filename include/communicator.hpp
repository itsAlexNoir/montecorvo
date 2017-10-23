////////////
////  communicator.hpp
////
//// Declarations of the communicator class
////////////

#ifndef ____COMMUNICATOR__
#define ____COMMUNICATOR__

#include <iostream>
#include <math.h>
#include "mpi.h"
#include "constants.hpp"


class communicator{
  
private:

  int iprocessor;
  int iprocglobal;

  int numproc1dx;
  int numproc1dy;
  int numproc1dz;
  int maxproc1dx;
  int maxproc1dy;
  int maxproc1dz;

  int numprocessors;
  int ipx, ipy, ipz;
  int **iparray;
  int *iparrayx;
  int *iparrayy;
  int *iparrayz;
  
  
  
public:
  communicator(int argc, char **argv,
	       const int numprocessorsx,
	       const int numprocessorsy);
  ~communicator()
  { delete[] iparrayx;
    delete[] iparrayy;
    for(int iprocx=0; iprocx<numproc1dx; iprocx++)
      delete[] iparray[iprocx];
    delete[] iparray;  
    MPI_Finalize();}

  // Set / Get routines
  int get_iprocessor() const {return iprocessor;}
  int get_numproc1dx() const {return numproc1dx;}
  int get_numproc1dy() const {return numproc1dy;}
  int get_maxproc1dx() const {return maxproc1dx;}
  int get_maxproc1dy() const {return maxproc1dy;}

  int get_ipx() const {return ipx;} 
  int get_ipy() const {return ipy;}
  
  int** get_iparray() const {return iparray;}
  int* get_iparrayx() const {return iparrayy;}
  int* get_iparrayy() const {return iparrayy;}
  
  void print_communications();
  
  // Print out communicator parameters on the screen
  void print_communicator_parameters();
  
  // Print message and exit program with error
  void parallel_stop(const string message) const;
  
  // Just print a message to the screen
  void print_message(const string message) const;
  
  // Passing-message routines
  void send_data(const dcomplex* data, const int length,
		 const int iproto) const;
  
  void get_data(dcomplex* data, const int length,
		const int iprofrom) const;
 
};


#endif // ____COMMUNICATOR_ //
