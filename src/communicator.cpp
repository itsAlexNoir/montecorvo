////////////
////  communicator.cpp
////
////  Contains the class functions for the communicator class
////////////
#include <iostream>
#include <math.h>
#include "communicator.hpp"


communicator::communicator(int argc, char **argv,
			   const int numprocessorsx,
			   const int numprocessorsy):
  numproc1dx{numprocessorsx}, numproc1dy{numprocessorsy}
{

  int ipro, mpi_size, ierr;
  
  // Start up MPI
  ierr = MPI_Init(&argc, &argv);
  
  // Find out number of processors.
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Check if the number of processors requested are OK for the calculation
  if(mpi_size == numprocessors)
    {
      string message = "Number of processors requested for calculation (N="
        + to_string(mpi_size)
	+ ") does not equal the number of processors specified in the input file (N = "
	+ to_string(numprocessors) + ") ";
      parallel_stop(message);
    }
  
  // Find out number of the processor we are working on.
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iprocessor);
  
  
  numprocessors = numproc1dx * numproc1dy; 
  maxproc1dx    = numproc1dx - 1;
  maxproc1dy    = numproc1dy - 1;

  iparray = new int*[numproc1dx];
  for(int iprocx=0; iprocx<numproc1dx; iprocx++)
    iparray[iprocx] = new int[numproc1dy];
  
  iparrayx = new int[numproc1dx];
  iparrayy = new int[numproc1dy];
  
  for(int iprocx=0; iprocx<numproc1dx; iprocx++)
    for(int iprocy=0; iprocy<numproc1dy; iprocy++)
      iparray[iprocx][iprocy] = 100.0;
  
  for(int iprocx=0; iprocx<numproc1dx; iprocx++)
    iparrayx[iprocx] = 100.0;
  
  for(int iprocy=0; iprocy<numproc1dy; iprocy++)
    iparrayy[iprocy] = 100.0;
  
  ipro = -1;
  for(int iprocx=0; iprocx<numproc1dx; iprocx++)
    for(int iprocy=0; iprocy<numproc1dy; iprocy++)
      {
	    ipro += 1;
	    
	    iparray[iprocx][iprocy] = ipro;
            iparrayx[ipro]          = iprocx;
	    iparrayy[ipro]          = iprocy;
	    
	    if(iprocessor == ipro)
	      {
		ipx = iprocx;
		ipy = iprocy;
	      }
      }
  
}

///////////////////////////////////////////////////////////////////

void communicator::print_communicator_parameters()
{
  cout << endl;
  cout << "----- Communication parameters ---------------------" << endl;
  cout << " Number of processors in x: " << numproc1dx << endl;
  cout << " Number of processors in y: " << numproc1dy << endl;
  cout << " Number of processors:      " << numprocessors << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl << endl; 

}

///////////////////////////////////////////////////////////////////

void communicator::parallel_stop(const string message ) const
{    
  
  if(iprocessor==0)
    cout << message << endl;

  exit(1);
  
}    

///////////////////////////////////////////////////////////////////
    

void communicator::print_message(const string message) const
{

  // Wait for all processors
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(iprocessor == 0)
    cout << message << endl;
  
}  

///////////////////////////////////////////////////////////////////
//
/// Send a data block from current processor to processor iproto.
//

void communicator::send_data(const dcomplex* data, const int length,
			     const int iproto) const
{
  int ierr;
  int itag = iprocessor;
  
  ierr = MPI_Send(data, length, MPI_C_DOUBLE_COMPLEX,
		  iproto, itag, MPI_COMM_WORLD);
  
}

///////////////////////////////////////////////////////////////////
//
/// Get a data block from processor iprofrom to current processor.
//

void communicator::get_data(dcomplex* data, const int length,
			    const int iprofrom) const
{
  int ierr;
  MPI_Status status;
  int itag = iprofrom;
  
  ierr = MPI_Recv(data, length, MPI_C_DOUBLE_COMPLEX,
		  iprofrom, itag, MPI_COMM_WORLD, &status);
  
}

///////////////////////////////////////////////////////////////////
