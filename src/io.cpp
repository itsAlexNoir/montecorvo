/*!
  @file io.cpp
  @brief Functions for I/O of the program.
  
  This file contains routines to write to disk
  datasets of different kind.
*/

#include "io.hpp"

//////////////////////////////////////////////////////////
/*!
  @function save2DMatrix_hdf5
  @brief Stores 2D real double matrices to HFD5 files.
  @params[in] std::string name 
  @params[out] arma::mat matrix Real double data matrix
*/
//////////////////////////////////////////////////////////
void save2DMatrix_hdf5(string name, arma::mat mymatrix)
{
  string filename {name+".h5"};
  hid_t  file_id; //> HDF5 File identifier */
  hid_t dataset_id, dataspace_id;  
  herr_t status;  //> Variable that checks error in the I/O of the file.
  hsize_t dims[2]; //> Array with matrix dimensions
  
  dims[0] = 4; 
  dims[1] = 6; 
  dataspace_id = H5Screate_simple(2, dims, NULL);

   /* Create a new file using default properties. */
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  /* Create the data space for the dataset. */
  dims[0] = mymatrix.n_cols; 
  dims[1] = mymatrix.n_rows; 
  dataspace_id = H5Screate_simple(2, dims, NULL);
  
  /* Create the dataset. */
  string setname = "/"+name;
  dataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Write the dataset. */
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		    mymatrix.memptr());
   
  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);
  
  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);

  /* Terminate access to the file. */
  status = H5Fclose(file_id); 

}
//////////////////////////////////////////////////////////
/*!
  @function save2DMatrix_hdf5
  @brief Stores 2D complex matrices to HFD5 files.
  @params[in] std::string name 
  @params[out] arma::cx_mat matrix Complex data matrix
*/
//////////////////////////////////////////////////////////
void save2DMatrix_hdf5(string name, arma::cx_mat mymatrix)
{
  string filename {name+".h5"};
  hid_t  file_id; //> HDF5 File identifier */
  hid_t dataset_id, dataspace_id;  
  herr_t status;  //> Variable that checks error in the I/O of the file.
  hsize_t dims[3]; //> Array with matrix dimensions
  
  
  /* Create a new file using default properties. */
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  /* Create the data space for the dataset. */
  dims[0] = mymatrix.n_cols; 
  dims[1] = mymatrix.n_rows; 
  dims[2] = 2;
  dataspace_id = H5Screate_simple(3, dims, NULL);
  
  /* Create the dataset. */
  string setname = "/"+name;
  dataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  // Transform n-D complex matrix in (n+1)-D real matrix
  
  
  /* Write the dataset. */
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		    mymatrix.memptr());
   
  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);
  
  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);

  /* Terminate access to the file. */
  status = H5Fclose(file_id); 

}
//////////////////////////////////////////////////////////
/*!
  @function save2DMatrix_hdf5
  @brief Stores 2D real matrices that are split 
  across processors to parallel HFD5 files.
  @params[in] std::string name 
  @params[out] arma::cx_mat matrix Real double data matrix
*/
//////////////////////////////////////////////////////////
void save2DMatrix_hdf5(string name, arma::mat mymatrix, fiber* mygrid, communicator* comm)
{
  string filename {name+".h5"};
  hid_t  file_id; //> HDF5 File identifier */
  hid_t	plist_id; //> Property list identifier(for accessing template) */
  hid_t dataset_id, dataspace_id, memspace; //> Dataspace, dataset and memory identifiers
  hid_t ds_id; //> Dimension scale identifier
  herr_t status;  //> Variable that checks error in the I/O of the file.
  hsize_t dims[2], count[2], offset[2]; //> Array with matrix dimensions
  MPI_Info info  = MPI_INFO_NULL; //> Control variable for MPI communications

  // Set up file access property list with parallel I/O access
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm->get_mpicomm(), info);
  
  /* Create a new file using default properties. */
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      plist_id);
  status = H5Pclose(plist_id);
  
  /* Create the data space for the dataset. */
  //dims[0] = mymatrix.n_cols; 
  //dims[1] = mymatrix.n_rows;
  dims[0] = mygrid->get_Nyglobal();
  dims[1] = mygrid->get_Nxglobal();
  // Create the dataspace for the dataset.
  dataspace_id = H5Screate_simple(2, dims, NULL);
  
  /* Create the dataset. */
  string setname = "/"+name;
  dataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, 
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(dataspace_id);
  
  // Each process defines dataset in memory and writes it to the hyperslab
  // in the file.
  count[0] = mygrid->get_Ny();
  count[1] = mygrid->get_Nx();
  offset[0] = comm->get_ipy() * mygrid->get_Ny();
  offset[1] = comm->get_ipx() * mygrid->get_Nx();;
  memspace = H5Screate_simple(2, count, NULL);
  
  // Select hyperslab in the file.
  dataspace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  
  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  //// To write dataset independently use
  // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); 
  
  /* Write the dataset. */
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace, dataspace_id,
		    plist_id, mymatrix.memptr() );

  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);
  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);
  // Close property list.
  status = H5Pclose(plist_id);

  ///////////////////////////////////
  /*           Axes                */
  ///////////////////////////////////
  /* Y axis (the first dimension) */
  
  dataspace_id = H5Screate_simple(1, &dims[0], NULL);
  setname = "/Y_AXIS";
  dataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, 
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(dataspace_id);
  memspace = H5Screate_simple(1, &count[0], NULL);
  
  dataspace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, &offset[0], NULL, &count[0], NULL);
  
  
  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  //// To write dataset independently use
  // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); 
  
  /* Write the dataset. */
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace, dataspace_id,
		    plist_id, mygrid->get_y_ax().memptr() );
  
  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);
  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);
  // Close property list.
  status = H5Pclose(plist_id);

  
  /* X axis (the first dimension) */
  dataspace_id = H5Screate_simple(1, &dims[1], NULL);
  setname = "/X_AXIS";
  dataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, 
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(dataspace_id);
  memspace = H5Screate_simple(1, &count[1], NULL);
  
  dataspace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, &offset[1], NULL, &count[1], NULL);
  
  
  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  //// To write dataset independently use
  // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); 
  
  /* Write the dataset. */
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace, dataspace_id,
		    plist_id, mygrid->get_x_ax().memptr() );
  
  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(dataset_id);
  /* Terminate access to the data space. */ 
  status = H5Sclose(dataspace_id);
  // Close property list.
  status = H5Pclose(plist_id);

  /////////// Attach Dimension scale to the main dataset
  ///////
  /* Attach first dimension to the dataset */
  setname = "/"+name;
  dataset_id = H5Dopen2(file_id,setname.c_str(), H5P_DEFAULT);
  setname = "/X_AXIS";
  ds_id = H5Dopen2(file_id,setname.c_str(), H5P_DEFAULT);
  H5DSattach_scale(dataset_id,ds_id,0);
  /* Close ds_id after use it */
  status = H5Dclose(ds_id);
  
  /* Attach second dimension to the dataset */
  setname = "/Y_AXIS";
  ds_id = H5Dopen2(file_id,setname.c_str(), H5P_DEFAULT);
  H5DSattach_scale(dataset_id,ds_id,1);
  
  /* Close ds_id after use it, and also dataset_di */
  status = H5Dclose(ds_id);
  status = H5Dclose(dataset_id);
  
  /* Terminate access to the file. */
  status = H5Fclose(file_id); 
  
}

//////////////////////////////////////////////////////////
/*!
  @function save2DMatrix_hdf5
  @brief Stores 2D complex matrices that are split 
  across processors to parallel HFD5 files.
  @params[in] std::string name 
  @params[out] arma::cx_mat matrix Real double data matrix
*/
//////////////////////////////////////////////////////////
void save2DMatrix_hdf5(string name, arma::cx_mat mymatrix, fiber* mygrid, communicator* comm)
{
  string filename {name+".h5"};
  hid_t  file_id; //> HDF5 File identifier */
  hid_t	plist_id; //> Property list identifier(for accessing template) */
  hid_t realdataset_id, imagdataset_id; //> Datasets identifiers
  hid_t realdataspace_id, imagdataspace_id; //> Dataspace identifiers
  hid_t realmemspace, imagmemspace; //>  Memory identifiers
  hid_t ds_id; //> Dimension scale identifier
  herr_t status;  //> Variable that checks error in the I/O of the file.
  hsize_t dims[2], count[2], offset[2]; //> Array with matrix dimensions
  MPI_Info info  = MPI_INFO_NULL; //> Control variable for MPI communications
  
  // Set up file access property list with parallel I/O access
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm->get_mpicomm(), info);
  
  /* Create a new file using default properties. */
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      plist_id);
  status = H5Pclose(plist_id);
  
  /* Create the data space for the dataset. */
  dims[0] = mygrid->get_Nyglobal();
  dims[1] = mygrid->get_Nxglobal();
  // Create the dataspace for the dataset.
  realdataspace_id = H5Screate_simple(2, dims, NULL);
  imagdataspace_id = H5Screate_simple(2, dims, NULL);
  
  /* Create the dataset. */
  /* Real dataset */
  string setname = "/"+name+".real";
  realdataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, realdataspace_id, 
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(realdataspace_id);
  /* Imaginary dataset */
  setname = "/"+name+".imag";
  imagdataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, imagdataspace_id, 
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(imagdataspace_id);
  
  // Each process defines dataset in memory and writes it to the hyperslab
  // in the file.
  count[0] = mygrid->get_Ny();
  count[1] = mygrid->get_Nx();
  offset[0] = comm->get_ipy() * mygrid->get_Ny();
  offset[1] = comm->get_ipx() * mygrid->get_Nx();;
  realmemspace = H5Screate_simple(2, count, NULL);
  imagmemspace = H5Screate_simple(2, count, NULL);
  
  // Select hyperslab in the file.
  realdataspace_id = H5Dget_space(realdataset_id);
  H5Sselect_hyperslab(realdataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  imagdataspace_id = H5Dget_space(imagdataset_id);
  H5Sselect_hyperslab(imagdataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  
  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  //// To write dataset independently use
  // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); 
  
  /* Write the dataset. */
  arma::mat realmat(real(mymatrix));
  status = H5Dwrite(realdataset_id, H5T_NATIVE_DOUBLE, realmemspace, realdataspace_id,
		    plist_id, realmat.memptr());
  arma::mat imagmat(imag(mymatrix));
  status = H5Dwrite(imagdataset_id, H5T_NATIVE_DOUBLE, imagmemspace, imagdataspace_id,
		    plist_id, imagmat.memptr());
  
  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(realdataset_id);
  status = H5Dclose(imagdataset_id);
  /* Terminate access to the data space. */ 
  status = H5Sclose(realdataspace_id);
  status = H5Sclose(imagdataspace_id);
  status = H5Sclose(realmemspace);
  status = H5Sclose(imagmemspace);
  // Close property list.
  status = H5Pclose(plist_id);


  ///////////////////////////////////
  /*           Axes                */
  ///////////////////////////////////
  /* Y axis (the first dimension) */
  
  realdataspace_id = H5Screate_simple(1, &dims[0], NULL);
  setname = "/Y_AXIS";
  realdataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, realdataspace_id, 
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(realdataspace_id);
  realmemspace = H5Screate_simple(1, &count[0], NULL);
  
  realdataspace_id = H5Dget_space(realdataset_id);
  H5Sselect_hyperslab(realdataspace_id, H5S_SELECT_SET, &offset[0], NULL, &count[0], NULL);
  
  
  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  //// To write dataset independently use
  // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); 
  
  /* Write the dataset. */
  status = H5Dwrite(realdataset_id, H5T_NATIVE_DOUBLE, realmemspace, realdataspace_id,
		    plist_id, mygrid->get_y_ax().memptr() );
  
  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(realdataset_id);
  /* Terminate access to the data space. */ 
  status = H5Sclose(realdataspace_id);
  status = H5Sclose(realmemspace);
  // Close property list.
  status = H5Pclose(plist_id);
  
  
  /* X axis (the first dimension) */
  realdataspace_id = H5Screate_simple(1, &dims[1], NULL);
  setname = "/X_AXIS";
  realdataset_id = H5Dcreate2(file_id, setname.c_str(), H5T_NATIVE_DOUBLE, realdataspace_id, 
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(realdataspace_id);
  realmemspace = H5Screate_simple(1, &count[1], NULL);
  
  realdataspace_id = H5Dget_space(realdataset_id);
  H5Sselect_hyperslab(realdataspace_id, H5S_SELECT_SET, &offset[1], NULL, &count[1], NULL);
  
  
  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  //// To write dataset independently use
  // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); 
  
  /* Write the dataset. */
  status = H5Dwrite(realdataset_id, H5T_NATIVE_DOUBLE, realmemspace, realdataspace_id,
		    plist_id, mygrid->get_x_ax().memptr() );
  
  /* End access to the dataset and release resources used by it. */
  status = H5Dclose(realdataset_id);
  /* Terminate access to the data space. */ 
  status = H5Sclose(realdataspace_id);
  status = H5Sclose(realmemspace);
  // Close property list.
  status = H5Pclose(plist_id);

  /////////// Attach Dimension scale to the main dataset
  ///////
  /* Attach first dimension to the dataset */
  setname =  "/"+name+".real";
  realdataset_id = H5Dopen2(file_id,setname.c_str(), H5P_DEFAULT);
  setname =  "/"+name+".imag";
  imagdataset_id = H5Dopen2(file_id,setname.c_str(), H5P_DEFAULT);
  setname = "/X_AXIS";
  ds_id = H5Dopen2(file_id,setname.c_str(), H5P_DEFAULT);
  H5DSattach_scale(realdataset_id,ds_id,0);
  H5DSattach_scale(imagdataset_id,ds_id,0);
  /* Close ds_id after use it */
  status = H5Dclose(ds_id);
  
  /* Attach second dimension to the dataset */
  setname = "/Y_AXIS";
  ds_id = H5Dopen2(file_id,setname.c_str(), H5P_DEFAULT);
  H5DSattach_scale(realdataset_id,ds_id,1);
  H5DSattach_scale(imagdataset_id,ds_id,1);
  
  /* Close ds_id after use it, and also dataset_di */
  status = H5Dclose(ds_id);
  status = H5Dclose(realdataset_id);
  status = H5Dclose(imagdataset_id);
  
  /* Terminate access to the file. */
  status = H5Fclose(file_id); 

}

//////////////////////////////////////////////////////////
/*!
  @function save2DMatrix_txt
  @brief Stores 2D real matrices in a text file
  @params[in] std::string name 
  @params[in] arma::cx_mat matrix Real double data matrix
*/
//////////////////////////////////////////////////////////
void save2DMatrix_txt(const string name, const arma::mat mymatrix)
{ 
  ofstream myfile;
  int dims[2];
  string filename = name + ".dat";
  myfile.open(filename);
  
  dims[0] = mymatrix.n_cols; 
  dims[1] = mymatrix.n_rows;
  
  for (int iy=0; iy<dims[1]; iy++)
    for (int ix=0; ix<dims[0]; ix++)
      myfile << mymatrix(ix,iy) << endl;
  
  myfile.close();
}

//////////////////////////////////////////////////////////
/*!
  @function save2DMatrix_txt
  @brief Stores 2D complex matrices in a text file
  @params[in] std::string name 
  @params[in] arma::cx_mat matrix Real double data matrix
*/
//////////////////////////////////////////////////////////
void save2DMatrix_txt(const string name, const arma::cx_mat mymatrix)
{ 
  ofstream myfile;
  int dims[2];
  string filename = name + ".dat";
  myfile.open(filename);
  
  dims[0] = mymatrix.n_cols; 
  dims[1] = mymatrix.n_rows;
  
  for (int iy=0; iy<dims[1]; iy++)
    for (int ix=0; ix<dims[0]; ix++)
      myfile << real(mymatrix(ix,iy)) << " "
	     << imag(mymatrix(ix,iy)) << endl;
  
  myfile.close();
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

// For complex observables
void save_observable(const arma::vec& obs, const string& name)
{
  ofstream file;
  
  string filename = name + ".dat";
  file.open(filename);
  
  for(int i=0; i<obs.size(); ++i)
    file << obs(i) << endl;
  
  file.close();
  
}
//////////////////////////////////
// For complex observables
void save_observable(const arma::cx_vec& obs, const string& name)
{
  ofstream file;
  
  string filename = name + ".dat";
  file.open(filename);
  
  for(int i=0; i<obs.size(); ++i)
    file << real(obs(i)) << " " << imag(obs(i)) << endl;
  
  file.close();
  
}

// For complex observables
void save_observable(const arma::vec& time, const arma::vec& obs,
			   const string& name)
{
  ofstream file;
  
  string filename = name + ".dat";
  file.open(filename);
  
  for(int i=0; i<obs.size(); ++i)
    file << time(i) << " " << obs(i) << endl;
  
  file.close();
  
}
//////////////////////////////////

void save_observable(const arma::vec& time, const arma::cx_vec& obs,
			    const string& name)
{
  ofstream file;
  
  string filename = name + ".dat";
  file.open(filename);
  
  for(int i=0; i<obs.size(); ++i)
    file << time(i) << " " << real(obs(i)) << " "
	 << imag(obs(i)) << endl;
  
  file.close();
  
}


//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
