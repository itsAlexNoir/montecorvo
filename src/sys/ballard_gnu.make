# -----------------------------------------------------------------------------#
#                                                                              #
# ballard gnu                                                                  #
#                							       #
# Set for intel i7 architecture with intel processors			       #
#                                                                              #
# -----------------------------------------------------------------------------#

# Paths to folder tree
MONTEC_ROOT    = ..
MONTEC_SRC     = ${MONTEC_ROOT}/src
MONTEC_INCLUDE = ${MONTEC_ROOT}/include
MONTEC_LIB     = ${MONTEC_ROOT}/lib
BIN            = ${MONTEC_ROOT}/bin

# Compiler
CC	 = gcc-7
CXX	 = g++-7
PCC 	 = mpicc
PCXX	 = mpicxx
LIBTOOL  = glibtool
compiler = gnu

# Debugging flags
#DEB         = -g -fprofile-arcs -ftest-coverage -fPIC -O0

# Optimisation flags
OPT         = -O3 -funroll-loops -ftree-vectorize -ffast-math -fcx-limited-range -fopenmp

# Vectorisation flags
#VEC = -funroll-loops -ftree-vectorize

CCFLAGS  = ${DEB} ${OPT} ${PETSC_CCFLAGS} ${SLEPC_CCFLAGS}
CXXFLAGS = ${DEB} ${OPT} -std=c++17 ${PETSC_CCPPFLAGS} ${SLEPC_CCPPFLAGS}

######### Paths to external libraries
### Eigen
EIGEN_ROOT    = /usr/local/opt/eigen
EIGEN_INCLUDE = ${EIGEN_ROOT}/include
EIGEN_LIBS    = ${EIGEN_ROOT}/lib

### Armadillo
ARMA_ROOT    = /usr/local/opt/armadillo
ARMA_INCLUDE = ${ARMA_ROOT}/include
ARMA_LIBS    = ${ARMA_ROOT}/lib

### Openblas
OPENBLAS_ROOT    = /usr/local/opt/openblas
OPENBLAS_INCLUDE = ${OPENBLAS_ROOT}/include
OPENBLAS_LIBS    = ${OPENBLAS_ROOT}/lib

### Boost
BOOST_ROOT    = /usr/local/opt/boost
BOOST_INCLUDE = ${BOOST_ROOT}/include
BOOST_LIB     = ${BOOST_ROOT}/lib

### PETSc
#PETSC_DIR = /opt/petsc/3.8.3/debug/complex/mpich/3.2.1/gcc/7.3.0/
PETSC_DIR = /opt/petsc/3.8.3/fast/complex/mpich/3.2.1/gcc/7.3.0/

### SLEPc
#SLEPC_DIR = /opt/slepc/3.8.2/debug/complex/mpich/3.2.1/gcc/7.3.0/
SLEPC_DIR = /opt/slepc/3.8.1/fast/complex/mpich/3.2.1/gcc/7.3.0/

### ZLIB
ZLIB_DIR     = /opt/zlib/1.2.11/gcc/7.3.0
ZLIB_INCLUDE = ${ZLIB_DIR}/include
ZLIB_LIB     = ${ZLIB_DIR}/lib

### SZIP
SZIP_DIR     = /opt/szip/2.1.1/gcc/7.3.0/
SZIP_INCLUDE = ${SZIP_DIR}/include
SZIP_LIB     = ${SZIP_DIR}/lib

### HDF5
HDF5_DIR     = /opt/hdf5/1.10.1/fast/mpich/3.2.1/gcc/7.3.0/
HDF5_INCLUDE = ${HDF5_DIR}/include
HDF5_LIB     = ${HDF5_DIR}/lib

### 
INCLUDE = ${MONTEC_INCLUDE} -I${HDF5_INCLUDE} \
	-I${ZLIB_INCLUDE} -I${SZIP_INCLUDE} # ${ARMA_INCLUDE} ${OPENBLAS_INCLUDE}
LIBS = ${OPENBLAS_LIBS} -lopenblas ${SLEPC_SYS_LIB} -L${HDF5_LIB} \
	${HDF5_LIB}/libhdf5_hl.a ${HDF5_LIB}/libhdf5.a\
	-L${ZLIB_LIB} -L${SZIP_LIB} -lsz -lz -ldl -lm

# -----------------------------------------------------------------------------#
#                                                                              #
# How to produce objects from source files                                     #
#                                                                              #
# -----------------------------------------------------------------------------#

# Fortran90

%.o:%.c
	${PCC} ${CCFLAGS} -c $< -I ${INCLUDE}

# Fortran77

%.o:%.cpp
	${PCXX} ${CXXFLAGS} -c $< -I ${INCLUDE}

# -----------------------------------------------------------------------------#
#                                                                              #
# End of architecture Makefile                                                 #
#                                                                              #
# -----------------------------------------------------------------------------#
