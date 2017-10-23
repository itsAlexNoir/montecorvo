# -----------------------------------------------------------------------------#
#                                                                              #
# ballard intel                                                                #
#                							       #
# Set for intel i7 architecture with intel processors			       #
#                                                                              #
# -----------------------------------------------------------------------------#

# Paths to folder tree
SFAPP_ROOT    = ..
SFAPP_SRC     = ${SFAPP_ROOT}/src
SFAPP_INCLUDE = ${SFAPP_ROOT}/include
SFAPP_LIB     = ${SFAPP_ROOT}/lib
BIN           = ${SFAPP_ROOT}/bin

# Compiler
CC	 = icc
CXX	 = icpc
PCC 	 = mpicc
LIBTOOL  = glibtool
compiler = intel

# Debugging flags
#DEB		= -no-vec
#DEB         = -check all
#DEB         = -g -traceback
#DEB         = -gen-interfaces -warn interfaces

# Optimisation flags
OPT		= -O2 -qopt-report-phase=vec -qopt-report=1
#OPT         = -O3

# Vectorisation flags
#VEC =	

CXXFLAGS = ${DEB} ${OPT} -std=c++11

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
OPENBLAS_LIBS    = ${OPENBLAS_ROOT}/lib -lopenblas

### MKL
MKL_ROOT = ${MKLROOT}
MKL_LIBS = ${MKLROOT}/lib -lmkl_intel_lp64 \
	 -lmkl_intel_thread -lmkl_core -qopenmp

### 
INCLUDE = ${SFAPP_INCLUDE} # ${ARMA_INCLUDE} ${OPENBLAS_INCLUDE}
LIBS = ${MKL_LIBS}

# -----------------------------------------------------------------------------#
#                                                                              #
# How to produce objects from source files                                     #
#                                                                              #
# -----------------------------------------------------------------------------#

# Fortran90

%.o:%.c
	$(CC) $(CCFLAGS) -c $< -I ${INCLUDE}

# Fortran77

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $< -I ${INCLUDE}

# -----------------------------------------------------------------------------#
#                                                                              #
# End of architecture Makefile                                                 #
#                                                                              #
# -----------------------------------------------------------------------------#
