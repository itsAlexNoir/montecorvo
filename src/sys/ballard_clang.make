# -----------------------------------------------------------------------------#
#                                                                              #
# ballard gnu                                                                  #
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
CC	 = gcc
CXX	 = g++
PCC 	 = /usr/local/opt/openmpi/bin/mpicc
PCXX	 = /usr/local/opt/openmpi/bin/mpicxx
LIBTOOL  = glibtool
compiler = clang

# Debugging flags
#DEB         = -g -fprofile-arcs -ftest-coverage -fPIC -O0

# Optimisation flags
#OPT         = -O3 -funroll-loops -ftree-vectorize -ffast-math -fcx-limited-range -fopenmp

# Vectorisation flags
#VEC = -funroll-loops -ftree-vectorize

CCFLAGS  = ${DEB} ${OPT} ${PETSC_CCFLAGS}	
CXXFLAGS = ${DEB} ${OPT} -std=c++11 ${PETSC_CCPPFLAGS}

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

### PETSc
PETSC_DIR = /usr/local/opt/petsc/complex

### SLEPc
SLEPC_DIR = /usr/local/opt/slepc/complex

### 
INCLUDE = ${SFAPP_INCLUDE} # ${ARMA_INCLUDE} ${OPENBLAS_INCLUDE}
LIBS = ${OPENBLAS_LIBS} -lopenblas ${SLEPC_SYS_LIB}

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
