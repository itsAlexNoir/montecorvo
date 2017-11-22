# -----------------------------------------------------------------------------#
#                                                                              #
# zahir gnu                                                                    #
#                							       #
# Set for linux machines with intel i7 architecture with intel processors      #
#                                                                              #
# -----------------------------------------------------------------------------#

# Paths to folder tree
HELM_ROOT    = ..
HELM_SRC     = ${HELM_ROOT}/src
HELM_INCLUDE = ${HELM_ROOT}/include
HELM_LIB     = ${HELM_ROOT}/lib
BIN          = ${HELM_ROOT}/bin

# Compiler
CC	 = gcc
CXX	 = g++
PCC 	 = mpicc
PCXX	 = mpicxx
LIBTOOL  = glibtool
compiler = gnu

# Debugging flags
#DEB         = -check all
#DEB         = -fprofile-arcs -ftest-coverage -fPIC -O0
#DEB          = -p -ffast-math
#DEB         =  -gen-interfaces -warn interfaces

# Optimisation flags
OPT         = -O3 -funroll-loops -ftree-vectorize -ffast-math -fcx-limited-range -fopenmp

# Vectorisation flags
#VEC = -funroll-loops -ftree-vectorize

CCFLAGS  = ${DEB} ${OPT} ${PETSC_CCFLAGS} ${SLEPC_CCFLAGS}
CXXFLAGS = ${DEB} ${OPT} -std=c++11 ${PETSC_CCPPFLAGS} ${SLEPC_CCPPFLAGS}

######### Paths to external libraries
### Eigen
#EIGEN_ROOT    = 
#EIGEN_INCLUDE = /usr/include/eigen3/Eigen
#EIGEN_LIBS    = ${EIGEN_ROOT}/lib

### Armadillo
ARMA_ROOT    = 
ARMA_INCLUDE = /usr/include/armadillo_bits
ARMA_LIBS    = /usr/lib

### Openblas
OPENBLAS_ROOT    = 
OPENBLAS_INCLUDE = /usr/includev
OPENBLAS_LIBS    = /usr/lib

### PETSc
#PETSC_DIR = /opt/petsc/3.8.0/openmpi-gnu/debug/complex
PETSC_DIR = /opt/petsc/3.8.0/openmpi-gnu/optim/complex

### SLEPc
#SLEPC_DIR = /opt/slepc/3.8.0/openmpi-gnu/debug/complex
SLEPC_DIR = /opt/slepc/3.8.0/openmpi-gnu/optim/complex

### 
INCLUDE = ${HELM_INCLUDE} #${ARMA_INCLUDE} ${OPENBLAS_INCLUDE}
LIBS = ${OPENBLAS_LIBS} -lopenblas ${SLEPC_EPS_LIB}

# -----------------------------------------------------------------------------#
#                                                                              #
# How to produce objects from source files                                     #
#                                                                              #
# -----------------------------------------------------------------------------#

# Fortran90

%.o:%.c
	${PCC} ${CCFLAGS} -c $< -I ${INCLUDE}

# Fortran77

%.o:%.cpp chkopts
	${PCXX} ${CXXFLAGS} -c $< -I ${INCLUDE}

# -----------------------------------------------------------------------------#
#                                                                              #
# End of architecture Makefile                                                 #
#                                                                              #
# -----------------------------------------------------------------------------#
