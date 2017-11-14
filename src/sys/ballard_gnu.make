# -----------------------------------------------------------------------------#
#                                                                              #
# ballard gnu                                                                  #
#                							       #
# Set for intel i7 architecture with intel processors			       #
#                                                                              #
# -----------------------------------------------------------------------------#

# Paths to folder tree
HELM_ROOT    = ..
HELM_SRC     = ${HELM_ROOT}/src
HELM_INCLUDE = ${HELM_ROOT}/include
HELM_LIB     = ${HELM_ROOT}/lib
BIN          = ${HELM_ROOT}/bin

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
CXXFLAGS = ${DEB} ${OPT} -std=c++11 ${PETSC_CCPPFLAGS} ${SLEPC_CCPPFLAGS}

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
PETSC_DIR = /opt/petsc/3.7.7/mpich-gnu/complex

### SLEPc
SLEPC_DIR = /opt/slepc/3.7.4/mpich-gnu/complex

### 
INCLUDE = ${HELM_INCLUDE} # ${ARMA_INCLUDE} ${OPENBLAS_INCLUDE}
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
