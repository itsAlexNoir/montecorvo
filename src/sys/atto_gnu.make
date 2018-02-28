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
BIN          = ${MONTEC_ROOT}/bin

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
EIGEN_ROOT    = 
EIGEN_INCLUDE = ${EIGEN_ROOT}/include
EIGEN_LIBS    = ${EIGEN_ROOT}/lib

### Armadillo
ARMA_ROOT    = /opt/armadillo/8.300.3/gcc/7.2.1/usr
ARMA_INCLUDE = ${ARMA_ROOT}/include/
ARMA_LIBS    = ${ARMA_ROOT}/lib

### Openblas
OPENBLAS_ROOT    = /opt/openblas/0.2.20/gnu-7
OPENBLAS_INCLUDE = ${OPENBLAS_ROOT}/include
OPENBLAS_LIBS    = ${OPENBLAS_ROOT}/lib

### PETSc
#PETSC_DIR = /opt/petsc/3.8.2/mpich-gnu-7/debug/complex
PETSC_DIR = /opt/petsc/3.8.3/fast/complex/openmpi/3.0.0/gcc/7.2.1/

### SLEPc
#SLEPC_DIR = /opt/slepc/3.8.1/mpich-gnu-7/debug/complex
SLEPC_DIR = /opt/slepc/3.8.2/fast/complex/openmpi/3.0.0/gcc/7.2.1/

### 
INCLUDE = ${MONTEC_INCLUDE} -I ${ARMA_INCLUDE} # ${OPENBLAS_INCLUDE}
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
