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

### 
INCLUDE = ${MONTEC_INCLUDE} -I ${BOOST_INCLUDE} # ${ARMA_INCLUDE} ${OPENBLAS_INCLUDE}
LIBS = ${OPENBLAS_LIBS} -lopenblas ${SLEPC_SYS_LIB} -L${BOOST_LIB}

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
