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
CC	 = nvcc
CXX	 = nvcc
PCC 	 = 
LIBTOOL  = 
compiler = cuda

# Debugging flags
#DEB         = -g

# Optimisation flags
#OPT         = -O3 -funroll-loops -ftree-vectorize -ffast-math -fcx-limited-range -fopenmp

# Vectorisation flags
#VEC = -funroll-loops -ftree-vectorize

ARCH = -arch=sm_20

CXXFLAGS = ${ARCH} -x cu -std=c++11 ${DEB}

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
OPENBLAS_INCLUDE = /usr/include
OPENBLAS_LIBS    = /usr/lib

### 
INCLUDE = ${SFAPP_INCLUDE} #${ARMA_INCLUDE} ${OPENBLAS_INCLUDE}
LIBS = ${OPENBLAS_LIBS} -lopenblas

# -----------------------------------------------------------------------------#
#                                                                              #
# How to produce objects from source files                                     #
#                                                                              #
# -----------------------------------------------------------------------------#

# C

%.o:%.c
	${CC} ${CCFLAGS} -c $< -I ${INCLUDE}

# C++

%.o: %.cpp
	${CC} ${CXXFLAGS} -I ${INCLUDE} -dc $< -o $@

# -----------------------------------------------------------------------------#
#                                                                              #
# End of architecture Makefile                                                 #
#                                                                              #
# -----------------------------------------------------------------------------#
