#####################################################
#####################################################
##              Makefile for SFApp                 ##
#####################################################
#####################################################
# Include settings for specific architecture

include sys/ballard_gnu.make

TARGET = ${BIN}/spectral_scan.exec
VPATH  = ${INCLUDE}

.PHONY: clean

default: ${TARGET}

# Include object files
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
include classes.lst

#####################################################
##             Module dependencies                 ##
#####################################################

params.o: 
communicator.o:
fiber.o: communicator.o
light.o: communicator.o fiber.o chkopts
spectralscan.o: params.o communicator.o fiber.o light.o 

#####################################################
##             		Rules                      ##
#####################################################


${TARGET}: ${OBJS}
	@echo 'Creating executable...'
	${PCXX} ${CXXFLAGS} -o ${TARGET} ${OBJS} \
	-I ${INCLUDE} -L ${LIBS}

clean::
	rm -f *.o *.optrpt *.gcno
	rm -f ${TARGET}
