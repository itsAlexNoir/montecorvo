#####################################################
#####################################################
##              Makefile for SFApp                 ##
#####################################################
#####################################################
# Include settings for specific architecture

include sys/ballard_gnu.make

TARGET = ${BIN}/helmholtz.exec
VPATH  = ${INCLUDE}

.PHONY: clean

default: ${TARGET}

# Include object files
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
include classes.lst

#####################################################
##             Module dependencies                 ##
#####################################################

communicator.o:
space.o: communicator.o
light.o: communicator.o space.o chkopts
helmholtz.o: space.o

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
