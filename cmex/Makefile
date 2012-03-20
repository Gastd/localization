TARGET = main

CC = gcc
CFLAGS = -Wall -O2 -g -DUNIX -DLINUX -D_REENTRANT
LDFLAGS =  -lm
INCLUDES =
LIBDIR = 
OBJS =  ${TARGET}.o
SRCS =  ${TARGET}.c lib/gmatrix_linalg.c lib/gmatrix_sparse.c lib/gmatrix_statistics.c lib/gmatrix.c lib/gps.c lib/imu.c lib/kalman.c lib/localization.c lib/magnetometer.c lib/rotation.c lib/sonar.c

# Rules:
all: clean ${TARGET}

# The variable $@ has the value of the target.
${TARGET}:
#	${CC} -I${INCLUDES} ${CFLAGS} ${SRCS} -o $@ ${LDFLAGS} -L${LIBDIR}
	${CC} -I${INCLUDES} ${CFLAGS} ${SRCS} -o $@ ${LDFLAGS}
clean:
	rm ${TARGET} -f

.PHONY: clean
