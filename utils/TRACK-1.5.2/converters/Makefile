# Makefile for bin2drs (binary to drs) and bin2nc (binary to netcdf)
#

CC = gcc

SOURCE_BIN2DRS = bin2drs.c
SOURCE_BIN2NC = bin2nc.c new_time.c

INCLUDES_DRS = -I/users/kih/local/include
INCLUDES_NC = -I/users/kih/local/NETCDF-3.6.1/include -I../include
LIBS_DRS = -L/users/kih/local/lib -ldrs /nerc/packages/forte/6.0/SUNWspro/WS6U2/lib/libF77.so.2 -lsunmath -lm
LIBS_NC  = -L/users/kih/local/NETCDF-3.6.1/lib -lnetcdf

.KEEP_STATE:


all: bin2drs bin2nc


bin2drs: ${SOURCE_BIN2DRS}
	${CC} ${INCLUDES_DRS} -o $@ ${SOURCE_BIN2DRS} ${LIBS_DRS}

bin2nc: ${SOURCE_BIN2NC}
	if [ ! -h new_time.c ]; then ln -s ../src/new_time.c new_time.c; fi
#	${CC} ${INCLUDES_NC} -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -o $@ ${SOURCE_BIN2NC} ${LIBS_NC}
	${CC} ${INCLUDES_NC} -o $@ ${SOURCE_BIN2NC} ${LIBS_NC}

