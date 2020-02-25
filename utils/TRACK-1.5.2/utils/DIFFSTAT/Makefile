SOURCES=malloc_initl.c \
       read_stats.c \
       statdmp.c

OBJECTS=$(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include ${ADDINCS} -D${SYSTEM} -D${HAVNETC} $(AD_C_FLAGS) 


LIBS = -lm ${ADLIBS}

.KEEP_STATE:

all: diff_stats persist statconv meanstd

diff_stats: $(OBJECTS)
	$(CC) $(CFLAGS) -o ../bin/$@ diff.c $(OBJECTS) $(LIBS)

persist: $(OBJECTS)
	$(CC) $(CFLAGS) -o ../bin/$@ persist.c $(OBJECTS) $(LIBS)

statconv: read_stats.o malloc_initl.o 
	$(CC) $(CFLAGS) -o ../bin/$@ statconv.c netcdf_write_stats.c $(OBJECTS) $(LIBS)

meanstd: $(OBJECTS)
	$(CC) $(CFLAGS) -o ../bin/$@ meanstd.c $(OBJECTS) $(LIBS)
