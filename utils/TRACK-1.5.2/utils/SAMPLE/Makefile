SOURCES=sample.c \
        meantrd.c \
        read_tracks.c \
        ran3.c \
        maxmin.c \
        read_stats.c \
        netcdf_write_stats.c \
        malloc_initl.c \
        confidence.c \
        statdmp.c \
        pval.c \
        missconv.c \
        z0.c \
        cnorminv.c \
        acc.c



OBJECTS=$(SOURCES:.c=.o)

OBJECT1=sample.o meantrd.o read_tracks.o ran3.o
OBJECT2=maxmin.o read_stats.o malloc_initl.o
OBJECT3=confidence.o read_stats.o malloc_initl.o statdmp.o ran3.o
OBJECT4=pval.o read_stats.o statdmp.o malloc_initl.o netcdf_write_stats.o
OBJECT5=z0.o read_stats.o statdmp.o cnorminv.o malloc_initl.o
OBJECT6=acc.o read_stats.o statdmp.o cnorminv.o malloc_initl.o

#CFLAGS = -g -Xa -xsb -I../../include

CFLAGS = -g -ansi -Wall -I../../include ${ADDINCS} -D${SYSTEM} -D${HAVNETC} $(AD_C_FLAGS) -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

LIBS = -lm ${ADLIBS}


.KEEP_STATE:

all: sample maxmin confidence pval missconv z0 acc



sample: $(OBJECT1)
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT1) $(LIBS)

maxmin: $(OBJECT2)
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT2) $(LIBS)

confidence: $(OBJECT3)
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT3) $(LIBS)

pval: $(OBJECT4)
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT4) $(LIBS)

missconv: missconv.o realloc_n.o
	$(CC) $(CFLAGS) -o ../bin/$@ missconv.o realloc_n.o $(LIBS)

z0: $(OBJECT5)
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT5) $(LIBS)

acc: $(OBJECT6)
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT6) $(LIBS)
	
