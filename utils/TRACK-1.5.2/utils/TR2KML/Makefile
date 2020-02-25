SOURCES= eps2kml.c \
         graceconv.c \
	 tr2kml.c \
         read_tracks.c \
         new_time.c \
         realloc_n.c \
         malloc_initl.c

OBJECTS = $(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm

.KEEP_STATE:

all: eps2kml tr2kml graceconv

eps2kml: eps2kml.o read_tracks.o new_time.o realloc_n.o malloc_initl.o
	$(CC) $(CFLAGS) -o ../bin/$@ eps2kml.o read_tracks.o new_time.o realloc_n.o malloc_initl.o $(LIBS)
	
tr2kml: tr2kml.o read_tracks.o
	$(CC) $(CFLAGS) -o ../bin/$@ tr2kml.o read_tracks.o $(LIBS)

graceconv: graceconv.o read_tracks.o
	$(CC) $(CFLAGS) -o ../bin/$@ graceconv.o read_tracks.o $(LIBS)




