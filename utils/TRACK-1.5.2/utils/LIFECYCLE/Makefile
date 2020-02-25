SOURCES=meantrd.c \
       read_tracks.c \
       realloc_n.c

OBJECTS = $(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm

.KEEP_STATE:

all: lifecycle vinteg combine_add

lifecycle: read_tracks.o
	$(CC) $(CFLAGS) -c lifecycle.c
	$(CC) $(CFLAGS) -o ../bin/$@ lifecycle.o read_tracks.o $(LIBS)

vinteg: $(OBJECTS)
	$(CC) $(CFLAGS) -c vinteg.c
	$(CC) $(CFLAGS) -o ../bin/$@ vinteg.o $(OBJECTS) $(LIBS)

combine_add: $(OBJECTS)
	$(CC) $(CFLAGS) -c combine_add.c
	$(CC) $(CFLAGS) -o ../bin/$@ combine_add.o $(OBJECTS) $(LIBS)


