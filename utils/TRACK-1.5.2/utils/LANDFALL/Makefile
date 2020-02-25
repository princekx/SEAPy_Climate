SOURCES=read_tracks.c \
        new_time.c \
	realloc_n.c \
	convert_track_add.c \
	orog_test.c



OBJECTS=$(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm


.KEEP_STATE:

all: landfall_eps landfall


landfall_eps: $(OBJECTS) landfall_eps.c
	$(CC) $(CFLAGS) -o ../bin/$@ landfall_eps.c $(OBJECTS) $(LIBS)


landfall: $(OBJECTS) landfall.c
	$(CC) $(CFLAGS) -o ../bin/$@ landfall.c $(OBJECTS) $(LIBS)
