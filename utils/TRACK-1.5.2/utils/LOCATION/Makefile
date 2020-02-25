SOURCES=location.c \
        read_tracks.c \
        meantrd.c



OBJECTS=$(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb -I../../include

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm


.KEEP_STATE:

all: location



tilt: $(OBJECTS)
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECTS) $(LIBS)

