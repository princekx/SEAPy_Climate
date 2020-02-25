SOURCES=create_weights.c \
	meantrd_tanh.c \
	read_tracks.c


OBJECTS=$(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb -I../../include

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm


.KEEP_STATE:

all: tele_weights



tele_weights: $(OBJECTS)
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECTS) $(LIBS)

