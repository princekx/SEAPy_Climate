SOURCES= eps2pgsl.c \
         read_tracks.c \
         new_time.c

OBJECTS = $(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm

.KEEP_STATE:

all: eps2pgsl 

eps2pgsl: eps2pgsl.o read_tracks.o new_time.o
	$(CC) $(CFLAGS) -o ../bin/$@ eps2pgsl.o read_tracks.o new_time.o $(LIBS)




