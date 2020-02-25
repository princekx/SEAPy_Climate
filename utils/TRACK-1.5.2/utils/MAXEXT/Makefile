SOURCES=maxext.c \
       meantrd.c \
       read_tracks.c \
       realloc_n.c

OBJECTS = $(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm

.KEEP_STATE:

all: maxext 

maxext: $(OBJECTS) 
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECTS) $(LIBS)

