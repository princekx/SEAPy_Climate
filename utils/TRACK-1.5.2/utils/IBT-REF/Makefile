SOURCES=ibtref.c \
       meantrd.c \
       read_tracks.c \
       realloc_n.c

OBJECTS = $(SOURCES:.c=.o)

OBJECT=ibtref.o meantrd.o read_tracks.o realloc_n.o

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm


.KEEP_STATE:

all: ibtref

ibtref: ${OBJECTS}
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT) $(LIBS)





