SOURCES=multi_ident.c \
       meantrd.c \
       read_tracks.c 

OBJECTS = $(SOURCES:.c=.o)

OBJECT=multi_ident.o meantrd.o read_tracks.o 

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = 


.KEEP_STATE:

all: mident 

mident: ${OBJECTS}
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT) $(LIBS)
