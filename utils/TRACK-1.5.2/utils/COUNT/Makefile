SOURCES=count.c \
       box.c \
       poly.c \
       meantrd.c \
       read_tracks.c \
       new_time.c

OBJECTS = $(SOURCES:.c=.o)

#OBJECT1=count.o meantrd.o read_tracks.o new_time.o sincos.o
OBJECT1=count.o meantrd.o read_tracks.o new_time.o

OBJECT2=box.o meantrd.o read_tracks.o new_time.o

OBJECT3=poly.o meantrd.o read_tracks.o new_time.o

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm


.KEEP_STATE:

all: count box poly

count: ${OBJECTS}
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT1) $(LIBS)

box: ${OBJECTS}
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT2) $(LIBS)
	
poly: ${OBJECTS}
	$(CC) $(CFLAGS) -o ../bin/$@ $(OBJECT3) $(LIBS)




