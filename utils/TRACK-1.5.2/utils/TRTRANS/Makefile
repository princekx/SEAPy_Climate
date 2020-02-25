SOURCES=trtrans.c \
       trorigin.c \
       trarith.c \
       meantrd.c \
       read_tracks.c

OBJECTS = $(SOURCES:.c=.o)

#CFLAGS = -g -Xa -xsb

CFLAGS = -g -ansi -Wall -I../../include -D${SYSTEM} $(AD_C_FLAGS)

LIBS = -lm

.KEEP_STATE:

all: trtrans trorigin trarith tr2csv

trtrans: trtrans.o meantrd.o read_tracks.o
	$(CC) $(CFLAGS) -o ../bin/$@ trtrans.o meantrd.o read_tracks.o $(LIBS)
	
trorigin: trorigin.o meantrd.o read_tracks.o
	$(CC) $(CFLAGS) -o ../bin/$@ trorigin.o meantrd.o read_tracks.o $(LIBS)
	
trarith: trarith.o meantrd.o read_tracks.o
	$(CC) $(CFLAGS) -o ../bin/$@ trarith.o meantrd.o read_tracks.o $(LIBS)
	

tr2csv: tr2csv.o read_tracks.o
	$(CC) $(CFLAGS) -o ../bin/$@ tr2csv.o read_tracks.o $(LIBS)


