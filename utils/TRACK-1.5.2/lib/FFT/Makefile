#
# Makefile for FFT library
#

  RANLIB = ranlib



 LIBRARY = ../libfft.a


  FFLAGS = -c $(N_FFLAGS)


F_SOURCES = gpfa.f \
            gpfa2f.f \
            gpfa3f.f \
            gpfa5f.f \
            setgpfa.f


F_OBJECTS = $(F_SOURCES:.f=.o)

all: $(LIBRARY)

.KEEP_STATE:


F_OBJS: $(F_OBJECTS)

$(LIBRARY):F_OBJS FORCE
	rm -f $@
	ar $(ARFLAGS) rc $@ $(F_OBJECTS)
	$(RANLIB) $@
#	rm -f $(F_OBJECTS)

depend:
	../../sfmakedepend $(F_SOURCES)

FORCE:

