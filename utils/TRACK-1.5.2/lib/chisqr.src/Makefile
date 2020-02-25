#
# Makefile for chisqr library
#

  RANLIB = ranlib



 LIBRARY = ../libchisqr.a


  FFLAGS = -c $(N_FFLAGS)


F_SOURCES = cdfchi.f \
            cumchi.f \
            cumgam.f \
            dinvr.f \
            dzror.f \
            erf.f \
            erfc1.f \
            exparg.f \
            gam1.f \
            gamma.f \
            gratio.f \
            ipmpar.f \
            rexp.f \
            rlog.f \
            spmpar.f

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

