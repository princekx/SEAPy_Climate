#
# Makefile for libs
#
    MAKE = make
  RANLIB = ranlib

  SOURCE = src
INCLUDES = include

SRC_MAKEARGS = CC='$(CC)' G77=$(G77) MDEBUG=$(MDEBUG) RANLIB=$(RANLIB) SYSTEM=$(SYSTEM) ADDINCS=$(ADDINCS) USAPPEND=$(USAPPEND) N_FFLAGS='$(N_FFLAGS)'

all lopt: lib libchi libfft
surfit: libsurf libchi

.KEEP_STATE:

includes: 
	$(MAKE) $(INCLUDES)

$(INCLUDES): FORCE
	cd $@; $(MAKE); rm -f *%

lib:
	$(MAKE) $(SOURCE)

libsurf: includes FORCE
	cd $(SOURCE); $(MAKE) $(SRC_MAKEARGS)

$(SOURCE): includes FORCE
	cd $@; $(MAKE) $(SRC_MAKEARGS); \
	$(MAKE) -f Makefile.opt $(SRC_MAKEARGS); \
	rm -f *%

libchi: FORCE
	cd chisqr.src; $(MAKE) RANLIB=$(RANLIB) N_FFLAGS='$(N_FFLAGS)'
	rm -f *%

libfft: FORCE
	cd FFT; $(MAKE) RANLIB=$(RANLIB) N_FFLAGS='$(N_FFLAGS)'
	rm -f *%

depend: depend_src depend_chi depend_fft

depend_src:
	cd $(SOURCE); $(MAKE) $(SRC_MAKEARGS) depend; \
	$(MAKE) -f Makefile.opt $(SRC_MAKEARGS) depend

depend_chi:
	cd chisqr.src; $(MAKE) RANLIB=$(RANLIB) N_FFLAGS='$(N_FFLAGS)' depend

depend_fft:
	cd FFT; $(MAKE) RANLIB=$(RANLIB) N_FFLAGS='$(N_FFLAGS)' depend

FORCE:
