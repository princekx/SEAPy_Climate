#
# Makefile for smoopy and sphery
#
      CC = gcc

  RANLIB = ranlib

  INCGCC = $(GCCINC)
 ADDINCS = 
AD_C_FLAGS =

  ARFLAGS =

#Include for older gcc compilation
#INCLUDES = -I../include -I../../include -I$(GCCINC) -I/usr/include -nostdinc $(ADDINCS)
#
INCLUDES = -I../include -I../../include -I$(INCGCC) -I/usr/include $(ADDINCS)

#INCLUDES = -I../include -I../../include $(ADDINCS)

 LIBRARY = ../libsmp.a

#
# CFLAGS for gcc, add -Wno-unknown-pragmas for egcs or gcc-2.95
CFLAGS = -ansi -g -D$(MDEBUG) -D$(SYSTEM) -D$(USAPPEND) -D$(G77) $(INCLUDES) -Wall -Wpointer-arith -Wshadow -Wnested-externs $(AD_C_FLAGS)

#
# CFLAGS for SUN cc (ANSI)
# code does not work with SunC due to using adjustable arrays which are 
# a GNU C extension. These must be replaced with dynamic routines to work
# with SunC.
#CFLAGS = -c -g -D$(MDEBUG) -D$(SYSTEM) -D$(USAPPEND) $(INCLUDES) -I/usr/include -Xa $(AD_C_FLAGS) #-xc99 #-xO4


  FFLAGS = -c $(N_FFLAGS)


C_SOURCES = bsplin_c_.c \
            error_interp.c \
            grismo_c_.c \
            grisph_c_.c \
            smoopy_c.c \
            smoopy_setup.c \
            sphery_c.c \
            sphery_setup.c \
            sqgrid_c_.c


F_SOURCES = bisp.f \
            bsplin.f \
            bsplin_all.f \
            cossin.f \
            cytri1.f \
            cytri2.f \
            disco.f \
            grismo.f \
            grisph.f \
            nknot.f \
            optsph.f \
            ration.f \
            rotate.f \
            smoopy.f \
            sphery.f \
            sqgrid.f \
            symsys.f

C_OBJECTS = $(C_SOURCES:.c=.o)
F_OBJECTS = $(F_SOURCES:.f=.o)

all: $(LIBRARY)

.KEEP_STATE:

C_OBJS: $(C_OBJECTS)

F_OBJS: $(F_OBJECTS)

$(LIBRARY): C_OBJS F_OBJS FORCE
	rm -f $@
	ar $(ARFLAGS) rc $@  $(C_OBJECTS) $(F_OBJECTS)
	$(RANLIB) $@
#	rm -f $(C_OBJECTS) $(F_OBJECTS)

depend:
	makedepend -- $(CFLAGS) -- $(C_SOURCES)
	../../sfmakedepend $(F_SOURCES)

FORCE:


