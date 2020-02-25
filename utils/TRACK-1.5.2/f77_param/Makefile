# Makefile for FORTRAN include files


CPPFLAGS = -P -C -D$(MNAME) -nostdinc
CPP = cpp

FILES.ptr = paramf.ptr \
            paramc.ptr \
            params.ptr \
            params_lab.ptr \
            path.ptr \
            paramcol.ptr



.KEEP_STATE:

all: $(FILES.ptr)

paramf.ptr: paramf.cpp
	$(CPP) $(CPPFLAGS) paramf.cpp $@ 

paramc.ptr: paramc.cpp
	$(CPP) $(CPPFLAGS) paramc.cpp $@

params.ptr: params.cpp
	$(CPP) $(CPPFLAGS) params.cpp $@


params_lab.ptr: params_lab.cpp
	$(CPP) $(CPPFLAGS) params_lab.cpp $@

paramcol.ptr: paramcol.cpp
	$(CPP) $(CPPFLAGS) paramcol.cpp $@

