.KEEP_STATE:

WORKING_JAVA := $(shell which javac)

all: symlinks bfilt combine convolve count diffstat ensemble hart java ibtref landfall lifecycle mident sample split tc tele_weights tilt transadd trtrans tr2kml tr2nc tr2pgsl runs

symlinks: FORCE
	test -d bin || mkdir bin
	test -d CLASS || mkdir CLASS
	chmod u+x makelinks.csh
	makelinks.csh

combine: FORCE
	cd COMBINE; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'

convolve: FORCE
	cd CONVOLVE; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'

count: FORCE
	cd COUNT; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)' 

diffstat: FORCE
	cd DIFFSTAT; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)' HAVNETC='NETCDF' ADDINCS='${ADDINCS}' ADLIBS='${ADLIBS}'

ensemble: FORCE
	cd ENSEMBLE; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
hart: FORCE
	cd HART; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
ibtref: FORCE
	cd IBT-REF; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
landfall: FORCE
	cd LANDFALL; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'

lifecycle: FORCE
	cd LIFECYCLE; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
mident: FORCE
	cd MULTI-IDENT; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
bfilt: FORCE
	cd MISC; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'

sample: FORCE
	cd SAMPLE; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)' HAVNETC='NETCDF' ADDINCS='${ADDINCS}' ADLIBS='${ADLIBS}'

split: FORCE
	cd SPLIT; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)' 

tc: FORCE
	cd TC; make SYSTEM='${SYSTEM}' HAVNETC='NETCDF' ADDINCS='${ADDINCS}' ADLIBS='${ADLIBS}' AD_C_FLAGS='$(AD_C_FLAGS)'

tele_weights: FORCE
	cd TELE_WEIGHTS; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'

tilt: FORCE
	cd TILT; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'

transadd: FORCE
	cd TRANSADD; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'

trtrans: FORCE
	cd TRTRANS; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
tr2kml: FORCE
	cd TR2KML; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
tr2nc: FORCE
	cd TR2NC; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)' HAVNETC='NETCDF' ADDINCS='${ADDINCS}' ADLIBS='${ADLIBS}'
	
tr2pgsl: FORCE
	cd TR2PGSL; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
runs: FORCE
	cd SIEGEL; make SYSTEM='${SYSTEM}' AD_C_FLAGS='$(AD_C_FLAGS)'
	
java: FORCE
ifdef WORKING_JAVA
	cd JAVA; make
else
	@echo "****WARNING****, no working java, java utils not compiled"
endif

FORCE:
