
/* CARE IF EDITING, editing this file (files_in.in) may 
   cause errors when the configure shell script is run.    */

/* Old version of gcc (pre 3.0) don't like spaces in Add */

#include "file_path.h"

/* data input file name */

#define  DATIN   Add(USER,PATHI,linux)

/* input file for time averaged data */

#define  DTIMAVG Add(USER,PATHI,TIME_AVG)

/* input file for frame times */

#define  FRTIMES Add(USER,PATHI,times.linux)

#define  DATCM   Add(USER,DATPATH,CMAP.dat)  /* country map data file name */
#define  DIMCM   3000               /* dimension of any country map data */


/* input zonal dmax */

#define  DATZN   Add(USER,DATPATH,zone.dat)

/* input adaptive phimax */


#define  DATAD   Add(USER,DATPATH,adapt.dat)
