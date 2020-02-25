
/* CARE IF EDITING, editing this file (files_out.in) may 
   cause errors when the configure shell script is run.    */

/* Old version of gcc (pre 3.0) don't like spaces in Add */

#include "file_path.h"

#define  NOTRACK   0     /* set to 1 for track initialization but no tracking */
                         /* set to 0 otherwise.                               */

#define  MAXCHR  500

#define  EXTENSION  "linux"

/* initialization file */

#define  FINIT         Add(USER,PATHO,initial.linux)

/* output file for thresholded data                      */

#define  DATTHRO       Add(USER,PATHO,throut.linux)

/* output file for points which constitute objects       */

#define  DOUTOBJ       Add(USER,PATHO,objout.linux)

#define  COMOBJ        Add(USER,PATHO,comb_obj.linux)
     
/* object data output file after feature point filtering */

#define  NEWOBJF       Add(USER,PATHO,objout.new.linux)

/* output file for final track data                      */

#define  TDUMP         Add(USER,PATHO,tdump.linux)

/* output file for initial track data                    */
 
#define  IDUMP         Add(USER,PATHO,idump.linux)

/* time average file                                     */

#define  TAVGE         Add(USER,PATHO,user_tavg.linux)

/* tendency file                                         */

#define  TENDENCY      Add(USER,PATHO,field_tend.linux)

/* spectraly filtered fields output stub filename        */

#define  SPECTRAL      Add(USER,PATHO,specfil.linux)

/* spectral response for Lanczos filtering               */

#define  LANCZOS_RESP  Add(USER,PATHO,lanczos_resp.linux)

/* Lanczos weights                                       */

#define  LANCZOS_W     Add(USER,PATHO,lanczos_w.linux)

/* time domain filtered fields output stub filename      */

#define  TIME_FILT     Add(USER,PATHO,filt_time.linux)

/* data conversion output filename stub                  */

#define  CONVERT       Add(USER,PATHO,conv_bin.linux)

/* data extraction output filename stub                  */

#define  EXTRACT       Add(USER,PATHO,extract.linux)

/* data smoothing output filename stub                  */

#define  SMOOTH       Add(USER,PATHO,smooth.linux)

/* data masking output filename stub                    */

#define  MASK       Add(USER,PATHO,mask.linux)

/* data output from threshold on current projection     */

#define  INTERP_TH       Add(USER,PATHO,interp_th.linux)
