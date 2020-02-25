
/* CARE IF EDITING, editing this file (file_cat_out.in) may 
   cause errors when the configure shell script is run.    */

/* output files for spliced track data sets and mean tracks*/
/* Old version of gcc (pre 3.0) don't like spaces in Add */

#include "file_path.h"

#define  MAXCHR 500         /* number of characters in an input file name */

#define  EXTENSION  "linux"

/* spliced output file for track data            */

#define  FPTTRS      Add(USER,PATHO,tr_trs.linux)

/* filtered spliced track file                   */

#define  FILTRS      Add(USER,PATHO,ff_trs.linux)

/* grid point track data                         */

#define  GRTRS       Add(USER,PATHO,tr_grid.linux)

/* mean tracks output file                       */

#define  MNTRS       Add(USER,PATHO,mean_trs.linux)

/* phase speed output file                       */

#define  PHTRS       Add(USER,PATHO,phase_trs.linux)

/* statistics output file                        */

#define  STATTRS     Add(USER,PATHO,stat_trs.linux)

/* scaled statistics output file                */

#define  STATTRS_SCL Add(USER,PATHO,stat_trs_scl.linux)

/* combined statistics output file              */

#define  STATCOM     Add(USER,PATHO,stat_com.linux)

/* initiation output file                       */

#define  INITTRS     Add(USER,PATHO,init_trs.linux)

/* disapearance output file                     */

#define  DISPTRS     Add(USER,PATHO,disp_trs.linux)
