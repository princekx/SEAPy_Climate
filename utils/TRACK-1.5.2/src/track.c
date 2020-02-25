#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "mem_er.h"

#include "grid.h"

#include "st_fo.h"                /* data structure for feature points
                                     correlated with objects.           */
#include "st_obj.h"               /* data structure for objects         */
#include "st_track.h"             /* data structure for tracks correlated
                                     with objects.                      */
#include "file_path.h"
#include "files_in.h"             /* input files.                       */
#include "files_out.h"            /* output files.                      */
#include "splice.h"               /* data structure for combined track data. */

#include "file_handle.h"

#include "proj.h"
#include "boundary.h"


#include "netcdf_info.h"


#define  PERH    0
#define  PERV    0

#define  PERIOD  360.0

#define  MEXT    50


/* Program track.c is a general purpose tracking code.
   Field data for a series of frames is interpreted for objects
   of interest, these are then defined by a set of points which
   are tracked through the frames using an adapted optimization
   algorithm.

   Tracks can be subsequentely interpreted statisticaly and statistics for
   several periods combined.
 
                                                                        */ 

struct frame_objs *threshold(FILE * , int , int );
int search_area_check(int , int , int);
struct tot_tr *splice(int * );
struct track_ind *mge_tracks(struct frame_objs * , int );
void statistic(struct tot_tr *, int );
void std_read_header(GRID * , int * , FILE * , int * , int * );
void utf_read_header(float ,GRID * , int * , FILE * ,int * , int * );
void back_to_sphere(struct frame_objs * , int * , int, PRFP );
void objectd(struct frame_objs *, FILE * ,int , int );
void run_means(void);
void tendency(FILE * , int ,  int );
int pp_header(GRID * ,FILE * , int * , int * , int , int );
int pp_read(FILE * , int );
int time_avg(FILE * , int , int , int , int );
void spectral_filter(FILE * , int , int , int );
double **lanczos_create(int * , int * , int *);
void time_filt(FILE * , double ** , int , int , int , int , int , int );
void spec_filt(FILE *, int , int , int , int );

void convert(FILE * , int , int , int );
int extract(FILE * );
int parse_com(int , char ** , char * , char * , int , int );
void additional_fields(struct tot_tr * , int , off_t , FILE * , int , float );
void spline_smooth(FILE * );
void wtm_combine(FILE * );
void compute_vorticity(FILE * , off_t , int );
void compute_vorticity_fd(FILE * , off_t , int );
void compute_gradient(FILE * , off_t , int );
void time_tele_avg(FILE * );

void fast_spectral_filter(FILE * , int , int , int );
void interp_ng(FILE * , off_t , int );
void remo_rotated_correct(struct frame_objs * , float , int * , int , int );

void limited_area_filter(FILE * , int , int , int );

int x1u, x2u, y1u, y2u;                 /* search area grid integers */
int frnum;                              /* array dimensions for each frame */ 
int perh=PERH, perv=PERV;               /* switchs for periodic boudarys */
int form=1;                             /* format of input data file */
int gof=0;                              /* grid offset */

int tl='n';                             /* switch for grid offset */
int tom='e';                            /* switch for distance measure */
int init=0;                             /* initialization switch */
int psw='n';
int nw_ln='n';                          /* newline character switch for UNIX binary data */
int iext=0;                             /* flag for a change of output file extension.   */

char *fext=NULL;                        /* alternative file extension for output files */

float xmn, ymn, xmx, ymx;               /* search area numbers */

float period=PERIOD;

FILE *finit=NULL;

GRID *gr=NULL, *gt=NULL;                /* grid data pointer*/
CNTRY *cm=NULL, *ct=NULL;               /* country map data pointer */

GRID *gr1=NULL, *gr2=NULL;              /* subsidary grid pointers */
CNTRY *cm1=NULL, *cm2=NULL;             /* subsidary country map pointers */

PROJ *pp=NULL;

int aniso='n';

float *abuf=NULL;
void *databuf=NULL;

extern int fruser;
extern int tf;
extern int shp;
extern int track_num;
extern int fptwrt;

extern int nf, nfld;
extern int *nfwpos;

extern char *chrfld;
extern int objwr;


int main(int argc, char **argv)
{

   int i, snum=0, frame_num=0;
   int j, k;
   int idt;
   int sn1, sn2;
   int reo='n', scmp='n';
   int anty, ty=0;                   /* switchs for track analysis routines */
   int trackn=0;                   /* total number of tracks for analysis */
   int tf1=0, tf2=0;
   int nmo=0;
   int frtmp=0;
   int fstart, frate, fend;
   int trnm1=0, trnm2=0;
   int iper='n';
   int nfn;
   int nfb=0, forder=0;
   int nn;
   int tsamp=0;
   int ch;
   int nchr=0;
   int sp_typ=0;
   int ivty=0;

   off_t pl1=0, pl2=0, plw=0;

   float ival=0.0;
   float cset=0.0;
   float cmp_off=0;
   
   float gtr=0.0;

   double **lweights=NULL;   /* pointers to lanczos weights */


   FILE *fdatin=NULL, *fcm=NULL;        /* pointers to open data files */
   FILE *fobjo=NULL;
   FILE *lwght=NULL;
   FILE *finit_temp=NULL;

   char filnam[MAXCHR];
   char ifil[MAXCHR];
   char filtmp[MAXCHR];
   char ctybuf[MAXCHR];
   char lwfiln[MAXCHR];
   char initf[MAXCHR];
   char comobj[MAXCHR];
   char strtmp[MAXCHR];

   struct frame_objs *fo=NULL, *fs=NULL;      /* pointer to feature points */
   struct frame_objs *f1=NULL, *f2=NULL;
   struct object *ob1=NULL, *ob2=NULL;
   struct point *p1=NULL, *p2=NULL;
   struct feature_pts *fp1=NULL, *fp2=NULL;
   struct tot_tr *all_tr=NULL, *altr=NULL;    /* pointer for combined track data */
   struct track_ind *tind=NULL, *tinds=NULL;  /* pointer to track points */
   struct boundary_pt *b1=NULL, *b2=NULL;

   CNTRY *ctemp=NULL;

   PROJ proj={NULL, NULL};

   mem_debug(2);

   pp = &proj;


   setvbuf(stdout, NULL, _IONBF, 0);

/* parse command line arguments */

   strcpy(filnam, DATIN);
   strcpy(filtmp, DATIN);

   if(argc > 1) {

       fext = (char *)malloc_initl(MAXCHR*sizeof(char));
       mem_er((fext == NULL) ? 0 : 1, MAXCHR*sizeof(char));
       iext = parse_com(argc, argv, filnam, fext, MAXCHR, MEXT);
       strncpy(filtmp, filnam, MAXCHR);

       if(strstr(Add(USER,,), fext)  || 
          strstr(Add(PATHO,,), fext) || 
          strstr(Add(PATHI,,), fext))  {
          printf("****ERROR****, file extension %s may conflict with IO paths\r\n"
                 "               choose a different extension.               \n\n", EXTENSION);
          exit(1);
       }

   }

   if(!strstr(filnam, "http://")){

      if(!fexist(filnam, "r")){

         printf("***WARNING***, file %s                              \r\n"
                "               does not exist, input a valid filname\r\n"
                "               possibly with full path.               \n\n", filnam);
         scanf("%s", filnam);

         if(!fexist(filnam, "r")) {

            printf("***ERROR***, file %s does not exist or does not have the \r\n"
                   "             correct permissions, aborting.         \n\n", filnam);
            exit(1);

         }

         strncpy(filtmp, filnam, MAXCHR);

      }

   }


/* check file extension does not conflict with rest of file IO path */


   if(!iext && (strstr(Add(USER,,), EXTENSION)  || 
                strstr(Add(PATHO,,), EXTENSION) || 
                strstr(Add(PATHI,,), EXTENSION))) {
      printf("****ERROR****, file extension %s may conflict with IO paths\r\n"
             "               choose a different extension.               \n\n", EXTENSION);
      exit(1);
   }


   strncpy(strtmp, filnam, MAXCHR);

   if(strstr(strtmp, "http://") && strstr(strtmp, "@")){
      nchr = (int)(strstr(strtmp, "@") - (strstr(strtmp, "http://") + 6));
      for(i=0; i < nchr; i++) *(strstr(strtmp, "http://") + i + 7) = '*';

      printf("***INFORMATION***, data file:-                              \r\n"
             "                   %s                                       \r\n"
             "                   exists and will be opened for data input.\n\n", strtmp);

   }

   else {

      printf("***INFORMATION***, data file:-                              \r\n"
             "                   %s                                       \r\n"
             "                   exists and will be opened for data input.\n\n", filnam);

   }


   gr = (GRID * )malloc_initl(sizeof(GRID));
   mem_er((gr == NULL) ? 0 : 1, sizeof(GRID));

   cm = (CNTRY * )malloc_initl(sizeof(CNTRY));
   mem_er((cm == NULL) ? 0 : 1, sizeof(CNTRY));

/* if a country map data file exists, read data */

   printf("do you want to load a country map, 'y' or 'n'\n");
   scanf("\n");
   if((scmp=getchar()) == 'y'){

     if(!(fcm = fopen(DATCM, "r"))){  /* open country map data file */

        printf("***Warning*** error opening country map data file, check file exists\n");
	scmp ='n';

     }

     else {

        fscanf(fcm, "%d", &(cm->dcm));

        cm->cmi = (int * )calloc(cm->dcm, sizeof(int));
        mem_er((cm->cmi == NULL) ? 0 : 1, cm->dcm * sizeof(int));

        cm->cmxg = (float * )calloc(cm->dcm, sizeof(float));
        mem_er((cm->cmxg == NULL) ? 0 : 1, cm->dcm * sizeof(float));

        cm->cmyg = (float * )calloc(cm->dcm, sizeof(float));
        mem_er((cm->cmyg == NULL) ? 0 : 1, cm->dcm * sizeof(float));

        printf("do you want an x-offset applied to the country data \r\n"
               "perhaps to align with grid, 'y' or 'n'.             \n\n");

        scanf("\n");
        if(getchar() == 'y'){

           printf("input x-offset for country map in degrees\n");
           scanf("%f", &cmp_off);


        }

        for(i=0; i < cm->dcm; i++){

            fgets(ctybuf, MAXCHR, fcm);

            sscanf(ctybuf, "%d %f %f", cm->cmi+i, cm->cmxg+i, cm->cmyg+i);

            *(cm->cmxg + i) += cmp_off;
            if(*(cm->cmxg+i) < 0.) *(cm->cmxg+i) += period;
            else if(*(cm->cmxg+i) > period) *(cm->cmxg+i) -= period;

        }


     }

     if(fcm && fclose(fcm) != 0){

        printf("***error*** closing country map data file\n");

     }

   }
   
   ct = cm;
    
/* end of reading country map data, file closed */

/* initialize, data format, grid translation, projection */

   strncpy(initf, FINIT, MAXCHR);
   if(iext) strcpy(strstr(initf, EXTENSION), fext);

   printf("do you want to use an existing initialization, '0' for no, '1' for yes\n");

   scanf("%d", &init);
   if(init){

      printf("Which initialization file is required\n");
      scanf("%s", ifil);
      finit = open_file(ifil, "r");
      finit_temp = open_file(initf, "w");
      while((ch = fgetc(finit)) != EOF) fputc(ch, finit_temp);
      close_file(finit_temp, initf);
      fseeko (finit, (off_t)0, FSTART);

   }

   else

      finit = open_file(initf, "w");


/* read in the grid data and store in the arrays xgrid and ygrid.
   Application dependent, may need editing for different data file
   format. */

    printf("what format is the data\r\n"
           "Input '0' for UNIX binary format       \r\n"
           "Input '1' for standard ASCII format.   \r\n"
           "Input '2' for UTF format.              \r\n"
           "Input '3' for Met. Office PP format.   \r\n"
           "Input '4' for NETCDF.                  \n\n");

    if(init) fscanf(finit, "%d", &form);
    else{
      scanf("%d", &form);
      fprintf(finit, "%d\n", form);
    }
    
    if(form != 4) {
       fdatin = open_file(filnam, "r");

       pl1 = ftello(fdatin);

    }

    if(form == 0 || form == 1) {
       if(form == 0){
          printf("Are frames seperated by a newline character, 'y' or 'n'?\n");
          if(init){
             fscanf(finit, "\n");
             nw_ln = getc(finit);
          }
          else{
             scanf("\n");
             nw_ln = getchar();
             fprintf(finit, "%c\n", nw_ln);
          }
       }
       else nw_ln = 'y';
       
       std_read_header(gr, &frnum, fdatin, &tl, &gof);
    }
    else if(form == 2){

       printf("what is the time interval in days e.g. 0.25 day == 6hr\n");
       if(init) fscanf(finit, "%f", &ival);
       else{
         scanf("%f", &ival);
         fprintf(finit, "%f\n", ival);
       }
       utf_read_header(ival, gr, &frnum, fdatin, &tl, &gof);

    }

    else if(form == 3){

       pp_header(gr, fdatin, &tl, &gof, 1, 1);
       fseeko (fdatin, (off_t)0, FSTART);

       frnum = 0;

/* scan file for number of frames */

       printf("****INFORMATION****, scanning file for number of frames.\r\n"
              "                     Please Wait ........               \n\n");


       while(pp_read(fdatin, 0)) ++frnum; 
       clearerr(fdatin);
       fseeko (fdatin, (off_t)0, FSTART);

       printf("Number of frames in file is %d, not neccesarily all the same field.\n\n", frnum);


    }


    else if(form == 4){

       fdatin = (FILE *)netcdf_info(gr, &frnum, filnam, &tl, &gof);
       if(!fdatin){
          printf("****ERROR****, no NETCDF info.\n\n");
          exit(1);
       }

    }


    else {printf("***error*** incorrect data format index on input\n"); exit(1);}

    gt = gr;

    if(form != 4) pl2 = ftello(fdatin);

    if(scmp == 'y' && cm->dcm){

      if(gof > 0){
         cset = *(gr->xgrid + gr->ix - 1);
         for(i=0; i < cm->dcm; i++) if(*(cm->cmxg+i) > cset) *(cm->cmxg+i) -= period;
      }
      else if(gof < 0){
         cset = *(gr->xgrid);
         for(i=0; i < cm->dcm; i++) if(*(cm->cmxg+i) < cset) *(cm->cmxg+i) += period;
      }      

    }

/* translate x-axis to start from the origen */
 
   if(tl == 'y'){

       cset = *(gr->xgrid);
       if(cm->dcm){
          for(i=0; i < cm->dcm; i++) *(cm->cmxg+i) -= cset;
       }
       for(i=0; i < gr->ix; i++) *(gr->xgrid+i) -= cset;
       
       if(fabs(cset) > 0.0) gtr = cset;

       printf("****WARNING*****, translating grid by %f\n\n", cset);

    }

    if(gr->igtyp == 1){

       printf("The current grid is global but not periodic, do you\r\n"
              "want to make the grid periodic, 'y' or 'n'.        \r\n"
              "This is required for tracking on a global domain,  \r\n"
              "but is not required for spectral filtering.        \n\n");

       if(init){
          fscanf(finit, "\n");
          iper = getc(finit);
       }
       else{
          scanf("\n");
          iper = getchar();
          fprintf(finit, "%c\n", iper);
       }

       if(iper == 'y'){

          gr->ix += 1;
/*          gr->igtyp = 0; */
          gr->igc = 1;

          cset = *(gr->xgrid);
	  if(cm->dcm){
             for(i=0; i < cm->dcm; i++) {
                *(cm->cmxg+i) -= cset;
                if(*(cm->cmxg+i) > period) *(cm->cmxg+i) -= period;
                else if(*(cm->cmxg+i) < 0.) *(cm->cmxg+i) += period;
             }
	  }
          for(i=0; i < gr->ix - 1; i++) *(gr->xgrid+i) -= cset;
          *(gr->xgrid + gr->ix - 1) = period;
	  
	  if(fabs(cset) > 0.0) gtr = cset;

          printf("****WARNING*****, translating grid by %f\n\n", cset);

          printf(" the NEW grid dimensions are %d * %d\n\n", gr->ix, gr->iy);
       }

    }

    else if(gr->igtyp == 0){
      iper = 'y';
      gr->igc = 1;
    }

/* what kind of distance measure is required */

    printf("what distance norm is required, input 'e' for euclidean,\r\n"
           "input 'g' for geodesic norm (if working on a sphere)\n");

    if(init){
       fscanf(finit, "\n");
       tom = getc(finit);
    }
    else{
       scanf("\n");
       tom = getchar();
       fprintf(finit, "%c\n", tom);
    }

/*  change projection  */

    printf("***WARNING***,if the data is defined on a lat-long grid,       \r\n"
           "              a Plate Caree projection is used by default. \r\n\n");


    printf("Is a different projection required, 'y' or 'n'\r\n\n"
           "Use this if data are already on a projection different from Plate Caree.\n\n");

    if(init){
       fscanf(finit, "\n");
       psw = getc(finit);
    }
    else{
       scanf("\n");
       psw = getchar();
       fprintf(finit, "%c\n", psw);
    }

    if(psw == 'y' ) proj =  proj_group(-1, -1);
    

/* define the search area */

/* initialise extent */

    x1u = 1;
    x2u = gr->ix;
    y1u = 1;
    y2u = gr->iy;

    if(!gr->prgr){

       printf(" the maximum possible search area is\n ");

       printf("X= %f to %f \n", *(gr->xgrid), *(gr->xgrid+gr->ix-1));
       printf(" Y= %f to %f \n", *(gr->ygrid), *(gr->ygrid+gr->iy-1));
       printf("\n\n");

       printf("define a search area in terms of the grid numbering\n");

       printf("X grid numbers =\n");

       if(init) fscanf(finit, "%d %d", &x1u, &x2u);
       else{
          scanf("%d %d", &x1u, &x2u);
          fprintf(finit, "%d\n%d\n", x1u, x2u);
       }

       sn1 = search_area_check(x1u, x2u, gr->ix);

       printf("X= %f to %f\n\n", *(gr->xgrid+x1u-1), *(gr->xgrid+x2u-1));
    
       printf("Y grid numbers =\n");

       if(init) fscanf(finit, "%d %d", &y1u, &y2u);
       else{
          scanf("%d %d", &y1u, &y2u);
          fprintf(finit, "%d\n%d\n", y1u, y2u);
       }

       sn2 = search_area_check(y1u, y2u, gr->iy);

       if(sn1 && sn2) snum = 1;
   
       printf("Y= %f to %f\n\n", *(gr->ygrid+y1u-1), *(gr->ygrid+y2u-1));


    }



    xmn = *(gr->xgrid + x1u-1);
    ymn = *(gr->ygrid + y1u-1);
    xmx = *(gr->xgrid + x2u-1);
    ymx = *(gr->ygrid + y2u-1);


/* close initialization file and reset */

    if(init){close_file(finit, ifil); init = 0;}

    else close_file(finit, initf);
    finit = NULL;

/* is a collection of object and track data files required to be combined
   and displayed.                                                        */

   printf("do you want to combine existing sets of track data and display them,    \r\n"
          "or display exsisting sets of combined tracks,                           \r\n"
          "and/or perform statistical analysis and display,                        \r\n"
          "or compute time series analyses, e.g. mean, varience,                   \r\n"
          "filtered varience, filtered fields, or                                  \r\n"
          "combine weighted statistics from time series analysis,                  \r\n"
          "or perform spatial spectral filtering of fields,                        \r\n"
          "or convert data to binary format,                                       \r\n"
          "or extract fields from a file at a chosen sampling,                     \r\n"
          "or compute vorticity from wind fields,                                  \r\n"
          "or interpolate to a new grid.                                           \r\n"
          "y or n.                                                                 \n\n");

    scanf("\n");

    if(getchar() == 'y'){ 

analy:

       printf("do you want to use any of the analysis routines\r\n"
              "Input  '0' to exit\r\n"
              "Input  '1' to combine track data sets, analyse and display.          \r\n"
              "Input  '2' to display only existing statistical analyses             \r\n"
              "Input  '3' to compute time average, varience and                     \r\n"
              "           filtered varience fields, note this can also              \r\n"
              "           be done as part of the feature identification.            \r\n"
              "Input  '4' to perform spectral filtering of fields by,               \r\n"
              "           spherical harmonic decomposition.                         \r\n"
              "Input  '5' to perform a time domain filtering of fields using        \r\n"
              "           a Lanczos filter.                                         \r\n"
              "Input  '6' to perform a spectral domain (time) filtering of fields.  \r\n"
              "Input  '7' to convert data to standard binary format with            \r\n"
              "           header.                                                   \r\n"
              "Input  '8' to extract fields from a file at chosen sampling,         \r\n"
              "           fill missing data holes by simple interpolation,          \r\n"
              "           apply mathematical transformation and convert to          \r\n"
              "           binary format.                                            \r\n"
              "Input  '9' to perform field smoothing using B-splines.               \r\n"
              "Input '10' to combine weighted time analysis fields from option '3'  \r\n"
              "           seperate files for each period.                           \r\n"
              "Input '11' to combine weighted time analysis fields from option '3'  \r\n"
              "           single file for all periods.                              \r\n"
              "Input '12' to compute vorticty, EKE or wind speed from wind fields.  \r\n"
              "Input '13' to compute gradient fields and thermal front parameter.   \r\n"
              "Input '14' to interpolate data to a new grid using B-splines.        \n\n");


       scanf("%d", &anty);

       if(anty == 0) goto tidy;

       else if(anty ==1 || anty == 2) {

          if(anty == 1) all_tr = splice(&trackn);

          printf("do you want to perform any analyses, or display     \r\n"
                 "existing results:                                   \r\n"
                 "Input '0' to exit                                   \r\n"
                 "Input '1' for calculation of addtional field values \r\n"
                 "Input '2' for propog. speeds, and statistics        \r\n"
                 "Input '3' for combining statistical distributions   \n\n");

          scanf("%d", &ty);

          if(ty == 0) goto tidy;

          else if(ty == 1) additional_fields(all_tr, trackn, pl2, fdatin, iper, gtr);

          else if (ty == 2) statistic(all_tr, trackn);

          else if (ty == 3) run_means();

          else {

              printf("****ERROR****, no action for option %d\n\n", ty);

          }


       }

       else if(anty == 3){

           printf("What start, rate and end frames are required, fs, fr, fe\n\n");
           scanf("%d %d %d", &fstart, &frate, &fend);

           time_avg(fdatin, fstart, frate, fend, 0);

           goto tidy;

       }

       else if(anty == 4){


           printf("What start, rate and end frames are required, fs, fr, fe\n\n");
           scanf("%d %d %d", &fstart, &frate, &fend);

           printf("Use least squares, '0', fast spectral transform, '1' or limited area dct, '2' \n\n");
           scanf("%d", &sp_typ);

           if(!sp_typ){
	      printf("****WARNING****, data must be global.\n\n");
              spectral_filter(fdatin, fstart, frate, fend);
	   }
           else if(sp_typ == 1){
	      printf("****WARNING****, data must be global.\n\n");
              fast_spectral_filter(fdatin, fstart, frate, fend);
	   }
	   else if(sp_typ == 2) {
	      printf("****WARNING****, data must be on a uniform grid for a limited area of the globe.\n\n");
	      limited_area_filter(fdatin, fstart, frate, fend);
	   }
	   else {
	      printf("****ERROR****, option %d incorrect for choice of spectral filter.\n\n", sp_typ);
	      exit(1);
	   }

           goto tidy;


       }

       else if(anty == 5){

           if(form != 4) plw = ftello(fdatin);

           printf("Do you want a new set of Lanczos weights, or read weights from a file,\r\n"
                  "input 'c' for create new weights and 'r' to read from file.           \n\n");
           scanf("\n");
           if(getchar() == 'r'){

              printf("What is the file to read Lanczos weights from?\n\n");
              scanf("%s", lwfiln);

              lwght = open_file(lwfiln, "r");

              fscanf(lwght, "%d %d %d", &nfb, &forder, &tsamp);

              lweights = (double **)calloc(nfb, sizeof(double *));
              mem_er((lweights == NULL) ? 0 : 1, nfb * sizeof(double *));

              for(i=0; i < nfb; i++){

                 *(lweights + i) = (double *)calloc(forder+1, sizeof(double));
                 mem_er((*(lweights + i) == NULL) ? 0 : 1, (forder+1) * sizeof(double));

              }

              for(i=0; i<= forder; i++) {

                  fscanf(lwght, "%*d ");
                  for(j=0; j < nfb; j++)fscanf(lwght, "%le", *(lweights + j) + i);

                  printf("%3d ", i);
                  for(j=0; j < nfb; j++)printf("% 10.5e ", *(*(lweights + j) + i));
                  printf("\n");

 
              }

              close_file(lwght, lwfiln);

           }

           else {

               printf("****INFORMATION****, creating a new set of weights for Lanczos filtering.\n\n");

               lweights = lanczos_create(&nfb, &forder, &tsamp);

           }
           

           printf("****INFORMATION****, the  time average is computed first to determine the number of frames.\n\n");

           printf("What start, rate and end frames are required, fs, fr, fe\n\n");
           scanf("%d %d %d", &fstart, &frate, &fend);

           nn = time_avg(fdatin, fstart, frate, fend, 1);

           if(nn <= forder){

              printf("****ERROR****, insufficient data for choosen Lanczos order\n\n");
              goto tidy;


           }

           if(form != 4) fseeko(fdatin, plw, FSTART);

           time_filt(fdatin, lweights, nfb, forder, fstart, frate, fend, nn);
         
           goto tidy;

       }

       else if(anty == 6){

           if(form != 4) plw = ftello(fdatin);

           printf("****INFORMATION****, the  time average is computed first to determine the number of frames.\n\n");

           printf("What start, rate and end frames are required, fs, fr, fe\n\n");
           scanf("%d %d %d", &fstart, &frate, &fend);

           nn = time_avg(fdatin, fstart, frate, fend, 1);

           if(form != 4) fseeko(fdatin, plw, FSTART);

           spec_filt(fdatin, fstart, frate, fend, nn);
         
           goto tidy;

       }

       else if(anty == 7) {

           printf("What start, rate and end frames are required, fs, fr, fe\n\n");
           scanf("%d %d %d", &fstart, &frate, &fend);

           convert(fdatin, fstart, frate, fend);

           goto tidy;

       }

       else if (anty == 8) {

           extract(fdatin);
           goto tidy;

       }

       else if(anty == 9){

            spline_smooth(fdatin);
            goto tidy;

       }

       else if(anty == 10) wtm_combine(fdatin);

       else if(anty == 11) time_tele_avg(fdatin);

       else if(anty == 12) {
          printf("If computing vorticity use B-splines or finite differences, '0' for B-splines or '1' for finite differences.\n\n");
	  scanf("%d", &ivty);
	  if(ivty < 0 || ivty > 1){
	     printf("****ERROR****, wrong value for option.\n\n");
	     exit(1);
	  }
	  if(!ivty) compute_vorticity(fdatin, pl2, iper);
	  else compute_vorticity_fd(fdatin, pl2, iper);
       }

       else if(anty == 13) compute_gradient(fdatin, pl2, iper);
       
       else if(anty == 14) interp_ng(fdatin, pl2, iper);

       else {

          printf("****ERROR****, no action for chosen option %d\n\n", anty);

          goto tidy;

       }

       if(all_tr != NULL){             

          for(i=0; i < trackn; i++) {
	      altr = all_tr + i;
	      for(j=0; j < altr->num; j++) {if(altr->trpt) free((altr->trpt + j)->add_fld);}
	      free(altr->trpt);
	      
	  }

          free(all_tr);

       }

       goto analy;

    }
   

/* is an existing set of object data required to be used */

    printf("do you want to use an existing set of object and feature point\r\n"
          "data, y or n\n");

    scanf("\n");
    reo = getchar();

/* is tendency field data required calculated from the current field data */

    if(reo == 'n'){

       printf("Do you want to compute the tendency of the current data, 'y' or 'n'?\n\n");
       scanf("\n");
       if(getchar() == 'y'){

          if(form == 4) {
 
             printf("****WARNING****, tendency calculation does not support netcdf at present\r\n"
                    "                 use convert option to convert data to standard binary  \r\n"
                    "                 format and use this data to compute tendencies.        \n\n");
             exit(1);

          }

          close_file(fdatin, filnam);
          strncpy(filnam, TENDENCY, MAXCHR);
          if(iext) strcpy(strstr(filnam, EXTENSION), fext);

          if(fexist(filnam, "r")){

             fdatin = open_file(filnam, "r");

             fscanf(fdatin, "%d %d\n", &frtmp, &idt);
          
             printf("****WARNING****, the tendency file:              \r\n"
                    "                 %s                              \r\n"
                    "                 exists, with %d frames computed \r\n"
                    "                 using a frame seperation of %d .\n\n", 
                    filnam, frtmp, idt);

             printf("Do you want to use this file or compute a new data file for tendencies, 'y' or 'n'?\n");
             scanf("\n");

             if(getchar() == 'n'){

               close_file(fdatin, filnam);

               printf("What is the start and end frames, fs,  fe\n");
               scanf("%d %d", &fstart, &fend);


               fdatin = open_file(filtmp, "r");
               fseeko(fdatin, pl2, FSTART);
               tendency(fdatin, fstart, fend);
               close_file(fdatin, filtmp);
               fdatin = open_file(filnam, "r");
               fscanf(fdatin, "%d %d\n", &frnum, &idt);


             }

             else frnum = frtmp;

             
           }

           else{

              printf("What is the start and end frames, fs,  fe\n");
              scanf("%d %d", &fstart, &fend);

              fdatin = open_file(filtmp, "r");
              fseeko(fdatin, pl2, FSTART);
              tendency(fdatin, fstart, fend);
              close_file(fdatin, filtmp);
              fdatin = open_file(filnam, "r");
              fscanf(fdatin, "%d %d\n", &frnum, &idt);


           }

           printf("****WARNING****, data input now being read from the file:\r\n"
                  "                 %s \n\n", filnam);

       }

    }

/* produce thresholding of the user defined frames */

    if(gr->prgr){

       printf("\n**** Performing pass 1 through data ****\n\n");

       gt = gr;
       gr = gr1;
       ctemp = cm;
       cm = cm1;

       x1u = 1;
       x2u = gr->ix;
       y1u = 1;
       y2u = gr->iy;

       xmn = *(gr->xgrid + x1u-1);
       ymn = *(gr->ygrid + y1u-1);
       xmx = *(gr->xgrid + x2u-1);
       ymx = *(gr->ygrid + y2u-1);

       snum = 0;

       snum = (search_area_check(1, gr->ix, gr->ix) && 
               search_area_check(1, gr->iy, gr->iy)    ) ? 1 : 0;

       fo = threshold(fdatin, reo, snum);
       
       frame_num = fruser;
       if(tom == 'e' && tf){tind = mge_tracks(fo, frame_num); trnm1 = track_num;}

       tf1 = tf;

       if(gr2){

          printf("\n**** Performing pass 2 through data ****\n\n");

          gr = gr2;
          cm = cm2;

          x1u = 1;
          x2u = gr->ix;
          y1u = 1;
          y2u = gr->iy;

          xmn = *(gr->xgrid + x1u-1);
          ymn = *(gr->ygrid + y1u-1);
          xmx = *(gr->xgrid + x2u-1);
          ymx = *(gr->ygrid + y2u-1);

          snum = 0;

          snum = (search_area_check(1, gr->ix, gr->ix) && 
                  search_area_check(1, gr->iy, gr->iy)    ) ? 1 : 0;

          if(form != 4) fseeko(fdatin, pl2, FSTART);

          fs = threshold(fdatin, reo, snum);

          if(fruser != frame_num){
             printf("****ERROR****, the number of frames from both projections\r\n"
                    "               must be identical, exiting.               \n\n");
             exit(1);

          }

          frame_num = fruser;
          if(tom == 'e' && tf){tinds = mge_tracks(fs, frame_num); trnm2 = track_num;}

          tf2 = tf;

          if(tf1 != tf2){

             printf("***WARNING***, different feature detection methods have been used\n\n");

             if(!(tf1*tf2)){
                printf("***ERROR***, incompatable feature detection methods have been used, aborting program\n\n");
                exit(1);
             }

          }

       }
     

    }

    else  {fo = threshold(fdatin, reo, snum); tf1 = tf;}

    frame_num = fruser;

    if(form != 4)close_file(fdatin, filnam);
    else netcdf_close((NETCDF_INFO *)fdatin);


/* if geodesic measure is to be used project onto sphere */


    if(tom == 'g' && (gr->prty || gr->prgr == 3) && tf){

       switch(gr->prgr){
         case 0:

            proj = map_project(&gr->prty, 0);
            back_to_sphere(fo, &tf, frame_num, proj.prj1);
            break;

         case 1:
	 case 2:
	 case 3:

            gr = gr1;
            x1u = 1;
            y1u = 1;
            back_to_sphere(fo, &tf1, frame_num, NULL);
            tf = tf1;

            if(fs){gr = gr2; back_to_sphere(fs, &tf2, frame_num, NULL);}

            gr = gt; 
            gr->prty = 0;
            gr->prgr = 0;
            cm = ctemp;

            x1u = gr->ox1u = (gr2) ? ((gr1->ox1u < gr2->ox1u) ? gr1->ox1u : gr2->ox1u) : gr1->ox1u;
            x2u = gr->ox2u = (gr2) ? ((gr1->ox2u > gr2->ox2u) ? gr1->ox2u : gr2->ox2u) : gr1->ox2u;
            y1u = gr->oy1u = (gr2) ? ((gr1->oy1u < gr2->oy1u) ? gr1->oy1u : gr2->oy1u) : gr1->oy1u;
            y2u = gr->oy2u = (gr2) ? ((gr1->oy2u > gr2->oy2u) ? gr1->oy2u : gr2->oy2u) : gr1->oy2u;

            if(gr1){free(gr1->xgrid); free(gr1->ygrid); free(gr1);}
            if(gr2){free(gr2->xgrid); free(gr2->ygrid); free(gr2);}
            if(cm1){free(cm1->cmxg); free(cm1->cmyg); free(cm1->cmi); free(cm1);}
            if(cm2){free(cm2->cmxg); free(cm2->cmyg); free(cm2->cmi); free(cm2);} 

            break;
         }

/* merge fo and fs */

         if(fs){

            objwr = 0;

            strncpy(comobj,COMOBJ, MAXCHR);
            if(iext) strcpy(strstr(comobj, EXTENSION), fext);
            fobjo = open_file(comobj, "w");

            if(aniso == 'y') fprintf(fobjo, "%d %d\n", tf, 1);
            else fprintf(fobjo, "%d %d\n", tf, 0);

            fprintf(fobjo, "PROJ_DETAILS\n");

            fprintf(fobjo, "%d %d\n", gr->prgr, gr->prty);
            if(gr->prgr) fprintf(fobjo, "%f %f\n", gr->alat, gr->alng);

            fprintf(fobjo, "REGION_DETAILS\n");
            fprintf(fobjo, "%d %d %d %d\n", gr->ox1u, gr->ox2u, gr->oy1u, gr->oy2u);
	    
	    fprintf(fobjo, "ADDITIONAL_FIELDS\n");
	    fprintf(fobjo, "%3d %3d &", nf, nfld);
	    if(nfwpos){
	      for(i=0; i < nf; i++)fprintf(fobjo, "%1d", *(nfwpos + i));
	    }
	    fprintf(fobjo, "\n");

            for(i=0; i < frame_num; i++){

               f1 = fo + i;
               f2 = fs + i;

               ob2 = f2->objs;

               nmo = (f1->obj_num + f2->obj_num);

               (fo+i)->objs = (struct object * )realloc_n((fo+i)->objs, nmo * sizeof(struct object));
               mem_er(((fo+i)->objs == NULL) ? 0 : 1, nmo*sizeof(struct object));

               f1 = fo + i;
               ob1 = f1->objs;

               ob1 += f1->obj_num;          

               f1->obj_num += f2->obj_num;
               f1->tot_f_f_num += f2->tot_f_f_num;

               for(j=0; j < f2->obj_num; j++){
                  
                   ob1->b_or_i = ob2->b_or_i;
                   ob1->lab = ob2->lab;
                   ob1->mem_pt_size = ob2->mem_pt_size;
                   ob1->point_num = ob2->point_num;
                   ob1->bound_num = 0;
                   ob1->bound = NULL;
                   ob1->fet = NULL;
                   ob1->ext = NULL;

                   if(ob2->ext) {
                      ob1->ext = (struct extent * )malloc_initl(sizeof(struct extent));
                      mem_er((ob1->ext == NULL) ? 0 : 1, sizeof(struct extent));
                      *(ob1->ext) = *(ob2->ext);
                      free(ob2->ext);
                   }

                   ob1->pt = (struct point * )calloc(ob2->point_num, sizeof(struct point));
                   mem_er((ob1->pt == NULL) ? 0 : 1, ob2->point_num * sizeof(struct point));

                   p1 = ob1->pt;
                   p2 = ob2->pt;

                   for(k=0; k < ob2->point_num; k++) *(p1 + k) = *(p2 + k); 

                   free(ob2->pt);

                   if(ob2->bound_num){

                      ob1->bound_num = ob2->bound_num;
                      ob1->bound = (struct boundary_pt *)calloc(ob2->bound_num, sizeof(struct boundary_pt));
                      mem_er((ob1->bound == NULL) ? 0 : 1, ob2->bound_num * sizeof(struct boundary_pt));

                      b1 = ob1->bound;
                      b2 = ob2->bound;

                      for(k=0; k < ob2->bound_num; k++) *(b1 + k) = *(b2 + k);

                   }

                   free(ob2->bound);

                   if(ob2->fet){

                      ob1->fet = (struct features * )malloc_initl(sizeof(struct features));
                      mem_er((ob1->fet == NULL) ? 0 : 1, sizeof(struct features));

                      ob1->fet->fpt = (struct feature_pts * )calloc(ob2->fet->feature_num, sizeof(struct feature_pts));
                      mem_er((ob1->fet->fpt == NULL) ? 0 : 1, ob2->fet->feature_num * sizeof(struct feature_pts));

                      fp1 = ob1->fet->fpt;
                      fp2 = ob2->fet->fpt;
                      ob1->fet->feature_num = ob2->fet->feature_num;

                      for(k=0; k < ob2->fet->feature_num; k++) *(fp1 + k) = *(fp2 + k);

                      free(ob2->fet->fpt);
                      free(ob2->fet);

                   }

                   ++ob1;
                   ++ob2;

               }

               free(f2->objs);

/* write to file without object/boundary data */

               objectd(f1, fobjo, i+1, tf);

          

            }

            close_file(fobjo, comobj);

            free(fs);

         }

         xmn = *(gr->xgrid + x1u-1);
         ymn = *(gr->ygrid + y1u-1);
         xmx = *(gr->xgrid + x2u-1);
         ymx = *(gr->ygrid + y2u-1);

    }
    

/* undertake the tracking */

   tf = tf1;
   
   if(tom == 'g'){
     printf("Do you want to rotate feature locations to geographical locations \r\n"
            "to account for data relative to a rotated pole, 'y' or 'n'        \n\n");
     scanf("\n");
     if(getchar() == 'y') remo_rotated_correct(fo, cset, &tf, tl, frame_num);
   }

   if(!tind && tf) {tind = mge_tracks(fo, frame_num); trnm1 = track_num;}

/* tidy up and exit gracefully */

tidy:

   for(i=0; i < nfb; i++)free(*(lweights + i));
   free(lweights);

   if(form == 3) pp_header(NULL, NULL, NULL, NULL, -1, 0);

   if(gt) gr = gt;

   if(gr){
     free(gr->xgrid);
     if(shp == 'n') --(gr->ygrid); 
     free(gr->ygrid);
     free(gr->wrk);
     free(gr->aleng);
     if(gr->dleng) free(gr->dleng);
     free(gr->sico);
     free(gr);
   }
   
   if(ct) cm = ct;

   if(cm){

     free(cm->cmi);
     free(cm->cmxg);
     free(cm->cmyg);
     free(cm);

   }

   if(gr1){


      free(gr1->xgrid); 
      free(gr1->ygrid);
      free(gr1->aleng);
      if(gr1->dleng) free(gr1->dleng);
      free(gr1->sico);

      free(gr1);

      if(cm1){
         free(cm1->cmi);
         free(cm1->cmxg);
         free(cm1->cmyg);
         free(cm1);

      }

   }

   if(gr2){

      free(gr2->xgrid); 
      free(gr2->ygrid);
      free(gr2->aleng);
      if(gr2->dleng) free(gr2->dleng);
      free(gr2->sico);

      free(gr2);

      if(cm2){

         free(cm2->cmi);
         free(cm2->cmxg);
         free(cm2->cmyg);
         free(cm2);

      }

   }


   if(fo){


     nfn = (!tf || fptwrt) ? 1 : frame_num;

     for(i=0; i < nfn; i++){

         f1 = fo + i;

         if((f1->objs)){

            for(j=0; j < f1->obj_num; j++){

                ob1 = f1->objs + j;
                free(ob1->pt);
                if(ob1->fet) {
                  if(nf){
                    for(k=0; k < ob1->fet->feature_num; k++){
                       fp1 = (ob1->fet->fpt) + k;
                       free(fp1->add_fld);
                    }
		  }		
		  free(ob1->fet->fpt);
		}
                free(ob1->fet);
                free(ob1->ext);
                free(ob1->bound);

            }

            free(f1->objs);

         }


     } 

     free(fo);


   }

   if(tind){

      for(i=0; i<trnm1; i++) free((tind + i)->tp);
      free(tind);

   }

   if(tinds){

      for(i=0; i<trnm2; i++) free((tinds + i)->tp);
      free(tinds);

   }

   if(chrfld) free(chrfld);
   if(abuf) free(abuf);
   if(databuf) free(databuf);
   free(fext);


   return 0;
}

/* fix for linker bug */

void MAIN_(){}
