#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include <sys/types.h>
#include "mem_er.h"
#include "st_im.h"
#include "st_obj.h"
#include "st_fo.h"
#include "files_in.h"
#include "files_out.h"
#include "file_handle.h"
#include "grid.h"
#include "reg_dat.h"
#include "sphery_dat.h"
#include "bisp.h"
#include "proj.h"
#include "boundary.h"
#include "threshold.h"
#include "netcdf_info.h"
#include "m_values.h"

#define  TYPMIN   0
#define  TYPMAX   1
#define  MAXCHAR 20
#define  DELB    6
#define  OBR    -0.0
#define  FBR    20.0
#define  DIMSM   3

#define  INTERP_DAT  1


/* function to perform pre-processing, theshold and find feature
   points for each user defined frame                             */


float *ap=NULL;

float obr=OBR, fbr=FBR;
float offs, ilat;

int tf, fruser;
int delb=DELB, dbp=DELB+1, dbm=DELB-1;
int xl, xr;
int cc='e', bs='n', hemi='0', of;
int xdm, ydm;
int pb='n';
int objwr=1, fptwrt;
int ffirst=1, fint=1, frterm=0;
int frcnt=0;
int i_shape=0;                          /* indicator for re-read object data   */
int irmiss=0;                           /* missing frame check flag            */


long int frtime=0;                      /* frame time information              */

extern int aniso;
extern int x1u, x2u, y1u, y2u;          /* search area grid numbers            */
extern int frnum;                       /* total numer of frames in data file. */
extern float xmn, ymn, xmx, ymx;        /* region extent for identification.   */
extern GRID *gr, *gt, *gr1, *gr2;       /* pointers to grid data.              */
extern CNTRY *cm, *cm1, *cm2;           /* pointers to country map data        */
extern char *chrfld;
extern int dfil;
extern char *fext;
extern int iext;
extern int form;
extern int tom;

extern int nf, nfld;
extern int *nfwpos;

struct frame_objs *threshold(FILE *fdatin, int reo, int sn)

{


     int i=0, j=0;
     int frcount, frold=0;
     int type, level, lmax;
     int multi, multis, filt, pflag;
     int dimp, xdim=x2u-x1u+1, ydim=y2u-y1u+1;
     int c='n', frd, as='n', tscl, tsth='\0';
     int xddim=0, yddim=0;
     int trans[3]={0, 0, 0};
     int trn=0;
     int *llb=NULL, *rrb=NULL;
     int pgr=0, pty=0;
     int hty=0;
     int smsw=0;
     int rf=0, rpl=0;
     int dumi=0;
     int ssm='n';
     int mfilt='n', mfiltd='n';
     int unp_filt='n';
     int itpfc=0;
     int adj_ob='n';
     int nmiss=0, tstep;
     int delobd='n';

     int rp_s, rp_i, rp_e;
     int bf='n';
     int b_exc=0;
     
     off_t chrnum=0, place1=0, place2=0;
     long int frtim1, frtim2;

     float *app=NULL, *apt=NULL;               /* pointers to current field data  */
     float *aplw=NULL, *applw=NULL;            /* pointers to planetary wave data */
     float thresh=0.;                          /* user set threshold              */
     float scale=1.0;                          /* field scaling                   */
     float sign=1.0, *ata=NULL, *ataa=NULL;
     float alat, alng;
     float tavoff=0.;
     float dumf;
     float upint1=0.0, upint2=0.0;
     float rads=0.0;

     FILE *fthro=NULL, *fobjo=NULL, *fta=NULL;
     FILE *fplw=NULL;
     FILE *finp=NULL;

     char charin[MAXCHR], *filpt=NULL;
     char fout1[MAXCHR],fout2[MAXCHR], foutn[MAXCHR];
     char tavin[MAXCHR], filplw[MAXCHR];
     char fintp[MAXCHR];

     struct image *ip=NULL, **hierarchy=NULL;
     struct frame_objs *fo=NULL, *ff=NULL, *fp=NULL, *fr=NULL, *ft=NULL;
     
     struct object *ob=NULL;
     struct boundary_cntl bcntl = {0, 0, 0, 'n', 0};

/* pointers to data structures for surface fitting. May need to change to 
   pointer to pointer later.                                               */

     int smty[DIMSM]={0, 0, 0};         /* spline type, smoopy or sphery        */
     struct sp_dat *cokn[DIMSM];        /* spline data for any of the methods   */
     struct rspline *rs[DIMSM];         /* data structure for use with smoopy   */
     struct sspline *ss[DIMSM];         /* data structure for use with sphery   */
     struct savedat *sd[DIMSM];         /* data structure for use with smoopy   */
     struct savedat_sphy *ssd[DIMSM];   /* data structure for use with sphery   */
     double  sm[DIMSM]={0., 0., 0.};    /* surface fitting smoothing parameters */

     GRID *gtemp=NULL;

     aniso = 'n';
     irmiss = 0;

     fptwrt = FPTWRT;

     strncpy(fout1, DATTHRO, MAXCHR);
     strncpy(fout2, DOUTOBJ, MAXCHR);
     strncpy(foutn, NEWOBJF, MAXCHR);
     strncpy(fintp, INTERP_TH, MAXCHR);

     if(iext){
        strcpy(strstr(fout1, EXTENSION), fext);
        strcpy(strstr(fout2, EXTENSION), fext);
        strcpy(strstr(foutn, EXTENSION), fext);
	strcpy(strstr(fintp, EXTENSION), fext);
     }

     alat = gr->alat;
     alng = gr->alng;

     for(i=0; i < DIMSM; i++){
        cokn[i] = NULL; rs[i] = NULL; ss[i] = NULL; sd[i] = NULL; ssd[i] = NULL;
     }

     dimp = gr->ix * gr->iy;

     xdm = xdim;
     ydm = ydim;
     xl = delb - x1u;
     xr = gr->ix - x2u;

     llb = (int * )calloc(gr->iy, sizeof(int));
     mem_er((llb == NULL) ? 0 : 1, gr->iy * sizeof(int));

     rrb = (int * )calloc(gr->iy, sizeof(int));
     mem_er((rrb == NULL) ? 0 : 1, gr->iy * sizeof(int));

     if(reo == 'y'){

        printf("do you want to use the original object/feature data (1)   \r\n"
               "or do you want to use data which has been feature filtered\r\n"
               "according to the max. permisable displacement (2)         \r\n"
               "                                                          \n\n");

        scanf("%d", &frd);

        switch(frd){
               case 1:
                  printf("do you want to use a different file to the default 'y' or 'n'\n");
                  scanf("\n");
                  if(getchar() == 'y') {
                        printf("input file name\n"); scanf("%s", charin); filpt = charin;}

                  else {

                     filpt = fout2;
                     printf("***INFORMATION***, reading from default file %s\n\n", filpt);


                  }

                  fobjo = open_file(filpt, "r");
                  break;

               case 2:
                  printf("do you want to use a different file to the default 'y' or 'n'\n");
                  scanf("\n");
                  if(getchar() == 'y'){
                       printf("input file name\n"); scanf("%s", charin); filpt = charin;}

                  else {

                     filpt = foutn;
                     printf("***INFORMATION***, reading from default file %s\n\n", filpt);

                  }

                  fobjo = open_file(filpt, "r");
                  break;

               default:
                  fobjo = NULL;
                  printf("incorrect key!\n");
                  exit(1);

        }

        if(OBJWRT){

           printf("do you want to read the object data, boundary data as well as the feature data, 'y' or 'n'\n");
           scanf("\n");
           if(!(getchar() == 'y')) objwr = 0;
           else objwr = 1;

        }

        fr = read_obd(fobjo, &gt->alat, &gt->alng, &pgr, &pty);

        i_shape = 1;

        close_file(fobjo, filpt);

        frold = fruser;

        if(pgr){
           if(gr == gr1) {alat = gr1->alat; alng = gr1->alng;}
           else if(gr == gr2){alat = gr2->alat; alng = gr2->alng; hty = 1;}

        }

        if((pgr != gt->prgr) || (pty != gt->prty) ||
          (fabs(gt->alat - alat) > TOLPROJ)  ||
          (fabs(gt->alng - alng) > TOLPROJ)  ) {

           printf("***WARNING***, map projection and data do not match\r\n"
                  "               grid may be corrected if possible.  \r\n"
                  "               If errors or faults occur check that\r\n"
                  "               input data corresponds with original\r\n"
                  "               values for data file creation.      \n\n");

           printf("Do you want to correct the grid, 'y' or 'n'.        \r\n"
                  "This will depend on the object data read in and if  \r\n"
                  "shape attributes are to be re-computed.             \r\n"
                  "Or do you want to correct the data to the grid, 'g' \r\n"
                  "This last option only works currently for the       \r\n"
                  "feature points, so if object point are present they \r\n"
                  "must match the chosen projection.                   \n\n");
           scanf("\n");
           dumi = getchar();
           if(dumi == 'y'){


              proj_group(pgr, pty);

              if(pgr) { 

                 if(!hty){

                    gr = gr1;
                    cm = cm1;

                 }

                 else {

                    gr = gr2;
                    cm = cm2;

                 }

                 x1u = 1;
                 x2u = gr->ix;
                 y1u = 1;
                 y2u = gr->iy;

                 xmn = *(gr->xgrid + x1u-1);
                 ymn = *(gr->ygrid + y1u-1);
                 xmx = *(gr->xgrid + x2u-1);
                 ymx = *(gr->ygrid + y2u-1);


                 printf("***WARNING***, region extent must be the same as when\r\n"
                        "               data file was created, check original \r\n"
                        "               input. Otherwise errors may result.  \n\n");


              }

           }

           else if (dumi == 'g') correct_fpt(fr, frold);

        }

        if(objwr){

           if(frd == 2) correct_fpt(fr, frold);

           printf("Do you want to display the object data, 'y' or 'n'\n\n");
           scanf("\n");
           if(getchar() == 'y'){

              printf("****INFORMATION****, the number of frames is %d\n\n", fruser);

              printf("What frames do you want plotted, input start, interval, end\n\n");
              scanf("%d %d %d", &rp_s, &rp_i, &rp_e);

              if((rp_s < 1 || rp_s > fruser) || (rp_e < 1 || rp_e > fruser) ||
                 (rp_i < 1 || rp_i > fruser)){

                  printf("****ERROR****, no plotting available for these frame identifiers.\n\n");
                  return fr;

              }
      

              for(i=rp_s - 1; i < rp_e; i+=rp_i) {

                  dfil = 0;

                  for(j=0; j<(fr+i)->obj_num; j++){

                     dfil += ((fr+i)->objs + j)->point_num;

                  }

                  filtd(fr+i, tf, fruser, i+1, thresh);

              }

           }

        }


        if(aniso != 'y'){
           printf("Do you want to use the feature point data to compute area and  \r\n"
	          "shape attributes, or add planetary waves back into the object  \r\n" 
		  "and feature point data? 'y' or 'n'                             \r\n\n"
		  "****WARNING****, it is the users responsibility to choose the correct frame sampling.               \n\n");
           scanf("\n");
           if(getchar() == 'n') return fr;
        }
	
	else {
           printf("Do you want to re-compute area and shape attributes, possibly based on a different field. 'y' or 'n'\r\n\n"
	          "****WARNING****, it is the users responsibility to choose the correct frame sampling.               \n\n");
           scanf("\n");
           if(getchar() == 'n') return fr;	   
 /*          else return fr; */
	   
	}

     }

     printf("do you want to scale the field 'y' or 'n'\n");

     scanf("\n");
     if((tscl = getchar()) == 'y'){printf("input scaling\n"); scanf("%f", &scale);}

     printf("do you want any offset subtraction from the data 'y' or 'n'?\n");
     scanf("\n");
     if((of = getchar()) == 'y'){

        printf("input data offset\n"); 
        scanf("%f", &offs);

        printf("will field values be required to be returned to original values\n");

        scanf("\n");
        if(getchar() == 'y') {trans[1] = 1; trn = 1;}

     }


     printf("input the required threshold \n");

     scanf("%f",&thresh);            /* read in threshold from st. input */


/* assign memory for field values */

     ap=(float *)calloc(dimp, sizeof(float));
     mem_er((ap == NULL) ? 0: 1, dimp * sizeof(float));

     if(gr->prgr && !(gt->iaz)){

        app=(float *)calloc(gt->ix * gt->iy, sizeof(float));
        mem_er((app == NULL) ? 0: 1, gt->ix * gt->iy * sizeof(float));
	

/* assign memory for surface fitting for interpolation onto new grid */

        printf("***INFORMATION***, you are going to interpolate/smooth the\r\n"
               "                   field data onto a new grid, probably   \r\n"
               "                   associated with a new projection.      \n\n");

        smty[1] = query_interp(0);

        cokn[1] = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
        mem_er((cokn[1] == NULL) ? 0 : 1, sizeof(struct sp_dat));

        if(smty[1]){

          ss[1] = (struct sspline * )malloc_initl(sizeof(struct sspline));
          mem_er((ss[1] == NULL) ? 0 : 1, sizeof(struct sspline));

          ssd[1] = (struct savedat_sphy * )malloc_initl(sizeof(struct savedat_sphy));
          mem_er((ssd[1] == NULL) ? 0 : 1, sizeof(struct savedat_sphy));

        }

        else{

          rs[1] = (struct rspline * )malloc_initl(sizeof(struct rspline));
          mem_er((rs[1] == NULL) ? 0 : 1, sizeof(struct rspline));

          sd[1] = (struct savedat * )malloc_initl(sizeof(struct savedat));
          mem_er((sd[1] == NULL) ? 0 : 1, sizeof(struct savedat));

        }

     }

     else app = ap;
     
/* setup for adding back into the object and feature points another field */
    
     printf("Do you want to add another field onto the object and feature point data,     \r\n"
            "for example, planetary waves that have been removed by filtering, 'y' or 'n'.\n\n");
	   
     scanf("\n");
     if(getchar() == 'y'){
     
       printf("Do you want to apply the same field adjustment to the object data, 'y' or 'n'\n\n");
       scanf("\n");
       adj_ob = getchar();
    
       aplw=(float *)calloc(dimp, sizeof(float));
       mem_er((aplw == NULL) ? 0: 1, dimp * sizeof(float));
       
       if(gr->prgr){

          applw=(float *)calloc(gt->ix * gt->iy, sizeof(float));
          mem_er((applw == NULL) ? 0: 1, gt->ix * gt->iy * sizeof(float));
	
       }
       
       else applw = aplw;
    
       printf("What file contains the additional data?\n\n");
       scanf("%s", filplw);
       
       if(form != 4) fplw = open_file(filplw, "r");
       else {
          fplw = (FILE *)nc_clone((NETCDF_INFO *)fdatin, filplw, NC_OPEN_MODE);
          ((NETCDF_INFO *)fplw)->iframe = 0;
       }
       
       if(of == 'y') {trans[1] = 1; trn = 1;}      
    
     }

/* set number of levels for data hierarchy */

     level = 1;

     if(sn == 0) {

        i = 1;
        while(xdim > (xddim = powi(2, i++))){};
        i = 1;
        while(ydim > (yddim = powi(2, i++))){};

     }

     else {xddim = xdim; yddim = ydim;}

     lmax = (xddim >= yddim) ? xddim : yddim;

     while(powi(2, level++) != lmax){};

     printf("The total number of levels for the data hierarchy is = %d\n\n", level);

/*     printf("Do you wish to work at a resolution different from the image's\n");
     scanf("\n");
     if(getchar() == 'y'){

        printf("What resolution do you want to work < levels\n\n");
        scanf("%d", &lr);

     } */

/* assign tempory storage for binary map */

     ip=(struct image *)calloc(xddim*yddim, sizeof(struct image));
     mem_er((ip == NULL) ? 0 : 1, xddim*yddim * sizeof(struct image));

/* assign work space for image hierarchy */

     hierarchy = (struct image ** )calloc(level, sizeof(struct image *));
     mem_er((hierarchy == NULL) ? 0 : 1, level * sizeof(struct image *));

     printf("do you want MAX (input '1') or MIN (input '0') thresholding\n\n");
     scanf("%d", &type);

     switch(type){
            case TYPMIN: 

                  sign=-1.0;
                  break;

            case TYPMAX:

                  sign=1.0;
                  break;

            default:

                  printf("Key does not match for thresholding.\n");
                  exit(1);
                  break;                 

     }


     if(!type){


           if(!fplw){
              printf("will field values be required to be returned to original values\n");

              scanf("\n");
              if(getchar() == 'y') {trans[0] = 1; trn = 1;}
	   }
	   
	   else {trans[0] = 1; trn = 1;}

     }

/* is inversion of a hemispherical part of the field required e.g. relative vorticity */


     if(!(gr->prgr)){

        printf("do you want to invert the field for a hemisphere?\r\n"
               "required for fields such as the relative vorticity.\r\n"
               "input 'y' or 'n'\n");

        scanf("\n");

        if(getchar() == 'y'){

           printf("input 'n' for northern and 's' for southern\n");
           scanf("\n");
           hemi = getchar();
 
           printf("input the lat. about which to invert\n");
           scanf("%f", &ilat);
	   
	   if(!fplw){

              printf("will field values be required to be returned to original values\n");
 
              scanf("\n");
              if(getchar() == 'y') {trans[2] = 1; trn = 1;}
	   
	   }
	   else {trans[2] = 1; trn = 1;}

        }

     }

     

/* define starting frame and framing interval */

     printf("the available number of frames is = %d\n\n", frnum);

     printf("which frame do you wish to start from\n\n");
     scanf("%d", &ffirst);

     if(ffirst <= 0 || ffirst > frnum) {
        
        printf(" ***error*** no such frame number\n\n");
        exit(1);
        
     }


     printf("what frame interval do you require\n\n");
     scanf("%d", &fint);
     
     if(fint < 1 || fint > frnum){
        printf("***ERROR***, frame interval not valid.\n\n");
	exit(1);
     }

     printf("what frame do you wish to terminate on\n\n");
     scanf("%d", &frterm);

     if(frterm < ffirst) {
        printf("***ERROR*** termination frame index is less than starting frame index\n\n");
	exit(1);
     }

/* define what kind of connectivity search is required */

     printf("what type of connectivity search is required \r\n"
            "for vertex and edge connectivity input v (not V) \r\n"
            "for edge connectivity only input e (not E) \n");

     scanf("\n");
     cc = getchar();

/* is border searching required for rest of objects outside region of
   interest and/or to satisfy continuity at periodic boundarys        */

   if(!(gr->prgr)){

      printf("is a boundary search required for parts of objects?\r\n"
             "answer y or n\n");

      scanf("\n");
      bs = getchar();

   }

   else bs = 'n';

/* are periodic X-boudaries needed */

   if(bs == 'y') {

      printf("are the X-boundarys periodic, 'y' or 'n'\n");

      scanf("\n");

      pb=getchar();

   }

/* define the filtering threshold for object size */

    printf("what is the lower limiting size of objects, \r\n"
           "input as a number of points!\n");

    scanf("%d", &filt);


/* input the type of feature points required */

   printf("what kind of feature points are required \r\n"
          "   Input a '0' for no feature points, threshold display only \r\n"
          "   Input a '3' for local max's\r\n"
          "   Input a '4' for a connectivity search of local max's \r\n"
          "   Input a '7' for surface fitting of region of interest and local optimization  \r\n"
          "   Input a '8' for surface fitting to objects and local optimization             \r\n"
          "               and the determination of shape parameters.                      \r\n\n"
          " # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #       \b\r\n\n"
          "   ***WARNING***, last two options must have surface fitting and/or              \r\n"
          "                  local optimization libraries linked,                           \r\n"
          "                  this is the users responsibility.                              \n\n"
          " # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #           \n\n"
          "   Input a '9' for the distance transform method.                                \r\n" 
          "               (useful for objects without max. or min.).                        \n\n");


   scanf("%d", &tf);

/* are the shape measures, anisotropy and orientation required */

   if(!(tf == 7 || tf == 9)) {

      printf("Do you want boundary determination, 'y' or 'n'\r\n");

      scanf("\n");
      bf = getchar();

      if(bf == 'y') shape_setup(&bcntl, 1);

   }

   if(tf){
      printf("do you want to compute the anisotropy, orientation and area, 'y' or 'n'?\n\n");
      scanf("\n");
      aniso = getchar();
      if(fr) tf = 4;
   }

   if(!tf) fptwrt = 0;

/* assign work space if options >=4 are chosen. */

   if(tf >= 4){

      fp=(struct frame_objs *)malloc_initl(sizeof(struct frame_objs));
      mem_er((fp == NULL) ? 0 : 1, sizeof(struct frame_objs));

      fp->obj_num=0;
      fp->objs = NULL;
      fp->tot_f_f_num = 0;
      fp->b_state = -1;

    }

    if(tf == 7 || (fplw && tf > 3)){

       smty[0] = query_interp(1);

       if(gr->prty || gr->prgr) smty[0] = 0;

       cokn[0] = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
       mem_er((cokn[0] == NULL) ? 0 : 1, sizeof(struct sp_dat));

       if(smty[0]){

          ss[0] = (struct sspline * )malloc_initl(sizeof(struct sspline));
          mem_er((ss[0] == NULL) ? 0 : 1, sizeof(struct sspline));

          ssd[0] = (struct savedat_sphy * )malloc_initl(sizeof(struct savedat_sphy));
          mem_er((ssd[0] == NULL) ? 0 : 1, sizeof(struct savedat_sphy));

       }

       else{

          rs[0] = (struct rspline * )malloc_initl(sizeof(struct rspline));
          mem_er((rs[0] == NULL) ? 0 : 1, sizeof(struct rspline));

          sd[0] = (struct savedat * )malloc_initl(sizeof(struct savedat));
          mem_er((sd[0] == NULL) ? 0 : 1, sizeof(struct savedat));

       }

    }

    if(tf == 8){

        smty[2] = 0;
        cokn[2] = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
        mem_er((cokn[2] == NULL) ? 0 : 1, sizeof(struct sp_dat));

        cokn[2]->sm_type = 0;

        rs[2] = (struct rspline * )malloc_initl(sizeof(struct rspline));
        mem_er((rs[2] == NULL) ? 0 : 1, sizeof(struct rspline));

        sd[2] = (struct savedat * )malloc_initl(sizeof(struct savedat));
        mem_er((sd[2] == NULL) ? 0 : 1, sizeof(struct savedat));

        printf("Do you want object shape smoothing before smoothing \r\n"
               "over the objects, 'y' or 'n'.                       \n\n");

        scanf("\n");
        ssm = getchar();

    }

    if(tf){
       printf("Do you want to filter the object feature points for un-physical values? 'y' or 'n'\n\n");
       scanf("\n");
       unp_filt = getchar();

       if(unp_filt == 'y'){
         printf("What lower and upper intensity values do you want?\n\n");
         scanf("%f %f", &upint1, &upint2);
       }

       printf("Do you want to filter the object feature points, to retain the\r\n"
              "point with the largest value (depending on feature detection  \r\n"
              "scheme), 'y' or 'n'\n\n");
       scanf("\n");
       mfilt = getchar();

       if(mfilt != 'y'){
          printf("Do you want to filter the object feature points for points too close together, 'y' or 'n'.\n\n");
          scanf("\n");
          if((mfiltd = getchar()) == 'y'){
             printf("What is the search radius required to find feature points too close together?\n\n");
             scanf("%f", &rads);

             if(tom == 'g') rads *= FP_PI;
             
          }
       }

    }
    
/* perform time average subtraction or thresholding if required */

    printf("do you want time average subtraction from the data\r\n"
           "or time average thresholding                      \r\n"
           "input 'y' (not 'Y') for yes and                   \r\n"
           "input 'n' (not 'N') for no.                       \n\n");

    scanf("\n");
    as = getchar();
    
    if(as == 'y' && fplw){
       printf("****WARNING****, time average subtraction not available together with \r\n"
              "                 additional data addition. Time subtraction ignored.  \n\n");
       as = 'n';
    }

    if(as == 'y'){

      while(tsth != 's' || tsth != 't'){

            printf("which type of time average usage do you want?\r\n"
                   "input 's' for subtraction, or                \r\n"
                   "input 't' for thresholding.                  \n\n");

            scanf("\n");
            tsth = getchar();

            if(tsth == 't'){

               printf("What offset to the time average field is required     \r\n"
                      "for time average thresholding, input as a percentage. \r\n"
                      "The threshold is then TIME_AVERAGE +/- X%%.           \n\n");

               scanf("%f", &tavoff);

               tavoff /= 100.;

            }

            if((tsth == 's') || (tsth == 't')) break;

            else 
               printf("Incorrect input of descriptor try again\n\n");

      }

      ata = (float *)calloc(gr->ix * gr->iy, sizeof(float));
      mem_er((ata == NULL) ? 0 : 1, gr->ix * gr->iy * sizeof(float));

      if(gr->prgr && tsth == 't' && !(gt->iaz)){

        ataa=(float *)calloc(gt->ix * gt->iy, sizeof(float));
        mem_er((ataa == NULL) ? 0: 1, gt->ix * gt->iy * sizeof(float));
      }

      else ataa = ata;

      gtemp = gr;
      gr = gt;

      strncpy(tavin, DTIMAVG, MAXCHR);

      if(!fexist(tavin, "r")){

         strncpy(tavin, TAVGE, MAXCHR);
         if(iext) strcpy(strstr(tavin, EXTENSION), fext);
         if(fexist(tavin, "r")){ 

            printf("***WARNING***, the time average file %s exists    \r\n"
                   "               do you want to use this file       \r\n"
                   "               , 'y' or 'n'.                      \r\n",tavin);
            scanf("\n");
            if(getchar() == 'n') {  

               printf("***WARNING***, the time avearge file %s \r\n"
                      "               cannot be read.          \r\n"
                      "               Creating a new time average file\n\n", DTIMAVG);


              if(form != 4) place1=ftello(fdatin);


              time_avg(fdatin, ffirst, fint, frterm, 1);


              if(form != 4) fseeko(fdatin, place1, FSTART);

            }

         }

         else {


            printf("***WARNING***, the time avearge file %s \r\n"
                   "               cannot be read.          \r\n"
                   "               Creating a new time average file\n\n", DTIMAVG);


            if(form != 4) place1=ftello(fdatin);


            time_avg(fdatin, ffirst, fint, frterm, 1);


            if(form != 4) fseeko(fdatin, place1, FSTART);

         }

      }


      if(form != 4) fta = open_file(tavin, "r");
      else {
         fta = (FILE *)nc_clone((NETCDF_INFO *)fdatin, tavin, NC_OPEN_MODE);
         ((NETCDF_INFO *)fta)->iframe = 0;
      }

      switch(tsth){
          case 's':
             read_field(ata, ata, scale, fta, 1, 'n', 'n', '0', tscl);
             break;
          case 't':
             read_field(ataa, ata, scale, fta, sign, 'n', of, hemi, tscl);
             as = 'n';
             break;
      }

      if(form != 4)close_file(fta, DTIMAVG);
      else netcdf_close((NETCDF_INFO *)fta);

      gr = gtemp;
      gtemp = NULL;

      if(gr->prgr && tsth == 't' && !(gt->iaz)){

         printf("***INFORMATION***, interpolating time average field.\n\n");

         apt = ap;
         ap = ata;
         if(smty[1]) proj_interp(ataa, &sm[1], 0, smty[1], cokn[1], ss[1], ssd[1]);
         else proj_interp(ataa, &sm[1], 0, smty[1], cokn[1], rs[1], sd[1]);

         smsw = 1;

         ap = apt;
         apt = NULL;

      }

/* apply offset to time average */

      if(tsth == 't'){

         for(i=0; i< dimp; i++) *(ata + i) += *(ata + i) * tavoff;


      } 

    }

/* calculate data block size for each frame */

    if(form != 4){
       place1=ftello(fdatin);
       if(fplw) fseeko(fplw, place1, FSTART);  
    }
    else ((NETCDF_INFO *)fdatin)->iframe = 0;

    gtemp = gr;
    gr = gt;

    rf = read_field(app, ata, scale, fdatin, sign, as, of, hemi, tscl);
    frtim1 = frtime;
    if(frtime){
       printf("Specify time step for missing frame search.\n\n");
       scanf("%d", &tstep);
       irmiss = 1;
       printf("****INFORMATION****, a missing frame search will be performed.\n\n");
    }

    if(fplw) rpl = read_field(applw, ata, scale, fplw, 1, 'n', 'n', '0', tscl); 
    
    gr = gtemp;
    gtemp = NULL;

    if(gr->prgr && !(gt->iaz)){

       if(smty[1]) proj_interp(app, &sm[1], smsw, smty[1], cokn[1], ss[1], ssd[1]);
       else proj_interp(app, &sm[1], smsw, smty[1], cokn[1], rs[1], sd[1]);

       smsw = 1;
       
       if(fplw){
          apt = ap;
	  ap = aplw;
          if(smty[1]) proj_interp(applw, &sm[1], smsw, smty[1], cokn[1], ss[1], ssd[1]);
          else proj_interp(applw, &sm[1], smsw, smty[1], cokn[1], rs[1], sd[1]);
	  ap = apt;
	  apt = NULL;
       }

    }

    if(form != 4) {place2=ftello(fdatin); chrnum=place2-place1;}

/* find starting frame */

     if(ffirst > 1) {

       if(form != 4){
          fseeko(fdatin,(ffirst-2)*chrnum, ORIGIN);
	  if(fplw) fseeko(fplw,(ffirst-2)*chrnum, ORIGIN);  
       }
       else{
          ((NETCDF_INFO *)fdatin)->iframe = ffirst - 1;
          if(fplw)((NETCDF_INFO *)fplw)->iframe = ffirst - 1;	  
       }

       gtemp = gr;
       gr = gt;

       rf = read_field(app, ata, scale, fdatin, sign, as, of, hemi, tscl); 
       frtim1 = frtime;
       
       if(fplw) rpl = read_field(applw, ata, scale, fplw, 1, 'n', 'n', '0', tscl);
       
       gr = gtemp;       
       gtemp = NULL;
    
       if(gr->prgr && !(gt->iaz)){

          if(smty[1]) proj_interp(app, &sm[1], smsw, smty[1], cokn[1], ss[1], ssd[1]);
          else proj_interp(app, &sm[1], smsw, smty[1], cokn[1], rs[1], sd[1]);
	  
	  if(fplw){
             apt = ap;
	     ap = aplw;
             if(smty[1]) proj_interp(applw, &sm[1], smsw, smty[1], cokn[1], ss[1], ssd[1]);
             else proj_interp(applw, &sm[1], smsw, smty[1], cokn[1], rs[1], sd[1]);
	     ap = apt;
	     apt = NULL;
          }

       }

     }
     
     if(INTERP_DAT){
        printf("****INFORMATION****, data (no threshold) will be written to new file on current projection.\n\n");
	
	finp = open_file(fintp, "w");

/* write header */

        write_header(gr, finp);
     
     }

/* read in the frame data and write to output file */

     fruser=1;
     frcount=ffirst;

/* open output file for thresholded data */

/*     fthro = open_file(file_exist(fout1, "r", APPS), "w");

     fprintf(fthro, "%d %d %d %d\n", x1u, x2u, y1u, y2u);  */
     
     if(INTERP_DAT){
        fprintf(finp, "FRAME %6d\n", fruser);
	fwrite(ap, dimp * sizeof(float), 1, finp);
	fprintf(finp, "\n");
     }

/* begin main image processing for first frame */

     arrayd(ip, ap, ata, fthro, thresh, ffirst, xddim, tsth); 

     fo=(struct frame_objs *)malloc_initl(sizeof(struct frame_objs));
     mem_er((fo == NULL) ? 0 : 1, sizeof(struct frame_objs));

     fo->nmiss = 0;
     fo->obj_num=0;
     fo->objs = NULL;
     fo->tot_f_f_num = 0;
     fo->b_state = -1;
     fo->frame_id = frcount;

     hierarc_segment(hierarchy, ip, level, fo, xddim, yddim, cc);

/* perform boundary search */

     if(bs == 'y') {

        border_expand_object(fo, xdim, ydim, thresh);

        if(pb == 'y') periodic_boundary(fo, llb, rrb);

     }

/* filter objects */

     object_filter(fo, filt, xdim, ydim);

/* determine feature points */

     if(!fr){

/*        if(tf == 1 || tf == 2) centroid(fo); */
        if(tf >= 3){

          if(tf == 8){

             if(ssm == 'y') object_dist_trans(fo, fp, llb, rrb, 1, 0);
             object_smint(fo, &sm[2], sd[2], rs[2], cokn[2], 0);


          }

          else if(tf == 9)object_dist_trans(fo, fp, llb, rrb, 1, 1);


          else {


             printf("When finding grid local maxima do you want to exclude   \r\n"
                    "boundary maxima, e.g.. maxima adjacent to threshold     \r\n"
                    "boundary or maxima adjacent to a land boundary for ocean\r\n"
                    "problems?                                               \r\n"
                    "Input '1' for exlusion, '0' otherwize.                  \n\n"); 

             scanf("%d", &b_exc);

             object_local_maxs(fo, fp, tf, b_exc);

          }

/* fit surface to region of interest if required */

          if(tf == 7){

             if(pb == 'y') extnpts_to_domain(fo, llb, rrb);

             if(smty[0]) surfit(&sm[0], itpfc, smty[0], cokn[0], ss[0], ssd[0]);
             else surfit(&sm[0], itpfc, smty[0], cokn[0], rs[0], sd[0]);
	     itpfc = 1;

             non_lin_opt(fo, 0, smty[0]);

          }

          if(tf != 7 ){

             if(bf == 'y') frame_boundaries(fo, &bcntl, llb, rrb);

          }


        }

        else{

           if(bf == 'y') frame_boundaries(fo, &bcntl, llb, rrb);

        }

     }

     else {

        printf("When finding grid local maxima do you want to exclude   \r\n"
               "boundary maxima, e.g.. maxima adjacent to threshold     \r\n"
               "boundary or maxima adjacent to a land boundary for ocean\r\n"
               "problems?                                               \r\n"
               "Input '1' for exlusion, '0' otherwize.                  \n\n"); 

        scanf("%d", &b_exc);

        object_local_maxs(fo, fp, tf, b_exc);
     }

/* match feature points with objects */


     if(fr && tf && fruser <= frold) feature_to_object(fo, fr, frd);

     if(aniso == 'y') {
        apt = ap;
        ap = app;
        anisotropy(fo, 0);
	ap = apt;
	apt = NULL;
        if(fr) putback_area(fo, fr);
     }

     ft = (fr) ? fr : fo;

     if(tf != 7 && pb =='y') extnpts_to_domain(ft, llb, rrb);

/* open output file for object data */

     fobjo = open_file(file_exist(fout2, "r", APPS), "w");

     if(aniso == 'y') fprintf(fobjo, "%d %d %d\n", tf, 1, irmiss);
     else fprintf(fobjo, "%d %d %d\n", tf, 0, irmiss);

     fprintf(fobjo, "PROJ_DETAILS\n");
     fprintf(fobjo, "%d %d\n", gr->prgr, gr->prty);
     if(gr->prgr) fprintf(fobjo, "%f %f\n", gr->alat, gr->alng);
     fprintf(fobjo, "REGION_DETAILS\n");
     fprintf(fobjo, "%d %d %d %d\n", x1u, x2u, y1u, y2u);
     fprintf(fobjo, "ADDITIONAL_FIELDS\n");
     fprintf(fobjo, "%3d %3d &", nf, nfld);
     if(nfwpos){
        for(i=0; i < nf; i++)fprintf(fobjo, "%1d", *(nfwpos + i));
     }
     fprintf(fobjo, "\n");

     if(unp_filt == 'y') fpt_intensity(ft, upint1, upint2);

     if(mfilt == 'y') mfilt_obj(ft);
     else if(mfiltd == 'y') ft = mfilt_obj_dist(ft, rads, tf);

     if(trn) object_realf(ft, sign, of, hemi, tf, trans);
     
/* add additional field values back in to the feature points */
     
     if(fplw){
       if(tf > 3){
          apt = ap;
          ap = aplw;
          if(smty[0]) surfit(&sm[0], itpfc, smty[0], cokn[0], ss[0], ssd[0]);
          else surfit(&sm[0], itpfc, smty[0], cokn[0], rs[0], sd[0]);
          ap = apt;
          apt = NULL;
          itpfc = 1;
       }
       adjust_fpt(ft, cokn[0], aplw, adj_ob);    
     }

     if(OBJWRT){

        printf("Do you want object data written to file as well as feature data, 'y' or 'n',\r\n"
               "this include boundary data if computed.                                     \n\n");
        scanf("\n");
        if(!(getchar() == 'y')) objwr = 0;

     }

/* delay writing object data if there maybe missing frames */

     if(frtime){

        printf("do you want to delay writing object data till the end of the identification stage, 'y' or 'n'\n\n");
        scanf("\n");
        delobd = getchar();

     }

/* plot first user chosen frame if asked for */

     printf("do you want the first chosen frame plotted, answer y or n\n");

     scanf("\n");
     c = getchar();
 
     if(c == 'y') filtd(ft, tf, fruser, frcount, thresh);

     if(delobd == 'n') objectd(ft, fobjo, frcount, tf);

     printf(" do you want any other user chosen frames plotted?\n");
     printf("input the multiplication factor of the frame interval you require\n");
     printf("input a '0' if no further plotting is required\n");


     scanf("%d", &multi); 

/* free local feature points pointer if no further plotting is required */

     frcount += fint;

     while(frcount<=frnum && !rf && !rpl) {

          
          if(form != 4) fseeko(fdatin, (fint-1)*chrnum, ORIGIN); 
          else ((NETCDF_INFO *)fdatin)->iframe = frcount - 1;

          gtemp = gr;
          gr = gt;

          rf = read_field(app, ata, scale, fdatin, sign, as, of, hemi, tscl);
          frtim2 = frtime;
	  
	  if(fplw) rpl = read_field(applw, ata, scale, fplw, 1, 'n', 'n', '0', tscl);

          gr = gtemp;
          gtemp = NULL;

          if(rf) break;

          ++fruser;

          printf("Processing frame %d \n", frcount);
 
          if(gr->prgr && !(gt->iaz)){

            if(smty[1]) proj_interp(app, &sm[1], smsw, smty[1], cokn[1], ss[1], ssd[1]);
            else proj_interp(app, &sm[1], smsw, smty[1], cokn[1], rs[1], sd[1]);
	    
	    if(fplw){
               apt = ap;
	       ap = aplw;
               if(smty[1]) proj_interp(applw, &sm[1], smsw, smty[1], cokn[1], ss[1], ssd[1]);
               else proj_interp(applw, &sm[1], smsw, smty[1], cokn[1], rs[1], sd[1]);
	       ap = apt;
	       apt = NULL;
            }

          }

          pflag = 0;

          if(multi > 0){

             multis=ffirst;

             while(multis <= frnum){ 

                if(frcount == multis) pflag=1;

                multis+=multi*fint;

             }
          }
	  
	  
	  if(INTERP_DAT){
             fprintf(finp, "FRAME %6d\n", fruser);
	     fwrite(ap, dimp * sizeof(float), 1, finp);
	     fprintf(finp, "\n");
          }

/* process each chosen frame */

          arrayd(ip, ap, ata, fthro, thresh, frcount, xddim, tsth);

          if(tf) {

             if(!(fptwrt)){

                fo=(struct frame_objs *)realloc_n(fo, fruser*sizeof(struct frame_objs));
                mem_er((fo == NULL) ? 0 : 1, fruser*sizeof(struct frame_objs));

                ff = fo+fruser-1;
                ff->nmiss = 0;

/* missing frame search */

                if(frtime){
                   nmiss = 0;
                   while(frtim2 > (frtim1 = new_time(frtim1, tstep))) ++nmiss;
                   (ff - 1)->nmiss = nmiss;
                   if(nmiss) {
                      printf("****WARNING***, missing frames have been detected.\n\n");
                   }

                }

             }

             else ff = fo;

          }

          else ff = fo; 

          ff->obj_num=0;
          ff->objs = NULL;
          ff->tot_f_f_num = 0;
          ff->b_state = -1;
          ff->frame_id = frcount;

          hierarc_segment(hierarchy, ip, level, ff, xddim, yddim, cc);

          if(bs == 'y') {

             border_expand_object(ff, xdim, ydim, thresh);

             if(pb == 'y') periodic_boundary(ff, llb, rrb);

          }

          object_filter(ff, filt, xdim, ydim);

          if(!fr){

/*             if(tf == 1 || tf == 2) centroid(ff); */
             if(tf >= 3){
  
                if(tf == 8){

                  if(ssm == 'y') object_dist_trans(ff, fp, llb, rrb, 0, 0);
                  object_smint(ff, &sm[2], sd[2], rs[2], cokn[2], 1);

 
                }

                else if(tf == 9)object_dist_trans(ff, fp, llb, rrb, 0, 1);


                else object_local_maxs(ff, fp, tf, b_exc);

                if(tf == 7){

                  if(pb == 'y') extnpts_to_domain(ff, llb, rrb);

                  if(smty[0]) surfit(&sm[0], itpfc, smty[0], cokn[0], ss[0], ssd[0]);
                  else surfit(&sm[0], itpfc, smty[0], cokn[0], rs[0], sd[0]);

                  non_lin_opt(ff, 1, smty[0]);

                }

                if(tf != 7){

                   if(bf == 'y') frame_boundaries(ff, &bcntl, llb, rrb);

                }


             }

             else{

               if(bf == 'y') frame_boundaries(ff, &bcntl, llb, rrb);

             }

          }

          else object_local_maxs(ff, fp, tf, b_exc);

          if(fr && tf && fruser <= frold) feature_to_object(ff, fr + fruser - 1, frd);

          if(aniso == 'y') {
	    apt = ap;
            ap = app;
            anisotropy(ff, 1);
	    ap = apt;
	    apt = NULL;
            if(fr) putback_area(ff, fr + fruser - 1);
          }

          ft = (fr) ? fr + fruser - 1 : ff;

          if(tf != 7 && pb =='y') extnpts_to_domain(ft, llb, rrb);

          if(unp_filt == 'y') fpt_intensity(ft, upint1, upint2);

          if(mfilt == 'y') mfilt_obj(ft);
          else if(mfiltd == 'y') ft = mfilt_obj_dist(ft, rads, tf);

          if(trn) object_realf(ft, sign, of, hemi, tf, trans);
	  
	  if(fplw){
	     if(tf > 3){
                apt = ap;
                ap = aplw;
                if(smty[0]) surfit(&sm[0], itpfc, smty[0], cokn[0], ss[0], ssd[0]);
                else surfit(&sm[0], itpfc, smty[0], cokn[0], rs[0], sd[0]);
                ap = apt;
                apt = NULL;
	     }
	     adjust_fpt(ft, cokn[0], aplw, adj_ob);   
          }

          if(pflag == 1){

            filtd(ft, tf, fruser, frcount, thresh);

          }

          if(delobd == 'n') objectd(ft, fobjo, frcount, tf);

          frcount += fint;

          frtim1 = frtim2;

          if(frcount > frterm) break;
  
     }
     
     
     if(rf != rpl){
        printf("****WARNING****, number of frames in data file and additional data file do not match.\n\n");
     }
     
/* close file for writing gridded data if used */

     if(INTERP_DAT) {
        fseeko(finp, (off_t)0, FSTART);
	fprintf(finp, "%6d %6d %6d\n", gr->ix, gr->iy, fruser);
        close_file(finp, fintp);
     }

     frcnt = frcount - fint;

/* close output file for writing object data */ 

    if(delobd == 'y') {

       for(i=0; i < fruser; i++) {

           ft = (fr) ? fr + i : fo + i;
           objectd(ft, fobjo, ft->frame_id, tf);

       }

    }

    close_file(fobjo, fout2);

    if(fr && fruser != frold) {

       printf("****ERROR****, incompatable number of frames processed compared to number read from file, %d != %d\n\n", frold, fruser);
       exit(1);

    }
    
/* if file is open for reading additional fields then close and free memory */

    if(fplw){
    
       if(form != 4)close_file(fplw, filplw);
       else netcdf_close((NETCDF_INFO *)fplw);
          
       free(aplw);
       free(applw);
    }

/* close output file for writing thresholded data */ 

/*    close_file(fthro, fout1); */

    if(tf == 7) {

      if(smty[0]) {
        surfit(&sm[0], -1, smty[0], cokn[0], ss[0], ssd[0]);
        free(ss[0]);
        free(ssd[0]);
      }
      else {
        surfit(&sm[0], -1, smty[0], cokn[0], rs[0], sd[0]);
        free(rs[0]);
        free(sd[0]);
      }

      free(cokn[0]); 
  
      non_lin_opt(NULL, -1, smty[0]);

    }

    else if(tf == 8){

        object_smint(NULL, &sm[2], sd[2], rs[2], cokn[2], -1);

        free(rs[2]);
        free(sd[2]);
        free(cokn[2]);

    }

/* de-allocate storage for the time average */

    free(ata);
    free(ataa);

/* de-allocate hierarchy level storage */


    free(hierarchy);

/* remove binary map storage */

    free(ip);

/* remove field storage */

    free(ap);

    if(gr->prgr && !(gt->iaz)){


       if(smty[1]){

          surfit(&sm[1], -1, smty[1], cokn[1], ss[1], ssd[1]);
          free(ss[1]);
          free(ssd[1]);

       }

       else {

          surfit(&sm[1], -1, smty[1], cokn[1], rs[1], sd[1]);
          free(rs[1]);
          free(sd[1]);

       }

       free(app);
       free(cokn[1]);

     }

/* remove storage for tempory object/feature point data */

     if(fr){

        if(tf == 3) tf = 4;

        for(i=0; i < fruser; i++){

            ff = fo + i;

            if((ff->objs)){

               for(j=0; j < ff->obj_num; j++){

                   ob = ff->objs + j;
                   free(ob->pt);
                   if(ob->fet) free(ob->fet->fpt);
                   free(ob->fet);
                   free(ob->ext);
                   free(ob->bound);
               }

               free(ff->objs);

            }

        }

        fo = fr;

     }

/* remove storage for area calculation */

     sys_area(NULL, 0, -1, 0, 0.0, NULL);


/* remove remaining tempory storage */

     if(tf >= 4) free(fp);

/* on completion of segmetation periodic boundary condition
   should not be in force.                                   */

     if(pb == 'y') pb = 'n';

     free(llb);
     free(rrb);


     free(chrfld);

     if(fptwrt) {

       free(fo);

       fobjo = open_file(fout2, "r");

       fo = read_obd(fobjo, &dumf, &dumf, &dumi, &dumi);

       close_file(fobjo, fout2);
       fptwrt = 0;

     }


     return fo;

}

