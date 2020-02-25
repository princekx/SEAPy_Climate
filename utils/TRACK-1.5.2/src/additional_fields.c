#include <Stdio.h>
#include <stdlib.h>

#ifndef  REGULAR

int nf=0, nfld=0;
int *nfwpos=NULL;
void additional_fields()

{

   printf("***error***, surface fitting impossible unless correct libraries\r\n"
          "             are linked, see compilation options.               \n\n");
   exit(1);

}

#else

#include <Math.h>
#include <string.h>
#include <sys/types.h>
#include "splice.h"
#include "grid.h"
#include "file_handle.h"
#include "file_cat_out.h"
#include "files_out.h"
#include "mem_er.h"
#include "reg_dat.h"
#include "sphery_dat.h"
#include "bisp.h"
#include "netcdf_info.h"
#include "m_values.h"
#include "st_obj.h"
#include "st_fo.h"
#include "vec.h"
#include "proj.h"


#define  TOLPER   1.0e-6 /* tolerance on points on periodic boundary. */

/*#define  LARGE    1.0e+12*/
#define  LARGE    ADD_UNDEF
#define  TOLAREA  1.0e-6
#define  LNGOFF   20.0
#define  MISSV    1.0e+3
#define  LINT     100000000

/* function to compute values from alternative fields at track feature points */

#ifdef  NOUNDERSCORE

void bisp(double * , double * , double * , int * , double * , int * , 
          double * , int * , double * , double * , int * );
	  
#else

void bisp_(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#endif

void surfit(double * , int , int , ... );
int query_interp(int );
int read_field(float * ,float * , float , FILE * , int , int , int , int , int );
void meantrd(FILE * , struct tot_tr * , int , int );
long int new_time(long int , int );
void non_lin_opt(struct frame_objs * , int , int );
void detext(struct object * , float , int *);
int *nln(GRID * , float );
void geo_conv(struct feature_pts * );
float measure(struct feature_pts * , struct feature_pts * );
float fmaxmin(struct feature_pts * , struct feature_pts * , struct extent * , float , float * , float * , float * , float * , float , float , float , int , int , int , int , int , int , int, int , int , GRDP * , GRDP *);
float maxmin(struct feature_pts * , float , float , float * , float * , int );
void integ(double * , double * , float , int , int , int );
void convert_track(struct tot_tr * , int , int , int );
int find_last(struct fet_pt_tr * , struct feature_pts * , int * , int * , int , int );
double fvec(VEC * , VEC * , VEC * , VEC * );
double dirangle(struct tot_tr * , double * , double * , int , int , int , int );
void vec_offset(struct feature_pts * , double , double , double , double );
int missing(float, float, int );
void smoopy_setup(int );
void sphery_setup(int );

void projext(struct object * , double , double , double * , double * , int , GRDP * , GRDP *);

float nearest_neigbour_intrp(float * , float , float , float , int );
float bilinear_intrp(float * , float , float , float , int );
float bicubic_intrp(float * , float , float , float , int );

int create_mask(int * , float * , float , int);
void writef(FILE * , float * , int );
void write_header(GRID * , FILE * );
int mask_data(int * , float , float , int );

int true_maxmin(float , float, int , int , int , int , int , int );


extern GRID *gr, *gr1;
extern int form;
extern int tl, gof;
extern int frnum;
extern float *ap;
extern float period;
extern int tom;
extern int pb;
extern int x1u, x2u, y1u, y2u;
extern int tf;
extern int geo_init;
extern int iext;
extern char *fext;

extern float xmn, ymn, xmx, ymx;

extern PROJ *pp;

int nf=0, nfld=0;
int *nfwpos=NULL;

void additional_fields(struct tot_tr *trs, int trnum, off_t pl, FILE *fst, int iper, float cset)
{

   int i, j, k;
   int *nfst=NULL, *ninv=NULL, *nlst=NULL;
   int smty=0;
   int frc=0, rf=0;
   int itpfc=0;
   int idif=0;
   int ittype=0;
   int tstep=6;
   int iin=0;
   int inmax=0, *inm=NULL;
   int frst=0;
   int uscl='n';
   int *nl=NULL, *nl_l=NULL;
   int wnpos=0;
   int nwfld=0, nofld=0, nff=0, nnfld=0, nxfl=0;
   int liter=0;
   int ii=0, kk=0;
   int nffld=0, nadf=0, nmnf=0;
   int ifnd=0;
   int mm_type=0;
   int idiag=0, idir=0, idiradd=0, ifdir=0, ifdirp=0;
   int nr=0, ntheta=0, npcy=0;
   int pt_id=0;
   int imnfr=10000000, imxfr=0;
   int ifrch=0, isgn=1;
   int ia_cnt=0, i_av_var=0, ia_thr=0;
   int ifcnt=0, *ifncnt=NULL, *ifncntp=NULL;
   int itpadd=0;
   int ianom=0, ianom_l=0, iprojd=0;
   int ndsmth=0;
   int nloc=0, ilocd=0;
   int ipoff=0;
   int imiss=0, icmp=-1;
   int iintp=1;
   int inewcor=0, icswap=0;
   int itmxmn=0;
   int iarea=0;
   int igsmp_typ=0;
   int nll=0, nll2=0;
   int itrtr = 0;
   int nfc=0;
   
   int nfo=0, nfldo=0;
   int ns=0;

   int igm='n', igmn=0;
   int **imask=NULL, *imy=NULL, *imk=NULL;
   int nmask=0, imyy=0;
   int dim=0;

   int intrp_type=0;
   int icirc=1;
   int iproj=0;

   long int ptnum=0;

   long int ttime=0, stime=0;
   long int i1, i2, fruser=0;

   off_t place1=0, place2=0, chrnum=0;
   off_t **fldpl=NULL, nfll=0;
   
   float *scale=NULL;
   float *app=NULL, *apn=NULL, *atmp=NULL;
   float *orog=NULL;
   float ltrans=0.0;
   float dist=0.0, disi=0.0, ndist=0.0;
   float dist_l=0.0, disi_l=0.0, ndist_l=0.0;
   float *p1=NULL, *p2=NULL;
   float fmxmn=0.0;
   float xt, yt;
   float thr_cnt=0.0, thrtmp=0.0;
   float var=0.0;
   float mval=0.0;
   float per2=period*0.5;
   float md=0.0, mdist=0.0, val=0.0;
   float aarea=0.0;
   float *fmask=NULL;

   float scl=1.0;
   float mpress=0.0;
   float mskval=0.0;
   float xll=0.0;

   float rrad=0.0, ravg=0.0;
   float *rr=NULL, *tht=NULL;
   float *scf=NULL, *slat=NULL, *slng=NULL;
   float *xxs=NULL, *yys=NULL;
   
   float circ4[4]={0.0, 90.0, 180.0, 270.0};
   float circ8[8]={0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0};
   
   double sm=0.0;
   double xx=0.0, yy=0.0, zz=0.0, dd[4];
   double s1=0.0, s2=0.0, c1=0.0, c2=0.0;
   double sn=0.0, cn=0.0;
   double norm=0.0, arot=0.0;
   double nearpoleS=0.0, nearpoleN=0.0;
   double rad_off=0.0, cosrad=0.0, sinrad=0.0;
   double rad_off_l=0.0, cosrad_l=0.0, sinrad_l=0.0;
   double cosang[8], sinang[8];
   double ang_off0=0.0, rad_off0, cosang0, sinang0, cosrad0=0.0, sinrad0=0.0;
   double dpN=0.0, dpS=0.0;

   VEC *vv=NULL, *vt=NULL;
   VEC rot[3]={{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
   VEC newv={0.0, 0.0, 0.0};
   VEC npole={0.0, 0.0, 1.0}, spole={0.0, 0.0, -1.0};
   VEC vsdir={0.0, 0.0, 0.0};
   VEC pt1, pt2, tvec;
   VEC vt1, vt2, vt3, vn;
   VEC ptt;
   
   GRID *gtmp=NULL;

   FILE **adf=NULL;
   FILE *tsf=NULL;
   FILE *fsdloc=NULL;
   FILE *flevd=NULL;
   FILE *forog=NULL;
   
   char **nfil=NULL;
   char trout[MAXCHR];
   char mskout[MAXCHR];
   char stdloc[MAXCHR];
   char line[MAXCHR];
   char fornm[MAXCHR];

   struct tot_tr *altr=NULL; 
   struct fet_pt_tr *atr=NULL;

   struct extent ext={0, 0, 0, 0}, ext2={0, 0, 0, 0};
   struct feature_pts fpt={{0}, {0}, {0}, {0}, {0.0, 0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, {0.0, 0.0}, 0, 0, 0, 0};
   struct feature_pts fptt={{0}, {0},{0}, {0}, {0.0, 0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, {0.0, 0.0}, 0, 0, 0, 0};
   struct feature_pts *fpt_loc=NULL;
   struct features fet={NULL, 1};
   struct object obj={NULL, NULL, NULL, NULL, 0, 0, 0, 0, '\0'}, *objt=NULL;
   struct frame_objs fo={NULL, 1, 0, 0, 0, 0};
      
   struct sp_dat *cokn=NULL;        /* spline data for any of the methods   */
   struct rspline *rs=NULL;         /* data structure for use with smoopy   */
   struct sspline *ss=NULL;         /* data structure for use with sphery   */
   struct savedat *sd=NULL;         /* data structure for use with smoopy   */
   struct savedat_sphy *ssd=NULL;   /* data structure for use with sphery   */

/* additional sline information for re-mapping longitudes */

   struct sp_dat *coknn=NULL, *coknt=NULL;
   
   GRDP g1={0, 0, 0.0, 0.0, 0.0, 0.0};
   GRDP g2={0, 0, 0.0, 0.0, 0.0, 0.0};

   tf = 4;
   if(tom == 'g') geo_init = 1;
   
   if(gr->prgr || gr->prty) { 
     gtmp = gr; gr = gr1; iproj = 1;
     g1.prgr = 0;
     g1.prty = 0;
     
     g2.prgr = gr->prgr;
     g2.prty = gr->prty;
     if(gr->prgr){
        g2.alat = gr1->alat;
	g2.alng = gr1->alng;
	g2.sp1 = gr1->sp1;
	g2.sp2 = gr1->sp2;
     }
   }

   nofld = nfld;
   nff = nf;

   fet.fpt = &fpt;
   obj.fet = &fet;
   obj.ext = &ext;
   fo.objs = &obj;

   objt = fo.objs;

   pb = iper;
   x1u = 1;
   x2u = gr->ix;
   y1u = 1;
   y2u = gr->iy;

   dim = gr->ix * gr->iy;
   
   
/* determine first and final frame id's from tracks */

   for(i=0; i < trnum; i++){
       altr = trs + i;
       if(altr->num){
         if(altr->trpt->fr_id < imnfr) imnfr = altr->trpt->fr_id;
         if((altr->trpt + altr->num - 1)->fr_id > imxfr) imxfr = (altr->trpt + altr->num - 1)->fr_id; 
       }
   }

   printf("****INFORMATION****, lowest frame Id is %d and highest frame Id is %d\n\n", imnfr, imxfr);
   
   strncpy(trout, FILTRS, MAXCHR);
   if(iext) strcpy(strstr(trout, EXTENSION), fext);
   strcat(trout, "_addfld");



   if(iproj){
      printf("***WARNING***, determining alternative field values at track feature points \r\n"
	     "               on a projection.                                             \n\n");
   }

  
   
   printf("***WARNING***, additional fields should be on the same grid as initial data file,\r\n"
          "               and in the same format.                                           \n\n");
   
   printf("How many additional fields are required to be sampled at \r\n"
          "the track feature points.                                \n\n");
	  
   scanf("%d", &nwfld);

   nfo = nf;
   nfldo = nfld;
   
   nfld += nwfld;
   nf += nwfld;
   nnfld = nfld;

   nfwpos = (int *)realloc_n(nfwpos, nf*sizeof(int));
   mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));

/* orographic masking */

   printf("Do you want to orographically mask data, 'y' or 'n'\n\n");
   scanf("\n");
   if((igm=getchar()) == 'y'){

/* read surface orography */
       printf("What is the name of the surface orography file?\n\n");
       scanf("%s", fornm);

       if(form != 4) {
          forog = open_file(fornm, "r");
	  fseeko(forog, pl, FSTART);
       }
       else {
          forog = (FILE *)nc_clone((NETCDF_INFO *)fst, fornm, NC_OPEN_MODE);
          ((NETCDF_INFO *)forog)->iframe = 0;

       }

       orog=(float *)calloc(gr->ix * gr->iy, sizeof(float));
       mem_er((orog == NULL) ? 0: 1, gr->ix * gr->iy * sizeof(float));

       printf("What scaling is required to get orography in to km?\n\n");
       scanf("%e", &scl);

       rf = read_field(orog, NULL, scl, forog, 1.0, 'n', 'n', '0', 'y');

       if(form != 4) close_file(forog, fornm);
       else  netcdf_close((NETCDF_INFO *)forog);

       printf("Will a single mask be applied to all fields or seperate masks for each field?\r\n"
              "Input '0' for single or '1' for seperate.                                    \n\n");
       scanf("%d", &igmn);
       if(igmn < 0 || igmn > 1){
          printf("****WARNING****. option not known, using single mask applied to all fields.\n\n");
          igmn = 0;
       }
 
       nmask = (!igmn) ? 1 : nwfld;

       printf("What value should be used to mask the data array?\n\n");
       scanf("%e", &mskval);

/* allocate memory for masks and create masks */

       strncpy(mskout, MASK, MAXCHR);
       if(iext) strcpy(strstr(mskout, EXTENSION), fext);

       forog = open_file(mskout, "w");
       write_header(gr, forog);

       imask = (int **)calloc(nmask, sizeof(int *));
       mem_er((imask == NULL) ? 0 : 1, nmask*sizeof(int *));
       imy = (int *)calloc(nmask, sizeof(int));
       mem_er((imy == NULL) ? 0 : 1, nmask*sizeof(int));
       fmask = (float *)calloc(dim, sizeof(float));
       mem_er((fmask == NULL) ? 0: 1, dim * sizeof(float));
       for(i=0; i < nmask; i++) {
           printf("What is the pressure for level %d in hPa?\n\n", i+1);
           scanf("%f", &mpress);
           *(imask + i) = (int *)calloc(dim, sizeof(int));
           mem_er((*(imask + i) == NULL) ? 0: 1, dim * sizeof(int));

           *(imy + i) = create_mask(*(imask + i), orog, mpress, dim);
           for(j=0; j < dim; j++) *(fmask + j) = *(*(imask + i) + j);
           fprintf(forog, "FRAME  %3d\n", i+1);
/*           fwrite(*(imask + i), sizeof(int) * dim, 1, forog);*/
           fwrite(fmask, sizeof(float) * dim, 1, forog);
           fprintf(forog, "\n");

       }
       free(fmask);

       fseeko(forog, (off_t)0, FSTART);
       fprintf(forog, "%6d %6d %6d\n", gr->ix, gr->iy, nmask);

       close_file(forog, mskout);

   }

/* convert tracks to Cartesian space */

   if(tom == 'g') convert_track(trs, trnum, 0, 0);

/* choose type of sampling */

   printf("Input '0' for basic interpolation, max/min      \r\n"
          "Input '1' for area averages or variance         \r\n"
          "Input '2' for regional radial or lat-lng values \r\n"
          "Input '3' for sampling at fixed locations.      \n\n");
   scanf("%d", &idiag);
   if(idiag < 0 || idiag > 3){
      printf("****ERROR****, diagnostic type not known, aborting.\n\n");
      exit(1);
   }
   else if(idiag == 1){
      printf("Do you want the area average or the area standard deviation, input '0' for average, '1' for STD or '2' for both.\n\n");
      scanf("%d", &i_av_var);
      if(!i_av_var && !iproj){
          printf("When performing area averaging, do you want cos(lat) weighted count above some threshold\r\n"
                 "instead of the actual area average? '0' for no and '1' for yes.                          \n\n");
          scanf("%d", &ia_cnt);
	  if(ia_cnt < 0 || ia_cnt > 1){
	     printf("****ERROR****, %d is incorrect input for this option.\n\n", ia_cnt);
	     exit(1);
	  }
          if(ia_cnt){
             printf("What is the required threshold?\n\n");
             scanf("%f", &thr_cnt);
          }
      }
      
      if(!ia_cnt){
         printf("Do you want area average and standard deviation relative to a threshold, input '0' for no, '1' for absolute threshold \r\n"
	        "or '2' for threshold relative to sampling centre value.                                                                \n\n");
	 scanf("%d", &ia_thr);
	 if(ia_thr < 0 || ia_thr > 2){
	     printf("****ERROR****, %d is incorrect input for this option.\n\n", ia_thr);
	     exit(1);
	 }
	 if(ia_thr > 0){
	    printf("What is the required threshold, input as a fraction if using relative to sampling centre value?\n\n");
            scanf("%f", &thrtmp);
	    if(ia_thr == 1) thr_cnt = thrtmp;
	 }
	 
	 printf("Do you want area as well as statistics, input '0' for no or '1' for yes.\n\n");
         scanf("%d", &iarea);
         if(iarea < 0 || iarea > 1){
	    printf("****ERROR****, %d is incorrect input for this option.\n\n", iarea);
	    exit(1);
         }
	 
      }

      if(i_av_var == 2){
         nfld += nwfld; nf += nwfld; nnfld = nfld;
         nfwpos = (int *)realloc_n(nfwpos, nf*sizeof(int));
         mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));
         for(i=0; i < 2*nwfld; i++) *(nfwpos + nff + i) = 0;
	 if(iarea){
	    nfld += nwfld; nf += nwfld; nnfld = nfld;
            nfwpos = (int *)realloc_n(nfwpos, nf*sizeof(int));
            mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));
            for(i=0; i < 3*nwfld; i++) *(nfwpos + nff + i) = 0;
	 }
      }
      else {
         for(i=0; i < nwfld; i++) *(nfwpos + nff + i) = 0;
	 if(iarea){
	    nfld += nwfld; nf += nwfld; nnfld = nfld;
            nfwpos = (int *)realloc_n(nfwpos, nf*sizeof(int));
            mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));
            for(i=0; i < 2*nwfld; i++) *(nfwpos + nff + i) = 0;
	 }
      }
	 	 
   }
   else if(idiag == 2){
      if(tom != 'g'){
         printf("****ERROR****, this option is only available on the sphere, choose geodesic distance measure\r\n"
                "               to make this implicit.                                                       \n\n");
         exit(1);
      }

      printf("Use system direction as the preferred direction when compositing, 'y' or 'n'.\n\n");
      scanf("\n");
      if(getchar() == 'y') {
         idir = 1;
         printf("****INFORMATION****, using direction of system motion to orientate the composite sampling grid.\n\n");
	 
	 printf("Do you want to use the default location for system direction or one of the added fields, input '0' for default or '1' to choose.\n\n");
	 scanf("%d", &idiradd);
	 
	 if(idiradd < 0 || idiradd > 1){
	    printf("****ERROR****, incorrect specifier, %d\n", idiradd);
	    exit(1);
	 }
	 
	 if(idiradd) {
	 
	    printf("****INFORMATION****, additional field for direction should have location and have no missing location values.\n\n");
	 
	    printf("Which additional field is required for direction?\n\n");
	    scanf("%d", &ifdir);
	    

            if(ifdir < 0 || ifdir > nff){
               printf("****ERROR****, additional field Id. does not exist, exiting.\n\n");
               exit(1);
            }

            if(ifdir){
               if(*(nfwpos + ifdir - 1)){
                  ifdirp = 0;
                  for(j=0; j < ifdir - 1; j++){
                     if(*(nfwpos + j)) ifdirp += 3;
                     else ifdirp += 1;
                  }
               }
               else{
                  printf("****ERROR****, selected field does not have positional information, exiting.\n\n");
                  exit(1);      
               }

            }
	    
	    convert_track(trs, trnum, ifdir, ifdirp);
	  
	 }

         printf("Do you want to use the current point to determine system direction or average over several points\r\n"
                "to improve direction smoothness, input '1' for single point or 'n' the number of points to use.  \r\n"
                "This will depend on the track lifetime and minimum lifetime threshold.                           \n\n");
         scanf("%d", &ndsmth);
         if(ndsmth < 1){
            printf("****ERROR****, number of points cannot be zero or negative for direction smoothing, exiting.\n\n");
            exit(1);
         }
      }
      else {
         printf("****INFORMATION****, using composite grid with standard orientation, i.e. no motion dependence.\n\n");
      }


      if(!iproj){
         nearpoleN = cos(FP_PI2 - *(gr->ygrid + gr->iy - 1) * FP_PI);
         nearpoleS =  cos(FP_PI2 - *(gr->ygrid) * FP_PI);
         printf("%f %f\n", nearpoleS, nearpoleN);
      }
      
      nf = nfo;
      nfld = nfldo;
      nfwpos = (int *)realloc_n(nfwpos, nf*sizeof(int));
      mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));

   }
   else if(idiag == 3){
/* read in data for fixed location sampling */
      printf("what is the file with the standard locations?\n\n");
      scanf("%s", stdloc);
      fsdloc = open_file(stdloc, "r");
      fgets(line, MAXCHR, fsdloc);
      sscanf(line, "%d", &nloc);

      fpt_loc=(struct feature_pts *)calloc(nloc, sizeof(struct feature_pts));
      mem_er((fpt_loc == NULL) ? 0 : 1, nloc * sizeof(struct feature_pts));

      for(i=0; i < nloc; i++){
        fgets(line, MAXCHR, fsdloc);
        sscanf(line, "%f %f", &(((fpt_loc + i)->x).xy), &(((fpt_loc + i)->y).xy));
        if(tom == 'g') geo_conv(fpt_loc + i);
      }
      close_file(fsdloc, stdloc);

      printf("Do you want distance from storm center to fixed locations, 'y' or 'n'\n\n");
      scanf("\n");
      if(getchar() == 'y'){
         ilocd=1;
         nfld += nloc * nwfld;
         nf += nloc * nwfld;
      }


      nfld += nloc * nwfld - 1;
      nf += nloc * nwfld - 1;
      nfwpos = (int *)realloc_n(nfwpos, nf*sizeof(int));
      mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));

      for(i=0; i < trnum; i++){
   
          altr = trs + i;
       
          for(j=0; j < altr->num; j++){
       
              atr = altr->trpt + j;
	      atr->add_fld = (float *)realloc_n(atr->add_fld, nfld * sizeof(float));
	      mem_er((atr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));

              for(k=nofld; k < nfld; k++) *(atr->add_fld + k) = ADD_UNDEF;
       
          }
   
       }

       if(ilocd)
          for(i=0; i < 2*nloc; i++) *(nfwpos + nff + i) = 0;
       else
          for(i=0; i < nloc; i++) *(nfwpos + nff + i) = 0;
      
   }

/* For sampling along track choose type of search region location */

   if(idiag < 3 && nff){
      printf("Do you want to use the default field position to center the sampling region or choose \r\n"
             "from the additional fields. Input '0' for default or '1' to choose.                   \r\n"
             "Note for interpolation, max/min this only provides the default value when a value     \r\n"
             "cannot be determined.                                                                 \n\n");
      scanf("%d", &ifcnt);
      if(ifcnt < 0 || ifcnt > 1){
         printf("****ERROR****, option unknown for centering of sampling region, exiting.\n\n");
         exit(1);
      }
      else if(ifcnt){
         printf("If a missing value is encountered should one of the other chosen additional fields be used          \r\n"
                "(recursive backward search), revert to the default or leave as missing. Input '0' for default, '1' for search and '2' for leave as misssing.\n\n");
         scanf("%d", &itpadd);
         if(itpadd < 0 || itpadd > 2){
            printf("****ERROR****, incorrect identifier, exiting.\n\n");
            exit(1);
         }

         if(idiag == 2){
            printf("Do you want to project the addtional field location onto the direction of the storm, '0' for no or '1' for yes.\n\n");
            scanf("%d", &iprojd);
            if(iprojd < 0 || iprojd > 1){
               printf("****ERROR****, incorrect identifier, exiting.\n\n");
               exit(1);
            }
         }
      }
                  
   }

   if(idiag < 2){
   
      for(i=0; i < trnum; i++){
   
          altr = trs + i;
       
          for(j=0; j < altr->num; j++){
       
              atr = altr->trpt + j;
	      atr->add_fld = (float *)realloc_n(atr->add_fld, nfld*sizeof(float));
	      mem_er((atr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));

              for(k=nofld; k < nfld; k++) *(atr->add_fld + k) = ADD_UNDEF;
       
          }
   
      }

      printf("Sample at point offset from actual point, 'y' or 'n'\n\n");
      scanf("\n");
      if(getchar() == 'y') {
         if(tom != 'g'){
            printf("****ERROR****, this option is only available on the sphere, choose geodesic distance measure\r\n"
                   "               to make this implicit.                                                       \n\n");
            exit(1);
         }
         ipoff = 1;
         printf("What is the angle (degrees) relative to North and the radial distance required for offset sampling?\n\n");
         scanf("%lf %lf", &ang_off0, &rad_off0);

         ang_off0 *= FP_PI;
         rad_off0 *= FP_PI;

         sincos(rad_off0, &sinrad0, &cosrad0);
         sincos(ang_off0, &sinang0, &cosang0);

      }

   }

   if(idiag == 0){
      printf("Do you want to remove the regional average, '1' for yes or '0' for no.\n\n");
      scanf("%d", &ianom);
      if(ianom < 0 || ianom > 1){
	 printf("****ERROR****, incorrect option for anomaly field, ianom = %d.\n\n", ianom);
	 exit(1);
      }
      if(ianom && tom == 'g' && !iproj){
         gr->coslat = (float *)calloc(gr->iy, sizeof(float));
         mem_er((gr->coslat == NULL) ? 0 : 1, (gr->iy)*sizeof(float));
         for(i=0; i < gr->iy; i++) *(gr->coslat + i) = cos(*(gr->ygrid +i) * FP_PI);
      }
   }

   else if(idiag == 1){
      printf("Do you want to first remove the regional average computed from a different region before computing regional average, '1' for yes or '0' for no.\n\n");
      scanf("%d", &ianom);
      if(ianom < 0 || ianom > 1){
	 printf("****ERROR****, incorrect option for anomaly field, ianom = %d.\n\n", ianom);
	 exit(1);
      }
      if(tom == 'g' && !iproj){
         gr->coslat = (float *)calloc(gr->iy, sizeof(float));
         mem_er((gr->coslat == NULL) ? 0 : 1, (gr->iy)*sizeof(float));
         for(i=0; i < gr->iy; i++) *(gr->coslat + i) = cos(*(gr->ygrid +i) * FP_PI);
      }
   }
   else if(idiag == 2){
   
        printf("Do you want sampling on a rotated radial grid or rotated latitude-longitude, input '0' or '1'\n\n");
	scanf("%d", &igsmp_typ);
	if(igsmp_typ < 0 || igsmp_typ > 1){
	   printf("****ERROR****, incorect choice of option for sampling grid.\n\n");
	   exit(1);
	}

        if(igsmp_typ){
	
/* setup for rotated latitude-longitude grid sampling */

           printf("****INFORMATION****, this option samples the field in (theta, phi) coordinates    \r\n"
                  "                     centered on each storm center by seting up points centered   \r\n"
                  "                     on the equator (0, 180) for various theta and phi and then   \r\n"
                  "                     rotating to the storm center.                                \n\n");
		  
           printf("What latitude-longitude spacing is required, xll?\n\n");
	   scanf("%f", &xll);
	   printf("How many points in either direction,nll,are required?\n\n");
	   scanf("%d", &nll);
	   
	   if(!(nll % 2)) ++nll;
	   nll2 = nll / 2;
	   
	   nr = ntheta = nll;
	}
	
	else {
	
/* setup for radial sampling */

           printf("****INFORMATION****, this option samples the field in (r(theta), theta(phi))     \r\n"
                  "                     coordinates on a spherical cap centered on each storm center\r\n"
                  "                     by seting up points around the N. pole (theta = 0) for      \r\n"
                  "                     various theta and phi and then rotating to the storm center.\n\n");

           printf("What is the spherical cap radius in degrees?\n\n");
           scanf("%f", &rrad);

           printf("How many radii, nr, and angles, ntheta, are required?\n\n");
           scanf("%d %d", &nr, &ntheta);
	   
	}	

        npcy = nr * ntheta;

        rr = (float *)calloc(nr, sizeof(float));
        mem_er((rr == NULL) ? 0 : 1, nr*sizeof(float));

        tht = (float *)calloc(ntheta, sizeof(float));
        mem_er((tht == NULL) ? 0 : 1, ntheta*sizeof(float));

        scf = (float *)calloc(npcy, sizeof(float));
        mem_er((scf == NULL) ? 0 : 1, npcy*sizeof(float));

        slat = (float *)calloc(nr, sizeof(float));
        mem_er((slat == NULL) ? 0 : 1, nr*sizeof(float));

        slng = (float *)calloc(ntheta, sizeof(float));
        mem_er((slng == NULL) ? 0 : 1, ntheta*sizeof(float));

        vv = (VEC *)calloc(npcy, sizeof(VEC));
        mem_er((vv == NULL) ? 0 : 1, npcy*sizeof(VEC));

        for(i=0; i < npcy; i++) *(scf + i) = ADD_UNDEF;
	
	if(igsmp_typ){	   
	
	   xxs = tht;
	   yys = rr;
	   
	   *(slat + nll2) = 0.0;
	   *(yys + nll2) = FP_PI2;
           for(i=0; i < nll2; i++) {
	       *(slat + nll2 - i - 1) = -(i + 1) * xll;
	       *(yys + nll2 - i - 1) = FP_PI2 - *(slat + nll2 - i - 1) * FP_PI;
	       *(slat + nll2 + i + 1) = (i + 1) * xll;
	       *(yys + nll2 + i + 1) = FP_PI2 - *(slat + nll2 + i + 1) * FP_PI;
	   }
	   
	   *(slng + nll2) = 0.0;
	   *(xxs + nll2) = 0.0;
           for(i=0; i < nll2; i++) {
	       *(slng + nll2 - i - 1) = - (i + 1) * xll;
	       *(xxs + nll2 - i - 1) = FPI2 + *(slng + nll2 - i - 1) * FP_PI;
	       *(slng + nll2 + i + 1) = (i + 1) * xll;
	       *(xxs + nll2 + i + 1) = *(slng + nll2 + i + 1) * FP_PI;
	   }	   
	
	}
	else {

           for(i=0; i < nr; i++) {
               *(rr + i) = rrad * FP_PI * (float) (1 + i) / nr;
               *(slat + i) = (FP_PI2 - *(rr + i)) / FP_PI;
           }

           for(i=0; i < ntheta; i++) {
              *(tht + i) = FPI2 * (float) i / ntheta;
              *(slng + i) = *(tht + i) / FP_PI;
           }
	
	}

        for(i=0; i < nr; i++){
            sincos(*(rr + i), &s1, &c1);
            for(j=0; j < ntheta; j++){
                sincos(*(tht + j), &s2, &c2);
                vt = vv + i * ntheta + j;
                vt->x = s1 * c2;
                vt->y = s1 * s2;
                vt->z = c1;
            }
	
	}

/* open file for writing data along track */
        fldpl = (off_t **)calloc(trnum, sizeof(off_t *));
        mem_er((fldpl == NULL) ? 0 : 1, trnum*sizeof(off_t *));

        strcat(trout, "_reg");
        tsf = open_file(trout, "w");
        fprintf(tsf, "%6d %10ld %6d %6d %2d %2d %3d %2d %2d %2d %2d %2d\n", trnum, ptnum, ntheta, nr, nwfld, idir, ndsmth, ifcnt, itpadd, iprojd, ifdir, igsmp_typ);

        fwrite(slng, ntheta * sizeof(float), 1, tsf);
        fprintf(tsf, "\n");
        fwrite(slat, nr * sizeof(float), 1, tsf);
        fprintf(tsf, "\n");

        free(slat); free(slng);

        ns = 0;
	
        for(i=0; i < trnum; i++){

            altr = trs + i;
	    
	    if(altr->num){

               *(fldpl + i) = (off_t *)calloc(altr->num, sizeof(off_t));
               mem_er((*(fldpl + i) == NULL) ? 0 : 1, (altr->num) * sizeof(off_t));

               for(j=0; j < altr->num; j++){
                  atr = altr->trpt + j;
                  *(*(fldpl + i) + j) = ftello(tsf);
                  place1 = *(*(fldpl + i) + j);
                  for(k=0; k < nwfld; k++){
                     fwrite(scf, npcy * sizeof(float), 1, tsf);
                     fprintf(tsf, "\n");
                     if(!k) place2 = ftello(tsf);
                  }

               }

               ptnum += altr->num;
	       ++ns;
	    
	    }

        }

        fseeko(tsf, (off_t)0, SEEK_SET);
        fprintf(tsf, "%6d %10ld %6d %6d %2d %2d %3d %2d %2d %2d %2d %2d\n", ns, ptnum, ntheta, nr, nwfld, idir, ndsmth, ifcnt, itpadd, iprojd, ifdir, igsmp_typ);

        nfll = place2 - place1;
        
   }
   
   adf = (FILE **)calloc(nwfld, sizeof(FILE *));
   mem_er((adf == NULL) ? 0 : 1, nwfld * sizeof(FILE *));
   
   nfil = (char **)calloc(nwfld, sizeof(char *));
   mem_er((nfil == NULL) ? 0 : 1, nwfld * sizeof(char *));
   
   nfst = (int *)calloc(nwfld, sizeof(int));
   mem_er((nfst == NULL) ? 0 : 1, nwfld * sizeof(int));
   
   ninv = (int *)calloc(nwfld, sizeof(int));
   mem_er((ninv == NULL) ? 0 : 1, nwfld * sizeof(int));
   
   nlst = (int *)calloc(nwfld, sizeof(int));
   mem_er((nlst == NULL) ? 0 : 1, nwfld * sizeof(int));
   
   scale = (float *)calloc(nwfld, sizeof(float));
   mem_er((scale == NULL) ? 0 : 1, nwfld * sizeof(float));

   inm = (int *)calloc(nwfld, sizeof(int));
   mem_er((inm == NULL) ? 0 : 1, nwfld * sizeof(int));

   if(ifcnt){
      ifncnt = (int *)calloc(nwfld, sizeof(int));
      mem_er((ifncnt == NULL) ? 0 : 1, nwfld * sizeof(int));
      ifncntp = (int *)calloc(nwfld, sizeof(int));
      mem_er((ifncntp == NULL) ? 0 : 1, nwfld * sizeof(int));
   }

   printf("Do you want to choose start and end time frames or use the calculated values, '0' or '1',    \r\n"
          "use '0' and -1 for frame interval for referencing longitudes and latitudes for regional data.\n\n");
   scanf("%d", &ifrch);
   if(ifrch < 0 || ifrch > 1){
      printf("****ERROR****, option for start, end times incorrect.\n\n");
      exit(1);
   }

/* open file for writing level id's for compositing */

   if(ifcnt && idiag == 2)
      flevd = open_file("level_ids.dat", "w");


   for(i=0; i < nwfld; i++){
   
       nfil[i] = (char *)calloc(MAXCHR, sizeof(char));
       mem_er((nfil[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
   
       printf("What is the file containing field ID %d\n\n", i+1);
       scanf("%s", nfil[i]);
       if(!ifrch){
          printf("What are the start, interval and final frame numbers for field ID %d?\r\n", i+1);      
          scanf("%d %d %d", nfst+i, ninv+i, nlst+i);
       }
       else{
          *(nfst+i) = imnfr;
          *(ninv+i) = 1;
          *(nlst+i) = imxfr;
       }

       if( ninv[i] < 1){
          printf("****WARNING****, sampling from the same frame for all track time steps.\n\n");
          if(!app){
             app=(float *)calloc(gr->ix * gr->iy, sizeof(float));
             mem_er((app == NULL) ? 0: 1, gr->ix * gr->iy * sizeof(float));
          }
       }
       
       printf("What scaling is required for field ID %d\n\n", i+1);
       scanf("%f", scale + i);

       if(scale[i] <= 0.0){
          printf("****WARNING****, scaling is negative, field will be inverted (scaled by -1)\n\n");
       }

       if(!idiag){

          printf("Do you want values at the feature point or try and find an associated max or min,\r\n"
                 "for this field, input '-1' for Min, '0' for at feature point or '1' for Max.     \n\n");
          scanf("%d", inm+i);
          if(*(inm + i) < -1 || *(inm + i) > 1){
             printf("****WARNING****, identifier not known defaulting to '0' option.\n\n");
             *(inm + i) = 0;
          }
          if(*(inm + i)) inmax = 1;

       }
       if(ifcnt) {
          printf("Which existing additional field do you want to use for positional information, use '0' for default.\n\n");
          scanf("%d", ifncnt + i);
          if(*(ifncnt + i) < 0 || *(ifncnt + i) > nff){
             printf("****WARNING****, additional field Id. does not exist, using default.\n\n");
             *(ifncnt + i) = 0;
          }

          if(idiag == 2) fprintf(flevd, "%3d\n", *(ifncnt + i));

          if(*(ifncnt + i)){
             if(*(nfwpos + *(ifncnt + i) - 1)){
                *(ifncntp + i) = 0;
                for(j=0; j < *(ifncnt + i) - 1; j++){
                   if(*(nfwpos + j)) *(ifncntp + i) += 3;
                   else *(ifncntp + i) += 1;
                }
             }
             else{
                printf("****WARNING****, selected field does not have positional information, setting to default.\n\n");
                *(ifncnt + i) = 0;       
             }

          }

       }

       if(form != 4) {
          adf[i] = open_file(nfil[i], "r");
	  fseeko(adf[i], pl, FSTART);
	  place1 = pl;
       }
       else {
          adf[i] = (FILE *)nc_clone((NETCDF_INFO *)fst, nfil[i], NC_OPEN_MODE);
          ((NETCDF_INFO *)adf[i])->iframe = 0;
          if(!imiss && ((NETCDF_INFO *)adf[i])->imiss) imiss = 1;
          else if(((NETCDF_INFO *)adf[i])->imiss != imiss){
             printf("****WARNING****, there are fields with and without missing values, do you want to continue, 'y' or 'n'.\n\n");
             scanf("\n");
             if(getchar() == 'n') exit(1);
          }
       }
   
   }

   if(flevd) close_file(flevd, "level_ids.dat");

   if(inmax){
      printf("When finding maxima or minima do you want to use the spline interpolation \r\n"
             "and steepest ascent/descent or a general search, input '0' for spline and \r\n"
             " steepest ascent or '1' for general search.                               \n\n");
      scanf("%d", &mm_type);
      
      if(mm_type){
         printf("Do want to search for a true max or min, '0' for no, '1' for yes.\n\n");
	 scanf("%d", &itmxmn);
	 if(itmxmn < 0 || itmxmn > 1){
	    printf("****ERROR****, incorrect option for true max or min.\n\n");
	    exit(1);
	 }
      }
   }

   printf("Are there missing values? 'y' or 'n'\n\n");
   scanf("\n");
   if(getchar() == 'y') imiss = 1;

   if(imiss){

     printf("What is the missing data value to test against?\n\n");
     scanf("%f", &mval);

     printf("How do you want to compare with missing value?\r\n"
            "'0' -- equality.                              \r\n"
            "'1' -- less than.                             \r\n"
            "'2' -- greater than.                          \n\n");
     scanf("%d", &icmp);

   }

   printf("Do you want to return additional fields to their unscaled values, 'y' or 'n'.\n\n");
   scanf("\n");
   uscl = getchar();

   printf("Do you want time in terms of frame Id's or actual time, '0' or '1'\n\n");
   scanf("%d", &ittype);

   if(ittype){
      printf("****NOTE****, this option assumes the track frame Id's have already \r\n"
             "              been converted to actual time.                        \n\n");
      printf("What is the starting time for the first field, format YYYYMMDDHH?   \n\n");
      scanf("%ld", &stime);

      printf("What is the time step, taking account of the sampling interval chosen for the fields.\n\n");
      scanf("%d", &tstep);
   }

   printf("If tracks are translated in longitude, convert back to actual longitude, \r\n"
          "or if grid has been translated translate tracks and fixed sampling       \r\n"
	  "locations accordingly? 'y' or 'n'                                        \n\n");
   scanf("\n");
   if(getchar() == 'y'){

      if(fabs(cset) > 0.0){
         itrtr = 1;
         ltrans = -cset;
      }
      else{
        printf("Translate tracks back after sampling, 'y' or 'n'\n\n");
        scanf("\n");
        if(getchar() == 'y') {
           ltrans = gr->f_off;
           itrtr = 1;
        }
        else{
           printf("Input the amount to translate by.\n\n");
           scanf("%f", &ltrans);
        }
      }

      for(i=0; i < trnum; i++){
          altr = trs + i;
          for(j=0; j < altr->num; j++){
              atr = altr->trpt + j;
              atr->xf += ltrans;
              if(atr->xf < 0) atr->xf += period;
              else if(atr->xf > period) atr->xf -= period;

              nfc = 0;
              for(k=0; k < nfo; k++){
                 if(*(nfwpos + k)){
                    *(atr->add_fld + nfc) += ltrans;
                    if(*(atr->add_fld + nfc) < 0) *(atr->add_fld + nfc) += period;
                    else if(*(atr->add_fld + nfc) > period) *(atr->add_fld + nfc) -= period;

                    nfc += 3;
                 }
                 else nfc += 1;
              }

          }

      }
      
      if(fpt_loc && itrtr){
        for(i=0; i < nloc; i++){
           ((fpt_loc + i)->x).xy += ltrans;;
           if(tom == 'g') geo_conv(fpt_loc + i);
        }
      }

   }

   if(inmax || idiag == 1 || ianom){

      printf("What is the maximum distance a new Max or Min can be from the feature point \r\n"
             "or the region to average over in the correct units, input an outer and inner\r\n"
             "value?\n\n");
      scanf("%f %f", &dist, &disi);
      if(dist < disi) {
         printf("****ERROR****, outer distance must be greater than inner distance.\r\n"
                "               exiting.                                           \n\n");
         exit(1);
      }
      ndist = disi;
      
      if(ianom){
         printf("Do you want to compute the anomaly field over a larger region than the search region, '0' for no or '1' for yes\n\n");
	 scanf("%d", &ianom_l);
	 if(ianom_l < 0 || ianom_l > 1){
	    printf("****ERROR****, incorrect option for anomaly field.\n\n");
	    exit(1);
	 }
	 if(ianom_l){
	    printf("Input an outer and inner region, inner different from search region.\n\n");
	    scanf("%f %f", &dist_l, &disi_l);
            if(dist_l < disi_l) {
               printf("****ERROR****, outer distance must be greater than inner distance.\r\n"
                      "               exiting.                                           \n\n");
               exit(1);
            }
	    if(disi_l <= disi) printf("****WARNING****, region for computing anomaly is equal to or smaller than sampling region.\n\n");
	    
	    if(tom =='g'){
	      if(!iproj) nl_l = nln(gr, dist_l); 
	      ndist_l = disi_l * FP_PI;
	      
	      if(iproj){
	        rad_off_l = dist_l * FP_PI;
	        sincos(rad_off_l, &sinrad_l, &cosrad_l);
	      }
	    }	    
            else ndist_l = disi_l;
	 }
      }

      if(tom == 'g') {
         if(!iproj) nl = nln(gr, dist); 
	 ndist = disi * FP_PI;
	 
	 
	 if(iproj){
	 
	    printf("Data is on a projection, use 4 or 8 circle samples for determining region extent, \r\n"
	           "input '0' for 4 or '1' for 8.                                                     \n\n");
            scanf("%d", &icirc);
	    if(icirc < 0 || icirc > 1){
	       printf("****ERROR****, identifier %d not know for circle sampling, defaulting to 8.\n\n", icirc);
	       icirc = 1;
	    }
	 
	    rad_off = dist * FP_PI;
	    sincos(rad_off, &sinrad, &cosrad);
	    
	    if(!icirc){
	       for(i=0; i < 4; i++) {
	           circ8[i] = circ4[i];
		   circ8[i] *= FP_PI;
	           sincos(circ8[i], sinang + i, cosang + i);
	       }
	    }
	    
	    else {
	       for(i=0; i < 8; i++){
	           circ8[i] *= FP_PI;
	           sincos(circ8[i], sinang + i, cosang + i);
	       }
	    
	    }
	 
	 }
      }

   }

   if(inmax){

      printf("Output new position of Min or Max, '1' for yes, '0' for no.\n\n");
      scanf("%d", &wnpos);

      if(wnpos < 0 || wnpos > 1){
         printf("****ERROR****, identifier for new position not recognised, defaulting to \r\n"
                "               current position, no new Max or Min.                      \n\n");
         wnpos = 0;
      }

      if(wnpos){

         nnfld = nfld;

         for(i=0; i < nwfld; i++){
             if(*(inm + i)) {
                nnfld += 2;
                *(nfwpos + nff + i) = 1;
             }
             else *(nfwpos + nff + i) = 0;
         }

         for(i=0; i < trnum; i++){
   
            altr = trs + i;
       
            for(j=0; j < altr->num; j++){
       
                atr = altr->trpt + j;
	        atr->add_fld = (float *)realloc_n(atr->add_fld, nnfld * sizeof(float));
	        mem_er((atr->add_fld == NULL) ? 0 : 1, nnfld * sizeof(float));

                for(k=nofld; k < nnfld; k++) *(atr->add_fld + k) = ADD_UNDEF;
       
            }
   
         }

         printf("Do you want last identified new position to be used as the \r\n"
                "starting point for the determination of the next new point.\r\n"
                "Input '1' for yes and '0' for no.                          \n\n");
         scanf("%d", &liter);
         if(liter < 0 || liter > 1){
            printf("****ERROR****, identifier for level interation not recognised, \r\n"
                   "               defaulting to no interation.                    \n\n");
            liter = 0;
         }   

         if(nff) {
            printf("Input additional field ID to use for starting field for last position.\r\n"
                   "Input '0' for default field use.                                      \n\n");
            scanf("%d", &nffld);
            if(nffld < 0 || nffld > nff){
               printf("****ERROR****, field ID not known, using default.\n\n");
               nffld = 0;
            }
            if(nffld && !(*(nfwpos + nffld - 1))){
               printf("****ERROR****, this field ID does not have positional information\r\n"
                      "               reverting to default.                             \n\n");
               nffld = 0;
            }

            nadf = 0;
            for(i=0; i < nffld - 1; i++){
                if(*(nfwpos + i)) nadf += 3;
                else nadf += 1;
            } 

         }      

      }
      else {
	 
	 nnfld = nfld;
         for(i=0; i < nwfld; i++) *(nfwpos + nff + i) = 0;
	    
      }      

   }

   else if(idiag == 3){nnfld = nfld;}
   else {
     nnfld = nfld;
     if(nfwpos && idiag != 2){ 
        for(i=0; i < nwfld; i++) *(nfwpos + nff + i) = 0;
     }
   }
   
   ap=(float *)calloc(gr->ix * gr->iy, sizeof(float));
   mem_er((ap == NULL) ? 0: 1, gr->ix * gr->iy * sizeof(float)); 


   if((idiag == 0 || !inmax) || idiag >= 2){
      printf("What form of interpolation is required,                       \r\n"
             "Input '0' for bi-cubic spline, no missing values.             \r\n"
             "Input '1' for nearest neigbour interpolation, missing values, \r\n"
             "Input '2' for bi-linear interpolation, missing values.        \r\n"
             "Input '3' for bi-cubic interpolation, missing values.         \n\n");
      scanf("%d", &intrp_type);
      if(intrp_type < 0 || intrp_type > 3){
         printf("****ERROR****, interpolation type unknown, exiting.\n\n");
         exit(1);
      }
      printf("****INFORMATION****, using interpolation option %d\n", intrp_type);
   }

   if(!intrp_type){
   
      smty = query_interp(0);
      
      if(iproj && smty){
         printf("****ERROR****, can't use SPHERY for data on a projection.\n\n");
	 exit(1);
      }

      cokn = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
      mem_er((cokn == NULL) ? 0 : 1, sizeof(struct sp_dat));

      if(smty){

         ss = (struct sspline * )malloc_initl(sizeof(struct sspline));
         mem_er((ss == NULL) ? 0 : 1, sizeof(struct sspline));

         ssd = (struct savedat_sphy * )malloc_initl(sizeof(struct savedat_sphy));
         mem_er((ssd == NULL) ? 0 : 1, sizeof(struct savedat_sphy));

       }

       else{

         rs = (struct rspline * )malloc_initl(sizeof(struct rspline));
         mem_er((rs == NULL) ? 0 : 1, sizeof(struct rspline));

         sd = (struct savedat * )malloc_initl(sizeof(struct savedat));
         mem_er((sd == NULL) ? 0 : 1, sizeof(struct savedat));

      }
   
   }

   nxfl = nofld;

   for(i=0; i < nwfld; i++){


       if(!trnum) continue;

/* setup masking */

       if(igm == 'y'){
          if(igmn) {imyy = *(imy + i); imk = *(imask + i);}
          else {imyy = *imy; imk = *imask;}
       }

       if(wnpos && liter){
          ifnd = 0;
          for(k=i-1; k >= 0; k--) {
              if(*(nfwpos + k)) {ifnd = 1; break;}
          }

          if(!ifnd){
             printf("****WARNING****, no previous additional field positional information can be \r\n"
                    "                 found for this additional field, reverting to default.     \n\n");
          }

       }
       
       frc = 0;
       rf = 0;
       iintp = 1;

       if(*(inm + i) < 0 && !mm_type) scale[i] *= (float)*(inm + i);

       ttime = stime;

       isgn = 1;
       if(scale[i] < 0.0){
          isgn = -1;
          scale[i] = fabs(scale[i]);
       }

       if(ninv[i] < 1 && !intrp_type){
          inewcor = 0;
          printf("Field %d, is this referencing to new longitude coordinates, e.g. from a projection, 'y' or 'n'\n\n", i+1);
          scanf("\n");
          if(getchar() == 'y') {
             inewcor = 1;
             apn=(float *)calloc(gr->ix * gr->iy, sizeof(float));
             mem_er((apn == NULL) ? 0: 1, gr->ix * gr->iy * sizeof(float));
             coknn = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
             mem_er((coknn == NULL) ? 0 : 1, sizeof(struct sp_dat));          
          }
       }
           
       while(frc <= frnum && !rf){
       
          if(!frc){
	    rf = read_field(ap, NULL, scale[i], adf[i], isgn, 'n', 'n', '0', 'y');
       
            if(form != 4) {place2 = ftello(adf[i]); chrnum = place2 - place1;}
       
            if(nfst[i] > 1){
              if(form != 4) fseeko(adf[i],(nfst[i]-2)*chrnum, ORIGIN);
	      else ((NETCDF_INFO *)adf[i])->iframe = nfst[i] - 1;	  
	      rf = read_field(ap, NULL, scale[i], adf[i], isgn, 'n', 'n', '0', 'y');
            }

            if( ninv[i] < 1) {
               memcpy(app, ap, gr->ix * gr->iy * sizeof(float));
               if(inewcor) {
                  memcpy(apn, ap, gr->ix * gr->iy * sizeof(float));
                  for(j=0; j < gr->ix * gr->iy; j++){ 
                      if(*(apn + j) > per2) *(apn + j) -= period;
                  }
		  
               }

            }
	    
	    frc = nfst[i];
	    fruser = (!ifrch) ? 1 : imnfr;

	  }
	  else {

            if( ninv[i] > 0 ){
               if(form != 4) fseeko(adf[i], (ninv[i]-1)*chrnum, ORIGIN); 
               else ((NETCDF_INFO *)adf[i])->iframe = frc - 1;
	  
	       rf = read_field(ap, NULL, scale[i], adf[i], isgn, 'n', 'n', '0', 'y');

            }

            else {
               memcpy(ap, app, gr->ix * gr->iy * sizeof(float));
               iintp = 0;
            }

	    ++fruser;

            if(ittype) ttime = new_time(ttime, tstep);
	  
	  }
	  
          if(iintp && !intrp_type){ 
             if(smty) surfit(&sm, itpfc, smty, cokn, ss, ssd);
             else surfit(&sm, itpfc, smty, cokn, rs, sd);
             itpfc = 1;

             if(inewcor){
                if(smty) sphery_setup(1);
                else smoopy_setup(1);
                itpfc = 0;
                atmp = ap;
                ap = apn;
                if(smty) surfit(&sm, itpfc, smty, coknn, ss, ssd);
                else surfit(&sm, itpfc, smty, coknn, rs, sd);
                ap = atmp;
                itpfc = 1;
             }
          } 
  
	  for(j=0; j < trnum; j++){
	 
	      altr = trs + j;
	      
	      if(! altr->num ) continue;

              iin = 0;

              if(!ittype){

	         i1 = altr->trpt->fr_id;
	         i2 = (altr->trpt + altr->num - 1)->fr_id;
	      
	         if(fruser >= i1 && fruser <= i2){	      
	      
	            atr = altr->trpt + (int)fruser - i1;
                    pt_id = (int)fruser - i1;
                    iin = 1;

                 }

               }

               else {

                 if(!altr->time){
                    printf("****ERROR****, actual time is not specified for tracks.\n\n");
                    exit(1);
                 }

	         i1 = altr->trpt->time;
	         i2 = (altr->trpt + altr->num - 1)->time; 

	         if(ttime >= i1 && ttime <= i2){	      
	      
                    for(k=0; k < altr->num; k++){
                        atr = altr->trpt + k;
                        if(atr->time == ttime) {iin = 1; pt_id = k; break;}
                    }

                 }  


               }

               if(iin){

                if(ifcnt){
                   if(find_last(atr, &fpt, ifncnt, ifncntp, i, itpadd)) continue;
                }
                else {

                  (fpt.x).xy = atr->xf;
                  (fpt.y).xy = atr->yf;

                }

                if(ipoff) {
                   geo_conv(&fpt);
                   vec_offset(&fpt, cosang0, sinang0, cosrad0, sinrad0);
                } 
		
/*		if(iproj) proj_point(*pp, &((fpt.x).xy), &((fpt.y).xy), 0, 0, gr->prgr, gr->prty, &(gr1->alat), &(gr1->alng)); */
		if(iproj) proj_point(*pp, &((fpt.x).xy), &((fpt.y).xy), &g1, &g2);

		if(!intrp_type){
	 
	           if(cokn->sm_type){
		       yy = (fpt.x).xy;
		       xx = (fpt.y).xy;		 
		       yy *= FP_PI;
		       xx = FP_PI2 - xx * FP_PI;
		       if(xx < 0.) xx = 0.0;
		       else if(xx > FPI) xx = FPI;		 
	            }
		 
	            else {
		       xx = (fpt.x).xy;
		       yy = (fpt.y).xy; 
	            }

/* re-mapping longitudes */

                    if(inewcor){
                      mdist = LARGE;
                      for(k=0; k < gr->ix; k++){
                          if((md = fabs((fpt.x).xy - *(gr->xgrid + k))) < mdist){mdist = md; kk = k;}
                      }
                      mdist = LARGE;
                      for(k=0; k < gr->iy; k++){
                          if((md = fabs((fpt.y).xy - *(gr->ygrid + k))) < mdist){mdist = md; ii = k;}
                      }

                      val = *(ap + ii * gr->ix + kk);

                      icswap = 0;
                      if(val <= LNGOFF || val >= period - LNGOFF){
                         coknt = cokn;
                         cokn = coknn;
                         icswap = 1;
                      } 

                    }

#ifdef  NOUNDERSCORE

                    bisp(&zz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#else

                    bisp_(&zz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#endif

                    if(icswap) cokn = coknt;

                 }

                 else {
		    xx = (fpt.x).xy;
		    yy = (fpt.y).xy;
		    if(intrp_type == 1) zz = nearest_neigbour_intrp(ap, xx, yy, mval, icmp);
                    else if(intrp_type == 2) zz = bilinear_intrp(ap, xx, yy, mval, icmp);
                    else if(intrp_type == 3) zz = bicubic_intrp(ap, xx, yy, mval, icmp);
                 }


/* mask data array */

                 if(igm == 'y' && imyy){
                    if(idiag < 2 || ianom){
                       for(k=0; k < dim; k++){
                           if(*(imk + k)) *(ap + k) = mskval;
                       }
                    }

                    if(mask_data(imk, (fpt.x).xy, (fpt.y).xy, imyy)) zz = ADD_UNDEF;

                 }
		 
                 if(!idiag){

                    if(ianom){
                       if(ifcnt){
                         if(find_last(atr, &fpt, ifncnt, ifncntp, i, itpadd)) continue;
                       }
                       else {
                         (fpt.x).xy = atr->xf;
                         (fpt.y).xy = atr->yf;
                       }

                       if(tom == 'g') {
                          geo_conv(&fpt);
                          if(ipoff) vec_offset(&fpt, cosang0, sinang0, cosrad0, sinrad0);
                       }
		       
		       if(ianom_l){
                          if(tom == 'g' && iproj) projext(objt, sinrad_l, cosrad_l, sinang, cosang, icirc, &g1, &g2);
                          else detext(objt, dist_l, nl_l);
                          ravg = fmaxmin(&fpt, &fptt, &ext, ndist_l, &xt, &yt, &var, &aarea, 0.0, 0.0, mval, 1, 0, 0, 0, 0, imiss, icmp, 0, iproj, &g1, &g2);		       
		       }
		       else {
                          if(tom == 'g' && iproj) projext(objt, sinrad, cosrad, sinang, cosang, icirc, &g1, &g2);
                          else detext(objt, dist, nl);
                          ravg = fmaxmin(&fpt, &fptt, &ext, ndist, &xt, &yt, &var, &aarea, 0.0, 0.0, mval, 1, 0, 0, 0, 0, imiss, icmp, 0, iproj, &g1, &g2);
		       }
                       if(zz < ADD_CHECK) zz -= ravg;
                    }	    

                    if(*(inm + i)){

                       if(wnpos) *(atr->add_fld + nxfl + 2) = zz;
                       else *(atr->add_fld + nxfl) = zz;

                       if(wnpos && liter){

                          if(i > 0){
                             nmnf = nxfl;
                             ii = i - 1;
                             p1 = p2 = NULL;
                             for(k=ii; k >= 0; k--) {
                                 if(*(nfwpos + k)) {
                                    nmnf -= 3; 
                                    p1 = atr->add_fld + nmnf;
                                    p2 = atr->add_fld + nmnf + 1;
                                    if(*p1 > ADD_CHECK) continue;
                                    break;
                                 }
                                 else --nmnf;
                             }

                             if(!p1){ p1 = &(atr->xf); p2 = &(atr->yf);}

                          }
                          else if(nffld){
                             p1 = atr->add_fld + nadf;
                             p2 = atr->add_fld + nadf + 1;
                          }
                          else{
                            p1 = &(atr->xf);
                            p2 = &(atr->yf);
                          }

                          if(*p1 < ADD_CHECK){
                             (fpt.x).xy = *p1;
                             (fpt.y).xy = *p2;
                          }
                          else {
                             (fpt.x).xy = atr->xf;
                             (fpt.y).xy = atr->yf;
                          } 

                       }
                       else {
                           if(ifcnt){
                             if(find_last(atr, &fpt, ifncnt, ifncntp, i, itpadd)) continue;
                           }
                           else {
                             (fpt.x).xy = atr->xf;
                             (fpt.y).xy = atr->yf;
                           }
                       }

                       fpt.str = zz;
                       if(tom == 'g') geo_conv(&fpt);
                       if(ipoff) vec_offset(&fpt, cosang0, sinang0, cosrad0, sinrad0);

                       if(!mm_type){

                          memcpy(&fptt, &fpt, sizeof(struct feature_pts));

                          if(tom == 'g') geo_conv(&fptt);

/* assign object extent */

                          if(tom == 'g' && iproj) projext(objt, sinrad, cosrad, sinang, cosang, icirc, &g1, &g2);
                          else detext(objt, dist, nl);
			  
/*			  if(iproj) proj_point(*pp, &((fpt.x).xy), &((fpt.y).xy), 0, 0, gr->prgr, gr->prty, &(gr1->alat), &(gr1->alng)); */
			  if(iproj) proj_point(*pp, &((fpt.x).xy), &((fpt.y).xy), &g1, &g2);

                          non_lin_opt(&fo, frst, smty);
                          frst = 1;

                          if(pb == 'y'){
                             if(ext.x1 < 1 && (fpt.x).xy - *(gr->xgrid) < TOLPER){
   
                                (fpt.x).xy = *(gr->xgrid + gr->ix - 1);
                                ext.x1 += gr->ix - 1;
                                ext.x2 = gr->ix;
                                non_lin_opt(&fo, frst, smty);

                             }
                             else if(ext.x2 > gr->ix && *(gr->xgrid + gr->ix - 1) - (fpt.x).xy < TOLPER){

                                (fpt.x).xy = *(gr->xgrid);
                                ext.x1 = 1;
                                ext.x2 -= gr->ix - 1;
                                non_lin_opt(&fo, frst, smty);

                             }

                          }

                          if(igm == 'y' && imyy){
                             if(mask_data(imk, (fpt.x).xy, (fpt.y).xy, imyy)){
                                (fpt.x).xy = ADD_UNDEF;
                                (fpt.y).xy = ADD_UNDEF;
                                fpt.str = ADD_UNDEF;
                             }
                          }
			  
/* 			  if(iproj) proj_point(*pp, &((fpt.x).xy), &((fpt.y).xy), gr->prgr, gr->prty, 0, 0, &(gr1->alat), &(gr1->alng)); */
			  if(iproj && fpt.str < ADD_CHECK) proj_point(*pp, &((fpt.x).xy), &((fpt.y).xy), &g2, &g1);
			  
			  if(tom == 'g') geo_conv(&fpt);

                          if(measure(&fpt, &fptt) <= ndist) {
                             if(ianom) {
			        if(ianom_l){
			           if(tom == 'g' && iproj) projext(objt, sinrad_l, cosrad_l, sinang, cosang, icirc, &g1, &g2);
                                   else detext(objt, dist_l, nl_l);
                                   ravg = fmaxmin(&fpt, &fptt, &ext, ndist_l, &xt, &yt, &var, &aarea, 0.0, 0.0, mval, 1, 0, 0, 0, 0, imiss, icmp, 0, iproj, &g1, &g2);				
				}
				else{		       
			           if(tom == 'g' && iproj) projext(objt, sinrad, cosrad, sinang, cosang, icirc, &g1, &g2);
                                   else detext(objt, dist, nl);
                                   ravg = fmaxmin(&fpt, &fptt, &ext, ndist, &xt, &yt, &var, &aarea, 0.0, 0.0, mval, 1, 0, 0, 0, 0, imiss, icmp, 0, iproj, &g1, &g2);
				}
                                if(fpt.str < ADD_CHECK) fpt.str -= ravg;
                             }
                             if(wnpos){
                                *(atr->add_fld + nxfl) = (fpt.x).xy;
                                *(atr->add_fld + nxfl + 1) = (fpt.y).xy;
                                *(atr->add_fld + nxfl + 2) = fpt.str;
                                if(fpt.str < ADD_CHECK){
                                   if(uscl == 'y') *(atr->add_fld + nxfl + 2) /= scale[i];
                                   if(*(inm + i) < 0) *(atr->add_fld + nxfl + 2) *= (float) *(inm + i);
                                }
                             }
                             else {
                                *(atr->add_fld + nxfl) = fpt.str;
                                if(fpt.str < ADD_CHECK){
                                   if(uscl == 'y') *(atr->add_fld + nxfl) /= scale[i];
                                   if(*(inm + i) < 0) *(atr->add_fld + nxfl) *= (float) *(inm + i);
                                }
                             }

                          }
                          else {
                            if(wnpos){
                              if(*(atr->add_fld + nxfl + 2) < ADD_CHECK){
                                 if(uscl == 'y') *(atr->add_fld + nxfl + 2) /= scale[i];
                                 if(*(inm + i) < 0) *(atr->add_fld + nxfl + 2) *= (float) *(inm + i);
                              }
                            }
                            else {
                              if(*(atr->add_fld + nxfl) < ADD_CHECK){
                                 if(uscl == 'y') *(atr->add_fld + nxfl) /= scale[i];
                                 if(*(inm + i) < 0) *(atr->add_fld + nxfl) *= (float) *(inm + i);
                              }
                            }
                          }

                       }

                       else {

                          if(tom == 'g') geo_conv(&fpt);
                          if(ipoff) vec_offset(&fpt, cosang0, sinang0, cosrad0, sinrad0);
/* assign object extent */

                          if(tom == 'g' && iproj) projext(objt, sinrad, cosrad, sinang, cosang, icirc, &g1, &g2);
                          else detext(objt, dist, nl);
			  ext2 = ext;

/* determine max or min. */

                          if(ianom) {
			  
			     if(ianom_l){
			        if(tom == 'g' && iproj) projext(objt, sinrad_l, cosrad_l, sinang, cosang, icirc, &g1, &g2);
                                else detext(objt, dist_l, nl_l);
			        ravg = fmaxmin(&fpt, &fptt, &ext, ndist_l, &xt, &yt, &var, &aarea, 0.0, 0.0, mval, 1, 0, 0, 0, 0, imiss, icmp, 0, iproj, &g1, &g2);
				ext = ext2;
			     }
			     else {
			        ravg = fmaxmin(&fpt, &fptt, &ext, ndist, &xt, &yt, &var, &aarea, 0.0, 0.0, mval, 1, 0, 0, 0, 0, imiss, icmp, 0, iproj, &g1, &g2); 
		             }
			     
			  }             

                          fmxmn = fmaxmin(&fpt, &fptt, &ext, ndist, &xt, &yt, NULL, NULL, 0.0, ravg, mval, idiag, *(inm + i), 0, 0, ianom, imiss, icmp, itmxmn, iproj, &g1, &g2);

                          if(igm == 'y' && imyy){
                             if(mask_data(imk, xt, yt, imyy)){
                                xt = ADD_UNDEF;
                                yt = ADD_UNDEF;
                                fmxmn = ADD_UNDEF;
                             }
                          }

                          if(wnpos){
                             *(atr->add_fld + nxfl) = xt;
                             *(atr->add_fld + nxfl + 1) = yt;
                             *(atr->add_fld + nxfl + 2) = fmxmn;
                             if(fmxmn < ADD_CHECK){
                                if(uscl == 'y') *(atr->add_fld + nxfl + 2) /= scale[i];
                             }
                          }
                          else {
                             *(atr->add_fld + nxfl) = fmxmn;
                             if(fmxmn < ADD_CHECK){
                                if(uscl == 'y') *(atr->add_fld + nxfl) /= scale[i];
                             }
                          }

                       }

                    }

                    else {
                       *(atr->add_fld + nxfl) = zz;
                       if(zz < ADD_CHECK){
                          if(uscl == 'y') *(atr->add_fld + nxfl) /= scale[i];
                       }
                    }

                 }

                 else if(idiag == 1){

                    if(ifcnt){
                        if(find_last(atr, &fpt, ifncnt, ifncntp, i, itpadd)) continue;
                    }
                    else {
                       (fpt.x).xy = atr->xf;
                       (fpt.y).xy = atr->yf;
		       fpt.str = atr->zf;
                    }

                    if(tom == 'g') geo_conv(&fpt);
                    if(ipoff) vec_offset(&fpt, cosang0, sinang0, cosrad0, sinrad0);

                    if(tom == 'g' && iproj) projext(objt, sinrad, cosrad, sinang, cosang, icirc, &g1, &g2);
                    else detext(objt, dist, nl);
		    ext2 = ext;
		    
		    if(ia_thr == 2) thr_cnt = fpt.str * thrtmp;
		    
		    if(ianom){
		       if(ianom_l){
			  if(tom == 'g' && iproj) projext(objt, sinrad_l, cosrad_l, sinang, cosang, icirc, &g1, &g2);
                          else detext(objt, dist_l, nl_l);
			  ravg = fmaxmin(&fpt, &fptt, &ext, ndist_l, &xt, &yt, &var, &aarea, 0.0, 0.0, mval, 1, 0, 0, 0, 0, imiss, icmp, 0, iproj, &g1, &g2);
		          ext = ext2;
		       }
		       else {
			  ravg = fmaxmin(&fpt, &fptt, &ext, ndist, &xt, &yt, &var, &aarea, 0.0, 0.0, mval, 1, 0, 0, 0, 0, imiss, icmp, 0, iproj, &g1, &g2); 
		       }		    
		    }

                    *(atr->add_fld + nxfl) = fmaxmin(&fpt, &fptt, &ext, ndist, &xt, &yt, &var, &aarea, thr_cnt, ravg, mval, idiag, 0, ia_cnt, ia_thr, ianom, imiss, icmp, 0, iproj, &g1, &g2); 
		    if(iarea) *(atr->add_fld + nxfl + 1) = aarea;

                    if(i_av_var == 1) {
		       *(atr->add_fld + nxfl) = var;
		       if(iarea) *(atr->add_fld + nxfl + 1) = aarea;
                    }
                    else if(i_av_var == 2) {
		       *(atr->add_fld + nxfl + 1) = var;
		       if(iarea) *(atr->add_fld + nxfl + 2) = aarea;
		    }

                 }

                 else if(idiag == 2){

                    if(idir){
		    
		       if(ifdir){
		          if(*(atr->add_fld + ifdirp) > ADD_CHECK) continue;  
		          xx = *(atr->add_fld + ifdirp) * FP_PI;
                          yy = FP_PI2 - *(atr->add_fld + ifdirp + 1) * FP_PI;
		       
		       }
		       else {
                          xx = atr->xf * FP_PI;
                          yy = FP_PI2 - atr->yf * FP_PI;		       
		       }

/* rotate standard vector (-1, 0, 0) */

                       if(yy < 0.) yy = 0.0;
                       sincos(xx, &s1, &c1);
                       sincos(yy, &s2, &c2);
                       vsdir.x = -c1 * c2;
                       vsdir.y = -s1 * c2;
                       vsdir.z = s2;

/* find rotation angle */

                       pt1.x = atr->pp[0];
                       pt1.y = atr->pp[1];
                       pt1.z = atr->pp[2];
                       arot =   dirangle(altr, &cn, &sn, pt_id, ndsmth, ifdir, ifdirp);
		       if(arot > ADD_CHECK) continue; 

                       norm = dotp(&vsdir, &pt1);
                       crosp(&vsdir, &pt1, &vt1);
                       mulv(&vt1, &vt1, sn);
                       mulv(&pt1, &vt2, ((1.0 - cn) * norm));
                       mulv(&vsdir, &vt3, cn);
                       addv(&vt1, &vt2, &vn);
                       addv(&vn, &vt3, &tvec);

                    }

                    if(ifcnt){
                       if(find_last(atr, &fpt, ifncnt, ifncntp, i, itpadd)) continue;
                       xx = (fpt.x).xy * FP_PI;
                       yy = FP_PI2 - (fpt.y).xy * FP_PI;
                    }
                    else {
                       xx = atr->xf * FP_PI;
                       yy = FP_PI2 - atr->yf * FP_PI;
                    }
                    if(yy < 0.) yy = 0.0;

/* project point onto direction */

                    if(iprojd){
                       sincos(xx, &s1, &c1);
                       sincos(yy, &s2, &c2);

                       ptt.x = c1 * s2;
                       ptt.y = s1 * s2;
                       ptt.z = c2;

                       if(idir) fvec(&pt1, &tvec, &ptt, &vn);
                       else {
                           if(pt_id < altr->num -1){
                              pt1.x = atr->pp[0];
                              pt1.y = atr->pp[1];
                              pt1.z = atr->pp[2];
                              pt2.x = (atr+1)->pp[0];
                              pt2.y = (atr+1)->pp[1];
                              pt2.z = (atr+1)->pp[2];
                           }
                           else {
                              pt1.x = (atr-1)->pp[0];
                              pt1.y = (atr-1)->pp[1];
                              pt1.z = (atr-1)->pp[2];
                              pt2.x = atr->pp[0];
                              pt2.y = atr->pp[1];
                              pt2.z = atr->pp[2];
                           }

                           fvec(&pt1, &pt2, &ptt, &vn);
                       }
                       yy = acos(vn.z);
                       xx = atan2(vn.y, vn.x);
                       if(xx < 0.0) xx = FPI2 + xx;
                    }

                    sincos(xx, &s1, &c1);
                    sincos(yy, &s2, &c2);

                    pt1.x = c1 * s2;
                    pt1.y = s1 * s2;
                    pt1.z = c2;
		    
                    if(igsmp_typ){
		       yy = -(FP_PI2 - yy);
                       sincos(yy, &s2, &c2);		    
		    }

                    rot[0].x = c1 * c2;
                    rot[0].y = -s1;
                    rot[0].z = c1 * s2;
                    rot[1].x = s1 * c2;
                    rot[1].y = c1;
                    rot[1].z = s1 * s2;
                    rot[2].x = -s2;
                    rot[2].z = c2;
      
                    for(k=0; k < npcy; k++){

                        newv.x = dotp(&rot[0], &vv[k]);
                        newv.y = dotp(&rot[1], &vv[k]);
                        newv.z = dotp(&rot[2], &vv[k]);

                        if(idir){
                           norm = dotp(&newv, &pt1);
                           crosp(&newv, &pt1, &vt1);
                           mulv(&vt1, &vt1, sn);
                           mulv(&pt1, &vt2, ((1.0 - cn) * norm));
                           mulv(&newv, &vt3, cn);
                           addv(&vt1, &vt2, &vn);
                           addv(&vn, &vt3, &newv);
                        }

                        if(!iproj){
                           dpN = dotp(&newv, &npole);
                           dpS = dotp(&newv, &spole);
                           if(dpN > nearpoleN || dpS > -nearpoleS){			   
/*			      printf("%f %f %f %f\n", dpS, dpN, nearpoleS, nearpoleN); */
                              *(scf + k) = ADD_UNDEF;
                              continue;
                           }
			}

                        if(!intrp_type){
                        
	                   if(cokn->sm_type){
		             yy = atan2(newv.y, newv.x);
                             if(yy < 0.0) yy = FPI2 + yy;
		             xx = acos(newv.z);
                             if(xx < 0.0) xx = 0.0;
                             else if(xx > FPI) xx = FPI;		 
 
	                   }
		 
	                   else {
			   
                             xx = atan2(newv.y, newv.x) / FP_PI;
                             if(xx < 0.0) xx = 360.0 + xx;
                             yy = (FP_PI2 - acos(newv.z)) / FP_PI;
                             if(yy > 90.0) yy = 90.0;
                             else if(yy < -90.0) yy = -90.0;
			     
                             if(iproj) {
			        xt = xx;
				yt = yy;
/*			        proj_point(*pp, &xt, &yt, 0, 0, gr->prgr, gr->prty, &(gr1->alat), &(gr1->alng)); */
				proj_point(*pp, &xt, &yt, &g1, &g2);
				xx = xt;
				yy = yt;
			     }
			     
			     if((xx - xmn) * (xmx - xx) < 0.0 || (yy - ymn) * (ymx - yy) < 0.0) {*(scf + k) = ADD_UNDEF; continue;}
			     
	                   }


                        
#ifdef  NOUNDERSCORE

                           bisp(&zz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#else

                           bisp_(&zz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#endif

                        }

                        else{
                           xx = atan2(newv.y, newv.x) / FP_PI;
                           if(xx < 0.0) xx = 360.0 + xx;
                           yy = (FP_PI2 - acos(newv.z)) / FP_PI;
                           if(yy > 90.0) yy = 90.0;
                           else if(yy < -90.0) yy = -90.0;
			   if(iproj) {
			      xt = xx;
			      yt = yy;
/*			      proj_point(*pp, &xt, &yt, 0, 0, gr->prgr, gr->prty, &(gr1->alat), &(gr1->alng)); */
			      proj_point(*pp, &xt, &yt, &g1, &g2);
			      xx = xt;
			      yy = yt;
			   }
			   
                           if(intrp_type == 1) zz = nearest_neigbour_intrp(ap, xx, yy, mval, icmp);
			   else if(intrp_type == 2) zz = bilinear_intrp(ap, xx, yy, mval, icmp);
                           else if(intrp_type == 3) zz = bicubic_intrp(ap, xx, yy, mval, icmp);
                        }

/* mask data array */

                        if(igm == 'y' && imyy){

                           xx = atan2(newv.y, newv.x) / FP_PI;
                           if(xx < 0.0) xx = 360.0 + xx;
                           yy = (FP_PI2 - acos(newv.z)) / FP_PI;
                           if(yy > 90.0) yy = 90.0;
                           else if(yy < -90.0) yy = -90.0;
			   if(iproj) {
			      xt = xx;
			      yt = yy;
/*			      proj_point(*pp, &xt, &yt, 0, 0, gr->prgr, gr->prty, &(gr1->alat), &(gr1->alng)); */
			      proj_point(*pp, &xt, &yt, &g1, &g2);
			      xx = xt;
			      yy = yt;
			   }
			   
                           if(mask_data(imk, xx, yy, imyy)) zz = ADD_UNDEF;
                        }  
			
                        *(scf + k) = zz;

                    }

/* Write to file */

                    fseeko(tsf, *(*(fldpl + j) + pt_id) + i * nfll, FSTART);

                    fwrite(scf, npcy * sizeof(float), 1, tsf);

                 } 
                 else if(idiag == 3){

                    for(k=0; k < nloc; k++){
		    
		       if(iproj) {
		          xt = ((fpt_loc + k)->x).xy;
			  yt = ((fpt_loc + k)->y).xy;
/*		          proj_point(*pp, &(((fpt_loc + k)->x).xy), &(((fpt_loc + k)->y).xy), 0, 0, gr->prgr, gr->prty, &(gr1->alat), &(gr1->alng)); */
			  proj_point(*pp, &(((fpt_loc + k)->x).xy), &(((fpt_loc + k)->y).xy), &g1, &g2);
		       }

                       if(!intrp_type){
	                  if(cokn->sm_type){
		             yy = ((fpt_loc + k)->x).xy;
		             xx = ((fpt_loc + k)->y).xy;		 
		             yy *= FP_PI;
		             xx = FP_PI2 - xx * FP_PI;
		             if(xx < 0.) xx = 0.0;
		             else if(xx > FPI) xx = FPI;		 
	                  } 
	                  else {
		             xx = ((fpt_loc + k)->x).xy;
		             yy = ((fpt_loc + k)->y).xy; 
	                  }

#ifdef  NOUNDERSCORE

                           bisp(&zz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#else

                           bisp_(&zz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#endif
                        }
                        else {
                           xx = ((fpt_loc + k)->x).xy;
		           yy = ((fpt_loc + k)->y).xy;
			   if(intrp_type == 1) zz = nearest_neigbour_intrp(ap, xx, yy, mval, icmp);
                           else if(intrp_type == 2) zz = bilinear_intrp(ap, xx, yy, mval, icmp);
                           else if(intrp_type == 3) zz = bicubic_intrp(ap, xx, yy, mval, icmp);
                        }

                        if(igm == 'y' && imyy){
		           xx = ((fpt_loc + k)->x).xy;
		           yy = ((fpt_loc + k)->y).xy;
                           if(mask_data(imk, xx, yy, imyy)) zz = ADD_UNDEF;
                        }
			
			if(iproj){
			   ((fpt_loc + k)->x).xy = xt;
			   ((fpt_loc + k)->y).xy = yt;
			}

                        if(ilocd){
                          (fpt.x).xy = atr->xf;
                          (fpt.y).xy = atr->yf;
                          *(atr->add_fld + nxfl + 2 * k) = zz;

                          if(zz < ADD_CHECK){
                             if(tom == 'g') {
                                geo_conv(&fpt);
                                *(atr->add_fld + nxfl + 2 * k + 1) = measure(fpt_loc + k, &fpt) / FP_PI;
                             }
                             else *(atr->add_fld + nxfl + 2 * k + 1) = measure(fpt_loc + k, &fpt);
                          }

                        }
                        else *(atr->add_fld + nxfl + k) = zz;

                    }

                 }

	      }

	  } 
	  
          if(*(ninv + i) > 0) frc += ninv[i];
          else ++frc;

          if(frc > nlst[i]) break;

       }
       
       if(form != 4) close_file(adf[i], nfil[i]);
       else  netcdf_close((NETCDF_INFO *)adf[i]);

       if(*(inm + i) && wnpos) nxfl += 3;
       else if(idiag == 3) {
          if(ilocd) nxfl += 2 * nloc;
          else nxfl += nloc;
       }
       else {
          if(i_av_var == 2) nxfl += 2;
          else nxfl += 1;
	  if(iarea) nxfl += 1;
       }

       if(inewcor) {free(apn); free(coknn);}

   }

   free(fpt_loc);

   if(idiag == 2) close_file(tsf, trout);

   if(gr->coslat) free(gr->coslat);
   if(idiag == 2){ 
      free(rr);
      free(tht);
      free(scf);
      free(vv);
      for(i=0; i < trnum; i++) free(*(fldpl + i));
      free(fldpl);
   }
   
   if(inmax) non_lin_opt(NULL, -1, smty);

   for(i=0; i < nwfld; i++) free(nfil[i]);
   free(nfil);

   nfld = nnfld;
   
   if(iproj) {
      gr = gtmp;
      proj_group(0, 0);
   }

   if(idiag != 2){
      if(itrtr){
        for(i=0; i < trnum; i++){
          altr = trs + i;
          for(j=0; j < altr->num; j++){
              atr = altr->trpt + j;
              atr->xf -= ltrans;
              if(atr->xf < 0) atr->xf += period;
              else if(atr->xf > period) atr->xf -= period;
	      
	      nfc = 0;
	      for(k=0; k < nf; k++){
	         if(*(nfwpos + k)){
		    *(atr->add_fld + nfc) -= ltrans;
		    if(*(atr->add_fld + nfc) < 0) *(atr->add_fld + nfc) += period;          
	            else if(*(atr->add_fld + nfc) > period) *(atr->add_fld + nfc) -= period;
		    
		    nfc += 3;
		 }
		 else nfc += 1;
	      }
          }
        }  
	    
      }   
      tsf = open_file(trout, "w");
      meantrd(tsf, trs, trnum, 's');   
      close_file(tsf, trout);
   }
   
   free(ap);
   if(app) free(app);

   if(!intrp_type){
      if(smty) surfit(&sm, -1, smty, cokn, ss, ssd);
      else surfit(&sm, -1, smty, cokn, rs, sd);
   
      free(cokn);
      free(ss);
      free(ssd);
      free(rs);
      free(sd);
   }

   free(nl);
   free(ifncntp); free(ifncnt);
   free(nfst); free(ninv); free(nlst);
   
   free(adf);

   return;
  
}

void detext(struct object *obj, float dis, int *nl)
{
   int i=0, ii=0;
   int in1=0, in2=0;
   int ic=0;

   float xx, yy, xy;

/* Y-extent */

   xy = (obj->fet->fpt->y).xy;

   i = gr->iy;
   in1 = 0;
   in2 = gr->iy - 1;
   yy = xy - dis;

   if(yy < *(gr->ygrid)) yy = *(gr->ygrid);

   while(i >= 1) {
      i = (in2 - in1) / 2;
      if(yy > *(gr->ygrid + in1 + i)) in1 = in1 + i;
      else in2 = in1 + i;
   }
   obj->ext->y1 = in1 + 1;

   i = gr->iy;
   in1 = 0;
   in2 = gr->iy - 1;
   yy = xy + dis;

   if(yy > *(gr->ygrid + gr->iy - 1)) yy = *(gr->ygrid + gr->iy - 1);

   while(i >= 1) {
	   i = (in2 - in1) / 2;
	   if(yy > *(gr->ygrid + in1 + i)) in1 = in1 + i;
	   else in2 = in1 + i;
   }
   obj->ext->y2 = in2 + 1;

   /* X-extent */

   if(tom == 'e'){

      xy = (obj->fet->fpt->x).xy;

      i = gr->ix - 1;
      in1 = 0;
      in2 = gr->ix - 1;
      xx = xy - dis;
      ic = 0;
      if(xx < *(gr->xgrid)) {
         if(pb == 'y') {xx += period; ic = 1;}
         else xx = *(gr->xgrid);
      }
      while(i >= 1) {
         i = (in2 - in1) / 2;
         if(xx > *(gr->xgrid + in1 + i)) in1 = in1 + i;
         else in2 = in1 + i;
      }
      if(ic) obj->ext->x1 = in1 - gr->ix + 2;
      else obj->ext->x1 = in1 + 1;

      i = gr->ix - 1;
      in1 = 0;
      in2 = gr->ix - 1;
      xx = xy + dis;
      ic = 0;
      if(xx > *(gr->xgrid + gr->ix - 1)) {
        if(pb == 'y') {xx -= period; ic = 1;}
        else xx = *(gr->xgrid + gr->ix - 1);
      }
      while(i >= 1) {
         i = (in2 - in1) / 2;
         if(xx > *(gr->xgrid + in1 + i)) in1 = in1 + i;
         else in2 = in1 + i;
      }
      if(ic) obj->ext->x2 = gr->ix + in2;
      else obj->ext->x2 = in2 + 1;

   }

   else if (tom == 'g') {
       xx = (obj->fet->fpt->x).xy;
       i = gr->ix - 1;
       in1 = 0;
       in2 = gr->ix - 1;
       while(i >= 1) {
          i = (in2 - in1) / 2;
          if(xx > *(gr->xgrid + in1 + i)) in1 = in1 + i;
          else in2 = in1 + i;
       }

       ii = (*(gr->xgrid + in2) - xx < xx - *(gr->xgrid + in1)) ? in2 : in1;

       ii += 1;

       if(yy > 0.){
          obj->ext->x1 = ii - *(nl + obj->ext->y2 - 1);
          obj->ext->x2 = ii + *(nl + obj->ext->y2 - 1);
       }
       else {
          obj->ext->x1 = ii - *(nl + obj->ext->y1 - 1);
          obj->ext->x2 = ii + *(nl + obj->ext->y1 - 1);
       }

       if(pb != 'y') {
         if(obj->ext->x1 < 1) obj->ext->x1 = 1;
         else if(obj->ext->x1 > gr->ix) obj->ext->x1 = gr->ix;
       }

   }

   else {
      printf("****ERROR****, distance measure descriptor '%c' not known in %s at %d\n\n", tom, __FILE__, __LINE__);
      exit(1);

   }


   return;
}

int *nln(GRID *grr, float dis)
{
    int i, j, ig;
    int ii=0;

    int *nl;

    float x, y, z;
    float x1, y1, z1;
    float xx, yy, xy;
    float sum=0.0;
    float cn=0.0;

    double s1, c1, s2, c2;

    cn = 1.0 / cos(dis * FP_PI);

    nl = (int *)calloc(grr->iy, sizeof(int));
    mem_er((nl == NULL) ? 0 : 1, grr->iy * sizeof(int));

    ig = grr->ix / 2;

    xx = *(grr->xgrid + ig) * FP_PI;

    for(i=0; i < grr->iy; i++){
        yy = *(grr->ygrid + i);
        yy = FP_PI2 - yy * FP_PI;
        if(yy < 0.) yy = 0.;

        sincos(xx, &s1, &c1);
        sincos(yy, &s2, &c2);

        x = s2 * c1;
        y = s2 * s1;
        z = c2;

        for(j=ig; j >= 0; j--){
            xy = *(grr->xgrid + j) * FP_PI;
            sincos(xy, &s1, &c1);
            x1 = s2 * c1;
            y1 = s2 * s1;
            z1 = c2;

            sum = cn * (x * x1 + y * y1 + z * z1) - 1.0;
            if(sum < 0.0) break;
            ii = j;
        }

        *(nl + i) = ig - ii + 1;

    }

    return nl;
}


float fmaxmin(struct feature_pts *fpt, struct feature_pts *fptt, struct extent *ext, float ndist, float *xx, float *yy, float *var, float *aarea, float thr_cnt, 
              float aavg, float mval, int idg, int imn, int ia_cnt, int ia_thr, int ianom, int imiss, int icmp, int itmxmn, int iproj, GRDP *g1, GRDP *g2)
{

     int i, j;
     int nn=0;
     int nmiss=0;

     float fmxmn=0.0, fval=0.0; 

     double area=0.0, avg=0.0; 
     double area2=0.0, avg2=0.0;  

     fmxmn = (imn > 0) ? -LARGE : LARGE;

     if(pb == 'y'){
       if(ext->x1 <= 1 && ext->x2 >= gr->ix) --(ext->x2);
     }

     for(i=ext->y1-1; i < ext->y2; i++){

        (fptt->y).xy = *(gr->ygrid + i);
        if(ext->x1 > 0 && ext->x2 <= gr->ix){
           for(j=ext->x1-1; j < ext->x2; j++){
              nmiss = 0;
              if(imiss){
                 if(missing(*(ap + i * gr->ix + j), mval, icmp)) nmiss = 1;
              }
              (fptt->x).xy = *(gr->xgrid + j);

	      if(iproj){
	      	 (fptt->y).xy = *(gr->ygrid + i);
/*	         proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), gr->prgr, gr->prty, 0, 0, &(gr1->alat), &(gr1->alng)); */
		 proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), g2, g1);
	      }
              if(tom == 'g') geo_conv(fptt);
              if(measure(fpt, fptt) <= ndist && !nmiss){                 
                 ++nn;
                 fval = (ianom) ? *(ap + i * gr->ix + j) - aavg: *(ap + i * gr->ix + j);
                 if(!idg) {
		    if(itmxmn) {
		       if(true_maxmin(mval, aavg, j, i, icmp, imiss, imn, ianom))  
		          fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);     
		    }
		    else fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);
		 }
                 else {
                    if(ia_cnt){
                       if(fval >= thr_cnt) {
                         integ(&avg, &area, fval, i, j, iproj);
                         *(ap + i * gr->ix + j) = 0.0;
                       }
                    }
                    else {
		       if(ia_thr && fval < thr_cnt) continue;
                       integ(&avg, &area, fval, i, j, iproj);
                       integ(&avg2, &area2, fval*fval, i, j, iproj);
                    }
                 }
              }
           }                   
        }

        else if(ext->x1 < 1 && ext->x2 <= gr->ix){
           
           for(j=0; j < ext->x2; j++){
              nmiss = 0;
              if(imiss){
                 if(missing(*(ap + i * gr->ix + j), mval, icmp)) nmiss = 1;
              }            
              (fptt->x).xy = *(gr->xgrid + j);
	      if(iproj){
	         (fptt->y).xy = *(gr->ygrid + i);
/*	         proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), gr->prgr, gr->prty, 0, 0, &(gr1->alat), &(gr1->alng)); */
		 proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), g2, g1);
	      }
              if(tom == 'g') geo_conv(fptt);
              if(measure(fpt, fptt) <= ndist && !nmiss){
                 ++nn;
                 fval = (ianom) ? *(ap + i * gr->ix + j) - aavg: *(ap + i * gr->ix + j);
                 if(!idg) {
		    if(itmxmn) {
		       if(true_maxmin(mval, aavg, j, i, icmp, imiss, imn, ianom))  
		          fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);     
		    }
		    else fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);
		 }
                 else {
                    if(ia_cnt){
                       if(fval >= thr_cnt) {
                          integ(&avg, &area, fval, i, j, iproj);
                          *(ap + i * gr->ix + j) = 0.0;
                       }
                    }
                    else {
		       if(ia_thr && fval < thr_cnt) continue;
                       integ(&avg, &area, fval, i, j, iproj);
                       integ(&avg2, &area2, fval*fval, i, j, iproj);
                    }
                 }
              }
           }

           for(j=gr->ix + ext->x1 - 2; j < gr->ix - 1 ; j++){
              nmiss = 0;
              if(imiss){
                 if(missing(*(ap + i * gr->ix + j), mval, icmp)) nmiss = 1;
              }
              (fptt->x).xy = *(gr->xgrid + j);
	      if(iproj){
	         (fptt->y).xy = *(gr->ygrid + i);
/*	         proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), gr->prgr, gr->prty, 0, 0, &(gr1->alat), &(gr1->alng)); */
		 proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), g2, g1);
	      }	      
              if(tom == 'g') geo_conv(fptt);
              if(measure(fpt, fptt) <= ndist && !nmiss){
                 ++nn;
                 fval = (ianom) ? *(ap + i * gr->ix + j) - aavg: *(ap + i * gr->ix + j);
                 if(!idg) {
		    if(itmxmn) {
		       if(true_maxmin(mval, aavg, j, i, icmp, imiss, imn, ianom))  
		          fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);     
		    }
		    else fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);
		 }
                 else {
                    if(ia_cnt){
                       if(fval >= thr_cnt) {
                          integ(&avg, &area, fval, i, j, iproj);
                          *(ap + i * gr->ix + j) = 0.0;
                       }
                    }
                    else {
		       if(ia_thr && fval < thr_cnt) continue;
                       integ(&avg, &area, fval, i, j, iproj);
                       integ(&avg2, &area2, fval*fval, i, j, iproj);
                    }
                 }
              }
           }
        }

        else if(ext->x1 > 0 && ext->x2 > gr->ix){

           for(j=ext->x1-1; j < gr->ix; j++){
              nmiss = 0;
              if(imiss){
                 if(missing(*(ap + i * gr->ix + j), mval, icmp)) nmiss = 1;
              }
              (fptt->x).xy = *(gr->xgrid + j);
	      if(iproj){
	         (fptt->y).xy = *(gr->ygrid + i);
/*	         proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), gr->prgr, gr->prty, 0, 0, &(gr1->alat), &(gr1->alng)); */
		 proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), g2, g1);
	      }
              if(tom == 'g') geo_conv(fptt);
              if(measure(fpt, fptt) <= ndist && !nmiss){
                 ++nn;
                 fval = (ianom) ? *(ap + i * gr->ix + j) - aavg: *(ap + i * gr->ix + j);
                 if(!idg) {
		    if(itmxmn) {
		       if(true_maxmin(mval, aavg, j, i, icmp, imiss, imn, ianom))  
		          fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);     
		    }
		    else fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);
		 }
                 else {
                    if(ia_cnt){
                       if(fval >= thr_cnt) {
                          integ(&avg, &area, fval, i, j, iproj);
                          *(ap + i * gr->ix + j) = 0.0;
                       }
                    }
                    else {		    		      
		       if(ia_thr && fval < thr_cnt) continue;
                       integ(&avg, &area, fval, i, j, iproj);
                       integ(&avg2, &area2, fval*fval, i, j, iproj);
                    }
                 }
              }
           }

           for(j=1; j < ext->x2 - gr->ix + 1; j++){
              nmiss = 0;
              if(imiss){
                 if(missing(*(ap + i * gr->ix + j), mval, icmp)) nmiss = 1;
              }
              (fptt->x).xy = *(gr->xgrid + j);
	      if(iproj){
	         (fptt->y).xy = *(gr->ygrid + i);
/*	         proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), gr->prgr, gr->prty, 0, 0, &(gr1->alat), &(gr1->alng)); */
		 proj_point(*pp, &((fptt->x).xy), &((fptt->y).xy), g2, g1);
	      }
              if(tom == 'g') geo_conv(fptt);
              if(measure(fpt, fptt) <= ndist && !nmiss){
                 ++nn;
                 fval = (ianom) ? *(ap + i * gr->ix + j) - aavg: *(ap + i * gr->ix + j);
                 if(!idg) {
		    if(itmxmn) {
		       if(true_maxmin(mval, aavg, j, i, icmp, imiss, imn, ianom))  
		          fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);     
		    }
		    else fmxmn = maxmin(fptt, fmxmn, fval, xx, yy, imn);
                 }
                 else {
                    if(ia_cnt){
                       if(fval >= thr_cnt) {
                          integ(&avg, &area, fval, i, j, iproj);
                          *(ap + i * gr->ix + j) = 0.0;
                       }
                    }
                    else {
		       if(ia_thr && fval < thr_cnt) continue;
                       integ(&avg, &area, fval, i, j, iproj);
                       integ(&avg2, &area2, fval*fval, i, j, iproj);
                    }
                 }
              }
           }
        }
     }

     if(!idg && !nn){
        fmxmn = ADD_UNDEF;
        *xx = ADD_UNDEF;
        *yy = ADD_UNDEF;
     }

     if(idg){
        if(area > TOLAREA){
           avg /= area;
           avg2 /= area2;
           *var = sqrt(avg2 - avg * avg);
           if(*var < 0.0) *var = 0.0;
	   *aarea = area;
        }
        else{
           avg = ADD_UNDEF;
           avg2 = ADD_UNDEF;
           *var = ADD_UNDEF;
	   *aarea = ADD_UNDEF;
        }
     }

     if(!idg){
       if(fmxmn > ADD_CHECK || fmxmn < -ADD_CHECK){
         fmxmn = ADD_UNDEF;
         *xx = ADD_UNDEF;
         *yy = ADD_UNDEF; 
       } 
     }

     return (!idg) ? fmxmn : (!ia_cnt) ? (float) avg : (float) area;
}

float maxmin(struct feature_pts *fptt, float fmxmn, float fval, float *xx, float *yy, int imn)
{

     if(imn > 0 && fval > fmxmn) {
        fmxmn = fval;
        *xx = (fptt->x).xy;
        *yy = (fptt->y).xy;
     }
     else if(imn < 0 && fval < fmxmn) {
        fmxmn = fval;
        *xx = (fptt->x).xy;
        *yy = (fptt->y).xy;
     }

     return fmxmn;
}

void integ(double *avg, double *area, float fval, int itt, int jtt, int iproj)
{
   int i1, i2, j1, j2;

   float dx1, dx2, dy1, dy2;
   float dx, dy;
   float aa;

   if(gr->ixfrm || gr->iyfrm) {
      j1 = jtt - 1;
      j2 = jtt + 1;

      if(j1 >= 0) dx1 = *(gr->xgrid + jtt) - *(gr->xgrid + j1);
      else if(j1 < 0 && pb == 'y') dx1 = *(gr->xgrid + jtt) - (*(gr->xgrid + gr->ix - 2) - period);
      else dx1 = *(gr->xgrid + j2) - *(gr->xgrid + jtt);

      if(j2 < gr->ix) dx2 = *(gr->xgrid + j2) - *(gr->xgrid + jtt);
      else if(j2 >= gr->ix && pb == 'y') dx2 =  *(gr->xgrid + gr->ix - 1) - *(gr->xgrid + jtt) + *(gr->xgrid + 1);
      else dx2 = *(gr->xgrid + jtt) - *(gr->xgrid + j1);

      dx = 0.5 * (dx1 + dx2);
      if(tom == 'g' && !iproj) dx *= FP_PI;

      i1 = itt - 1;
      i2 = itt + 1;

      if(i1 >= 0) dy1 = *(gr->ygrid + itt) - *(gr->ygrid + i1);
      else dy1 = *(gr->ygrid + i2) - *(gr->ygrid + itt);

      if(i2 < gr->iy) dy2 = *(gr->ygrid + i2) - *(gr->ygrid + itt);
      else dy2 = *(gr->ygrid + itt) - *(gr->ygrid + i1);

      dy = 0.5 * (dy1 + dy2);
      if(tom == 'g' && !iproj) dy *= FP_PI;

      if(tom == 'g' && !iproj) aa = *(gr->coslat + itt) * dx * dy;
      else aa = dx * dy;
      
   }

   else {
      if(tom == 'g' && !iproj) aa = *(gr->coslat + itt); 
      else aa = 1.0; 
   }

   *area += aa;
   *avg += fval * aa;

   return;

}

int find_last(struct fet_pt_tr *atr, struct feature_pts *fpt, int *ifncnt, int *ifncntp, int ic, int itpadd)
{

   int i, j, ii;
   int nmnf=0;
   int nxf=0;

   if(*(ifncnt + ic)){
     (fpt->x).xy = *(atr->add_fld + *(ifncntp + ic));
     (fpt->y).xy = *(atr->add_fld + *(ifncntp + ic) + 1);
     fpt->str = *(atr->add_fld + *(ifncntp + ic) + 2);

     if(itpadd == 2){
       if((fpt->x).xy > ADD_CHECK) return 1;
       else return 0;
     }

     if((fpt->x).xy > ADD_CHECK){ 

       if(itpadd && ic > 0){
          nmnf = *(ifncntp + ic);
          ii = *(ifncnt + ic) - 2;

          nxf = 0;
          for(i=ii; i >= 0; i--){
              if(*(nfwpos + i)) {
                 nmnf -= 3;
                 nxf = 0;
                 for(j=ic-1; j>=0; j--) {
                    if(i+1 == *(ifncnt + j)) nxf = 1;
                 }

                 if(*(atr->add_fld + nmnf) > ADD_CHECK || !nxf) continue;

                 break;
              }
              else --nmnf;            
          }

          if(nxf){
            (fpt->x).xy = *(atr->add_fld + nmnf);
            (fpt->y).xy = *(atr->add_fld + nmnf + 1);
	    fpt->str = *(atr->add_fld + nmnf + 2);
          }
          else {
            (fpt->x).xy = atr->xf;
            (fpt->y).xy = atr->yf;
	    fpt->str = atr->zf;
          }

        }
        else {
          (fpt->x).xy = atr->xf;
          (fpt->y).xy = atr->yf;
	  fpt->str = atr->zf;
        }
     }
   }
   else {
     (fpt->x).xy = atr->xf;
     (fpt->y).xy = atr->yf;
     fpt->str = atr->zf;
   }

   return 0;
}

double fvec(VEC *a, VEC *b, VEC *c, VEC *r)
{

    VEC an, bn;

    double aa, bb, dotab;
    double norm=0.0;

    dotab = dotp(a, b);
    aa = dotp(a, c) - dotp(b, c) * dotab;
    bb = dotp(b, c) - dotp(a, c) * dotab;

    mulv(a, &an, aa);
    mulv(b, &bn, bb);

    addv(&an, &bn, r);

    norm = sqrt(dotp(r, r));
    normv(r, norm);

    return dotab;

}

double dirangle(struct tot_tr *altr, double *cn, double *sn, int pt_id, int ndsmth, int ifdir, int ifdirp)
{
   int i=0;
   int n2, st=0, en=0;
   int narot=0;

   double aa=0.0, arot=0.0, norm=0.0;
   double dp=0.0;
   double xx, yy;
   double s1, s2, c1, c2;

   struct fet_pt_tr *atr=NULL;

   VEC vsdir={0.0, 0.0, 0.0};
   VEC pt1, pt2, tvec, vtmp;
   VEC vt1, vt2, vt3, vn;

   n2 = ndsmth / 2;

   atr = altr->trpt + pt_id;

   if(altr->num <= 1){
      printf("****ERROR****, insufficient number of points in this track to determine direction.\n\n");
      exit(1);
   }

   if(ndsmth % 2) {st = pt_id - n2; en = pt_id + n2;}
   else {st = pt_id - n2 - 1; en = pt_id + n2;}

   if(st < 0) st = 0;
   if(en > altr->num - 2) en = altr->num - 2;
   if(st > en) st = en;

   for(i = st; i <= en; i++){
       atr = altr->trpt + i;
       if(ifdir){
          xx = *(atr->add_fld + ifdirp) * FP_PI;
          yy = FP_PI2 - *(atr->add_fld + ifdirp + 1) * FP_PI;     
       }
       else {
          xx = atr->xf * FP_PI;
          yy = FP_PI2 - atr->yf * FP_PI;
       }      
       if(yy < 0.) yy = 0.0;
       
       sincos(xx, &s1, &c1);
       sincos(yy, &s2, &c2);
       vsdir.x = -c1 * c2;
       vsdir.y = -s1 * c2;
       vsdir.z = s2;

       pt1.x = atr->pp[0];
       pt1.y = atr->pp[1];
       pt1.z = atr->pp[2];
       pt2.x = (atr+1)->pp[0];
       pt2.y = (atr+1)->pp[1];
       pt2.z = (atr+1)->pp[2];
       crosp(&pt1, &pt2, &vtmp);
       crosp(&vtmp, &pt1, &tvec);
       norm = sqrt(dotp(&tvec, &tvec));
       normv(&tvec, norm);
       
       if(norm <= 0.0) continue;

/* compute rotation angle */
       dp = dotp(&vsdir, &tvec);
       if(dp > 1.0) dp = 1.0;
       else if(dp < -1.0) dp = -1.0;
       aa = acos(dp);

/* compute which side of tvec that vsdir is for anticlockwise rotation */

       norm = dotp(&vtmp, &vsdir);
       if(norm < 0.0) aa *= -1.0;
       sincos(aa, &s1, &c1);
/* rotate standard direction for check */
       norm = dotp(&vsdir, &pt1);
       crosp(&vsdir, &pt1, &vt1);
       mulv(&vt1, &vt1, s1);
       mulv(&pt1, &vt2, ((1.0 - c1) * norm));
       mulv(&vsdir, &vt3, c1);
       addv(&vt1, &vt2, &vn);
       addv(&vn, &vt3, &vn);

       if(fabs(dotp(&vn, &tvec) - 1.0) > 1.0e-4){
         printf("****ERROR****, orientating radial grid to storm direction incorrect, dotp=%e > 1.0e-6.\n\n", fabs(dotp(&vn, &tvec) - 1.0));
         exit(1);
       }

       arot += aa;
       ++narot;

   }
   
   if(narot == 0) return ADD_UNDEF;

   arot /= narot;

   sincos(arot, sn, cn);

   return arot;
}

void vec_offset(struct feature_pts *fpt, double cosang, double sinang, double cosrad, double sinrad)
{

   float xx, yy;
   double s1, c1, s2, c2;
   double norm=0.0;

   VEC vsdir={0.0, 0.0, 0.0};
   VEC pt1, vt1, vt2, vt3, vn;

   xx = (fpt->x).xy * FP_PI;
   yy = FP_PI2 - (fpt->y).xy * FP_PI;
   if(yy < 0.) yy = 0.0;
   sincos(xx, &s1, &c1);
   sincos(yy, &s2, &c2);
   vsdir.x = -c1 * c2;
   vsdir.y = -s1 * c2;
   vsdir.z = s2;

   pt1.x = fpt->gwk[0];
   pt1.y = fpt->gwk[1];
   pt1.z = fpt->gwk[2];
   
/* rotate vector */

   norm = dotp(&vsdir, &pt1);
   crosp(&vsdir, &pt1, &vt1);
   mulv(&vt1, &vt1, sinang);
   mulv(&pt1, &vt2, ((1.0 - cosang) * norm));
   mulv(&vsdir, &vt3, cosang);
   addv(&vt1, &vt2, &vn);
   addv(&vn, &vt3, &vsdir);

   mulv(&pt1, &pt1, cosrad);
   mulv(&vsdir, &vsdir, sinrad);

   addv(&pt1, &vsdir, &vn);

   fpt->gwk[0] = vn.x;
   fpt->gwk[1] = vn.y;
   fpt->gwk[2] = vn.z;

   (fpt->x).xy = atan2(vn.y, vn.x) / FP_PI;
   if((fpt->x).xy < 0.0) (fpt->x).xy = 360.0 + (fpt->x).xy;
   (fpt->y).xy = (FP_PI2 - acos(vn.z)) / FP_PI;
   if((fpt->y).xy > 90.0) (fpt->y).xy = 90.0;
   else if((fpt->y).xy < -90.0) (fpt->y).xy = -90.0;

   return;

}

int true_maxmin(float mval, float ravg, int ix, int iy, int icmp, int imiss, int imn, int ianom)
{

    int i=0 ,j=0, ii=0, jj=0;
    int isx=0, isy=0;
    int nxx=0;
    
    int dcount=0;
    
    float val=0.0, nval=0.0;
    
    nxx = (pb == 'y') ? gr->ix - 1 : gr->ix;

    if(ix < 0) val = *(ap + iy * gr->ix + ix + nxx); 
    else if(ix >= gr->ix) val = *(ap + iy * gr->ix + ix - nxx);
    else val = *(ap + iy * gr->ix + ix);

    if(ianom) val -= ravg;
    
    isx = ix - 1;
    isy = iy - 1;
    
    for(i=0; i < 3; i++){
       ii = isy + i;
       if(ii < 0 || ii >= gr->iy) continue; 
       else{
          for(j=0; j < 3; j++){
	      if(i == 1 && j == 1) continue;
	      jj = isx + j;
              if(jj < 0) nval = *(ap + ii * gr->ix + jj + nxx);
	      else if(jj >= gr->ix) nval = *(ap + ii * gr->ix + jj - nxx);
	      else nval = *(ap + ii * gr->ix + jj);
	      
	      if(ianom) nval -= ravg;
	      
	      if(imiss){
                 if(missing(nval, mval, icmp)) continue;
              }
	      if(imn == 1 && nval < val) ++dcount;
	      else if(imn == -1 && nval > val) ++dcount;
          }
       }
    }
    
    if(dcount == 8) return 1;
    else return 0;

}


/* determine extent on projection of spherical circle for data originally on a projection */

void projext(struct object *obj, double sinrad, double cosrad, double *sinang, double *cosang, int icirc, GRDP *g1, GRDP *g2)
{

   int i=0;
   int ncirc=0;
   int ii=0, ip=0;
   int in1=0, in2=0;
   
   float xx=0.0, yy=0.0;
   float dif1=0.0, dif2=0.0;
   
   struct feature_pts fpt={{0}, {0}, {0}, {0}, {0.0, 0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, {0.0, 0.0}, 0, 0, 0, 0};
   
   ncirc = (!icirc) ? 4 : 8;
   
   obj->ext->x1 = LINT;
   obj->ext->x2 = 0;
   obj->ext->y1 = LINT;
   obj->ext->y2 = 0;
   
   for(i=0; i < ncirc; i++){
       (fpt.x).xy = (obj->fet->fpt->x).xy;
       (fpt.y).xy = (obj->fet->fpt->y).xy;
       fpt.gwk[0] = obj->fet->fpt->gwk[0];
       fpt.gwk[1] = obj->fet->fpt->gwk[1];
       fpt.gwk[2] = obj->fet->fpt->gwk[2];
	 
       vec_offset(&fpt, cosang[i], sinang[i], cosrad, sinrad);
       
       xx = (fpt.x).xy;
       yy = (fpt.y).xy;
       
/*       proj_point(*pp, &xx, &yy, 0, 0, gr->prgr, gr->prty, &(gr1->alat), &(gr1->alng)); */
       proj_point(*pp, &xx, &yy, g1, g2);
       
       ii = gr->ix;
       in1 = 0;
       in2 = gr->ix - 1;

       if(xx < *(gr->xgrid)) xx = *(gr->xgrid);
       else if (xx > *(gr->xgrid + gr->ix - 1)) xx = *(gr->xgrid + gr->ix - 1);

       while(ii >= 1) {
          ii = (in2 - in1) / 2;
          if(xx > *(gr->xgrid + in1 + ii)) in1 = in1 + ii;
          else in2 = in1 + ii;
       }
       
       dif1 = fabs(xx - *(gr->xgrid + in1));
       dif2 = fabs(xx - *(gr->xgrid + in2));
       
       ip = (dif1 < dif2) ? in1 : in2;
       
       if(ip > obj->ext->x2) obj->ext->x2 = ip;
       if(ip < obj->ext->x1) obj->ext->x1 = ip;
       
       ii = gr->iy;
       in1 = 0;
       in2 = gr->iy - 1;

       if(yy < *(gr->ygrid)) yy = *(gr->ygrid);
       else if (yy > *(gr->ygrid + gr->iy - 1)) yy = *(gr->ygrid + gr->iy - 1);

       while(ii >= 1) {
          ii = (in2 - in1) / 2;
          if(yy > *(gr->ygrid + in1 + ii)) in1 = in1 + ii;
          else in2 = in1 + ii;
       }       
       
       dif1 = fabs(xx - *(gr->ygrid + in1));
       dif2 = fabs(xx - *(gr->ygrid + in2));
       
       ip = (dif1 < dif2) ? in1 : in2;
       
       if(ip > obj->ext->y2) obj->ext->y2 = ip;
       if(ip < obj->ext->y1) obj->ext->y1 = ip;
       
   }
   

   return;
   
}


#endif


