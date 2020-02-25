#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>
#include "proj.h"
#include "file_cat_out.h"
#include "st_fo.h"
#include "statistic.h"

#include "file_handle.h"
#include "mem_er.h"
#include "m_values.h"
#include "grid.h"
#include "splice.h"
#include "geo_values.h"
#include "sqt.h"
#include "tele.h"

#define CLKTCK sysconf(_SC_CLK_TCK)

#define  TOLBOUND  0.001    /* tolerance for boundary detection */
#define  FCTSEG    90.0     /* longitudinal anguler segment size */
#define  TTSTEP    0.0001   /* tolerance for lifetimes */
#define  TOLIPER   0.0001   /* tolerance for division by the number of time units */
#define  LPERC     0.333    /* percentage of track allowed outside chosen frame range */
#define  GXMN      0.0      /* maximum ranges for global data */
#define  GXMX      360.0
#define  GYMN      -90.0
#define  GYMX      90.0
#define  R5        2.236067977 /* square root of 5 */
#define  TOLCONIC  0.00001

/* function to compute the propogation speeds and frequencies of the
   trajectories and statistics.                                      */

float sp_kernal_estimate(double ** , struct dpt * , struct dpt * , float * , int , int , int , int * , int , float * , LEAF * , VEC * , int , int , TELE * );
void meantrd(FILE * , struct tot_tr * , int , int );
void plot_stats(struct tot_stat * , int , int , int );
int powi(int , int );
struct tot_tr *phase(struct tot_tr * , int , int * , int , int );
struct tot_stat *read_stats(FILE * );
struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * );
void splice_plot(struct tot_tr * , int , int, int , int , int );
void statdmp(FILE * , struct tot_stat * );

void select_region(float * , float * , float * , float * );
struct pt_stat *assign_sample_pt(struct tot_stat * , float , float , float , float );
struct dpt *extract_data(struct fet_pt_tr *at, struct dpt *, int * , float , float , float , float , int , float * , int );
double arc(struct fet_pt_tr * , struct dpt * );
float *weights(float * , int , float , int );
void available_stats(void);
double area_sp_triangle(SQT * , VEC * , int );
VEC *create_sqt(SQT ** , VEC * , int , int , int *, int );
void data_partition(SQT ** , VEC * , int , int , int , ... );
void sqt_neighbour_find(SQT ** , int , int , int , int );
struct dpt *grid_convert(struct tot_stat * );
void convert_track(struct tot_tr * , int , int , int );
void netcdf_write_stats(struct tot_stat * , GRID * , char * , char * , int , int );
void write_gstat(struct tot_stat * , GRID * , int , int , int );

extern float xmn, ymn, xmx, ymx;
extern GRID *gr, *gr1, *gr2;
extern int x1u, x2u, y1u, y2u;

int trd=0;
extern int tom;
extern int aniso;
extern int iper_num;
extern int gof;

extern float sum_wt, sum_per;

extern char *fext;

extern int iext;
extern char *trfil;

/* additional fields */

extern int nfld, nf;
extern int *nfwpos;

extern int trtyp;

extern int gppr, ippr;

void statistic(struct tot_tr *all_tr, int trackn)

{

   int i, j, k, ip=0, ifmp=0;
   int dcount=0;
   int trphn=0;                   /* total number of phase speed tracks */
   int tot_fet=0;
   int nb, nbx=0, nby=0;
   int f1, f2;
   int fb, fe;
   int ff1, ff2;
   int d1=0, d2=0, dd=0;
   int tpr, gpr;
   int prty;
   int cpty;
   int dtrn=0, dgnn=0, dlyn=0, dphn=0, dtdn=0;
   int dgrn=0, dann=0, dtnn=0, darn=0;
   int ddd;
   int scc=0, ii=0;
   int fnc=0;
   int nmiss=0;
   int nlif=0;
   int awt=0;
   int isc=0, styp=0;
   int imiss=0;
   int igf=0;
   int iadf=0, ipos=0;
   int iff[3]={0, 0, 0};
   int nftmp, nfltmp;
   int npav=0;
   int ift=0;

   int isamp=0;                   /* flag for sampling methodology  */
   int iplat=0, nlev=0, ntrang=0, nfc=0;
   int nvec=0, nelm=0;
   int mfpts=0;
   int ioff=0;
/*   int istmax=0; */
   int igrwth=0, ngrwth=0;
   int n_loc_add=0;

   clock_t tstartall, tstarttd, tendall;
   struct tms usage;

   float xr1, xr2, yr1, yr2;
   float xt1=0.0, yt1=0.0;
   float ta=0.0, bl=0.0, ypos=0.0;
   float *ltm=NULL;
   float ts;
   float arcm=0, arcl;
   float xaa=0., yaa=0.;
   float xat1=0., yat1=0., xat2=0., yat2=0.;
   float xa1=0., xa2=0., ya1=0., ya2=0.;
   float spav=0.;
   float nn=0.0;
   float xinit=0.0, yinit=0.0;
   float sc=1.0;
   float tarea=0.0, rtar=0.0;

   double **den=NULL;
   double sscale=0.0;

   float *wdtr=NULL, *wdgn=NULL, *wdly=NULL, *wdph=NULL, *wdtd=NULL;
   float *wdgr=NULL, *wdan=NULL, *wdar=NULL, *wdtn=NULL;

   float wt, tolwt=0.0;
/*   float strmax=0.0; */

   struct tot_tr *phsp=NULL, *altr=NULL, *altr2=NULL;
   struct tot_tr *tinit=NULL;
   struct fet_pt_tr *at, *at2, att;

   struct tot_stat *trsav=NULL;
   struct pt_stat *pst;

   struct dpt *dtr=NULL, *dgn=NULL, *dly=NULL, *dph=NULL, *dtd=NULL;
   struct dpt *dgr=NULL, *dan=NULL, *dtn=NULL, *dar=NULL;

   struct dpt *stat=NULL;       /* grid data converted to Cartesian (X, Y, Z) coordinates */
   struct dpt *stt=NULL;

   FILE *ph=NULL, *gend=NULL, *lysd=NULL;

   GRID *gstat1=NULL, *gstat2=NULL;

   PRFP prj=NULL;

   VEC oct_vecs[]={XP, YP, ZP, XPM, YPM, ZPM};
   VEC icos_vecs[]={NP, NA, NB, NC, ND, NE, SP, SA, SB, SC, SD, SE};
   VEC *gvecs=NULL;

   SQT **sqth=NULL;
   LEAF *lf=NULL;

   TELE *tele=NULL;

/* definition of polyhedra, 'u' up triangle, 'd' down triangle, numbering of
   verticies is from up or down vertex, anti-clockwize for 'u' triangle,
   clockwize for 'd' triangle.                                                */


   SQT oct[]={
       {'u', NULL, {NULL, NULL, NULL, NULL}, {2, 0, 1}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, {2, 1, 3}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, {2, 3, 4}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, {2, 4, 0}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, {5, 0, 1}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, {5, 1, 3}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, {5, 3, 4}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, {5, 4, 0}}
   };


   SQT icos[]={
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 0,  1,  2}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 0,  2,  3}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 0,  3,  4}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 0,  4,  5}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 0,  5,  1}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, {10,  1,  2}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, {11,  2,  3}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, { 7,  3,  4}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, { 8,  4,  5}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, { 9,  5,  1}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 2, 10, 11}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 3, 11,  7}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 4,  7,  8}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 5,  8,  9}},
       {'u', NULL, {NULL, NULL, NULL, NULL}, { 1,  9, 10}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, { 6, 10, 11}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, { 6, 11,  7}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, { 6,  7,  8}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, { 6,  8,  9}},
       {'d', NULL, {NULL, NULL, NULL, NULL}, { 6,  9, 10}}
   };

   char charin[MAXCHR];
   char phasef[MAXCHR];
   char statsf[MAXCHR];
   char statss[MAXCHR];

   char stcdff[MAXCHR];

   strncpy(phasef, PHTRS, MAXCHR);
   strncpy(statsf, STATTRS, MAXCHR);
   strncpy(statss, STATTRS_SCL, MAXCHR);
   
   if(iext){
      strcpy(strstr(phasef, EXTENSION), fext);
      strcpy(strstr(statsf, EXTENSION), fext);
      strcpy(strstr(statss, EXTENSION), fext);
   }


/* to compute phase speeds data needs to be defined on geographical coordinates 
   e.g. Plate Caree.                                                           */


   printf("*****     SS   TTTTTTTTTTTT      AA      TTTTTTTTTTTT     SS     *****\n"
          "*****   SS  SS      TT          AAAA          TT        SS  SS   *****\n"
          "*****  SS           TT         AA  AA         TT       SS        *****\n"
          "*****   SS          TT        AA    AA        TT        SS       *****\n"
          "*****     SS        TT       AA      AA       TT          SS     *****\n"
          "*****      SS       TT      AAAAAAAAAAAA      TT           SS    *****\n"
          "*****       SS      TT      AA        AA      TT            SS   *****\n"   
          "*****  SS  SS       TT      AA        AA      TT       SS  SS    *****\n"
          "*****    SS         TT      AA        AA      TT         SS      *****\n\n");


   if(trtyp != 's'){
      printf("****ERROR****, incorrect track file type for statistics calculations,\r\n"
             "               should be 's' actually have '%c'                      \n\n", trtyp);
      exit(1);
   }
 
   printf("***WARNING***, to compute phase speeds and statistics on a spherical \r\n"
          "               domain track data should be defined on a geographic   \r\n"
          "               (lat.-long.) coordinate system.                       \n\n");

   proj_report(gr->prty, gr->prgr);


   if(ippr > 0){

      printf("***ERROR***, cannot cope with track data on this projection \r\n"
             "             at the moment re-create data on a cylindrical  \r\n"
             "             (lat.-long.) projection.                       \n\n");
      exit(1);

   } 

   printf("what frame interval is required to compute statistics, f1 to f2\r\n"
          "ignored if displaying previous sets of statistics.              \n");
   scanf("%d %d", &f1, &f2);

   printf("Have the frame Id's been corrected for any missing time steps?   \r\n"
          "Input '0' for no change and '1' for corrected. Input '0' if there\r\n"
          "are no missing time steps in the data.                           \n\n");
   scanf("%d", &imiss);
   if(imiss > 1) imiss = 1;
   else if(imiss < 0) imiss = 0;

   if(nf){

      printf("****INFORMATION****, there are %d additional fields available, do \r\n"
             "                     you want to use any of these instead of the  \r\n"
             "                     default? Input additional field ID, zero is  \r\n"
             "                     default.                                     \n\n", nf);

      scanf("%d", &iadf);
      if(iadf < 0 || iadf > nf) {
        printf("****WARNING****, additional field ID %d not known using default.\n\n", iadf);
        iadf = 0;
      }


      if(iadf){
         iff[2] = 0;
         for(i=0; i < iadf; i++){
             if(*(nfwpos + i)) iff[2] += 3;
             else iff[2] += 1;
         }
         --(iff[2]);

         if(nfwpos[iadf-1]) {
            printf("****INFORMATION****, additional field has positional information,   \r\n"
                   "                     do you want to use this instead of the default.\r\n"
                   "                     Input '1' for yes and '0' for no.              \n\n");
            scanf("%d", &ipos);
            if(ipos < 0 || ipos > 1) {
               printf("****WARNING****, position option %d not known using default.\n\n", ipos);
               ipos = 0;
            }

            printf("If chosen additional value's have missing locational values do you want\r\n"
                   "to exclude them from the statistics, 'y' or 'n'.                       \n\n");
            scanf("\n");
            if(getchar() == 'y') n_loc_add = 1; 
            
            iff[0] = iff[2] - 2;
            iff[1] = iff[2] - 1;

            for(i=0; i < trackn; i++){
                altr = all_tr + i;
                if(altr->trpt != NULL){
                   at = altr->trpt;
                   for(j=0; j < altr->num; j++){
                       at = altr->trpt + j;
                       if(!(at->add_fld[iff[0]] > ADD_CHECK) && ipos){
                          at->xf = at->add_fld[iff[0]];
                          at->yf = at->add_fld[iff[1]];
                       }
                       if(n_loc_add && at->add_fld[iff[0]] > ADD_CHECK) at->zf = ADD_UNDEF;
                       else at->zf = at->add_fld[iff[2]];
                   }

                }

            }

         }
         else {

            for(i=0; i < trackn; i++){
                altr = all_tr + i;
                if(altr->trpt != NULL){
                   at = altr->trpt;
                   for(j=0; j < altr->num; j++){
                       at = altr->trpt + j;
                       at->zf = at->add_fld[iff[2]];
                   }

                }

            }

         }

      }

   }

   printf("do you wish to display phase speeds                   \r\n"
          "Input '0' for no                                      \r\n"
          "      '1' to compute a new set of phase speeds        \r\n"
          "      '2' to display an existing set of phase speeds  \r\n");

   scanf("%d", &ip);

   if(ip == 1){


      phsp = phase(all_tr, trackn, &trphn, f1, f2);

      if(phsp == NULL){

         printf("****WARNING****, no dynamic variables have been computed, \r\n"
                "                 e.g., speed, velocity, growth rates.     \n\n");


      }


/* write phase speed tracks to file */

      nftmp = nf;
      nfltmp = nfld;
      nf = 0;
      nfld = 0;
      ph = open_file(phasef, "w");
      meantrd(ph, phsp, trphn, 'v');
      close_file(ph, phasef);
      nf = nftmp;
      nfld = nfltmp;
   }

   else if(ip == 2) {

      printf("what is the name of the phase speed file to be read\n");

      scanf("%s", charin);

      ph = open_file(charin, "r");
      phsp = read_tracks(ph, &trphn, &gpr, &tpr, 'v', &gr->alat, &gr->alng);
      close_file(ph, charin);

      if(tpr){

         printf("***ERROR***, speed data defined on the wrong projection\r\n"
                "             re-create on the correct domain.          \n\n");
         exit(1);

      }

   }

   tpr = gr->prty;
   gpr = gr->prgr;

/* display phase speeds */

   if(phsp){

      printf("do you wish to display the phase speeds, 'y' or 'n'\n");

      scanf("\n");
      if(getchar() == 'y'){

         for(i=0; i < trphn; i++) tot_fet += (phsp + i)->num;

         splice_plot(phsp, tot_fet, trphn, 1, 0, 'v');

         ++dcount;

      }

   }


   printf("do you wish to compute a new set of statistics                \r\n"
          "or display an existing set of statistics                      \r\n"
          "Input '0' for no                                              \r\n"
          "Input '1' to compute a new set of statistics                  \r\n"
          "Input '2' to display an existing set of statistics            \r\n"
          "          and/or scale to number density and display          \n\n");

   scanf("%d", &ifmp);

   if(ifmp == 1){

      proj_report(gr->prty, gr->prgr);

      printf("do you wish to compute sample points on a different projection,\r\n"
             "'y' or 'n'                                                     \n\n");

      scanf("\n");
      prty = getchar();

      if(prty != 'y'){
         strncpy(gr->prgr_nm, "Cylindrical", MXPRCH);
         strncpy(gr->prty_nm, "Plate Caree", MXPRCH);

         printf("do you want use grid information already read in as the estimation grid, 'y' or 'n'\n\n");
         scanf("\n");
         if(getchar() == 'y'){
            igf = 1;
            nbx = gr->ix;
            nby = gr->iy;
            gstat1 = gr;
            xt1 = *(gstat1->xgrid);
            yt1 = *(gstat1->ygrid);
            xa1 = xr1 = *(gstat1->xgrid);
            xa2 = xr2 = *(gstat1->xgrid + gstat1->ix - 1);
            ya1 = yr1 = *(gstat1->ygrid);
            ya2 = yr2 = *(gstat1->ygrid + gstat1->iy - 1);

         }

      }

      if(!igf){

         printf("***WARNING***, asking for more sampling points than is neccessary\r\n"
                "               for a smooth contour plot can seriously affect    \r\n"
                "               the performance of the estimators, particulerly   \r\n"
                "               for track density and mean lifetime estimation.   \n\n");
   
         if(prty == 'y') proj_group(-1, -1);

         if(gr->prgr){

           xmn = *(gr1->xgrid);
           xmx = *(gr1->xgrid + gr1->ix - 1);
           ymn = *(gr1->ygrid);
           ymx = *(gr1->ygrid + gr1->iy - 1);

           printf("***Assign points for first region if using azimuthal projection***\n\n");

         }

         select_region(&xr1, &xr2, &yr1, &yr2);

         xa1 = xr1;
         xa2 = xr2;
         ya1 = yr1;
         ya2 = yr2;

         ta = (xr2 - xr1) * (yr2 - yr1);

         printf("total area is %f, how many sample points are required\n", ta);
         scanf("%d", &nb);

         bl = (float)sqrt((double)(ta/nb));

         if(bl > (xr2 - xr1) || bl > (yr2 - yr1)){
            printf("***error***, number of sample points is incompatable  \r\n"
                   "             with area, exiting from statistical analysis\n");
            return;
         }

         nbx = ((xr2 - xr1) / bl);
         nby = ((yr2 - yr1) / bl);

         xt1 = xr1 + (bl / 2.0);
         yt1 = yr1 + (bl / 2.0);

         gstat1 = (GRID * )malloc_initl(sizeof(GRID));
         mem_er((gstat1 == NULL) ? 0 : 1, sizeof(GRID));
         if(gr1)
            memcpy(gstat1, gr1, sizeof(GRID)); 
         else 
            memcpy(gstat1, gr, sizeof(GRID));

         gstat1->ix = nbx;
         gstat1->iy = nby;
         gstat1->xgrid = (float *)calloc(nbx, sizeof(float));
         mem_er((gstat1->xgrid == NULL) ? 0 :1, nbx * sizeof(float));
         gstat1->ygrid = (float *)calloc(nby, sizeof(float));
         mem_er((gstat1->ygrid == NULL) ? 0 :1, nby * sizeof(float));

         for(i=0; i < nbx; i++) *(gstat1->xgrid + i) = xt1 + i * bl;
         for(i=0; i < nby; i++) *(gstat1->ygrid + i) = yt1 + i * bl;

      }

/* assign memory for sample points */

      trsav = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
      mem_er((trsav == NULL) ? 0 : 1, sizeof(struct tot_stat));

      trsav->ptnum = nbx * nby;

      trsav->ptst = (struct pt_stat * )calloc(trsav->ptnum, sizeof(struct pt_stat));
      mem_er((trsav->ptst == NULL) ? 0 : 1, trsav->ptnum * sizeof(struct pt_stat));

/* transform back to sphere */

      if(gr->prty && !(gr->prgr)) prj = proj_assign(gr->prty);

      printf("****INFORMATION****, lower left point for statistical estimation grid is\r\n"
             "                     X = %f, Y = %f\n\n", xt1, yt1);

      pst = trsav->ptst;
  
      for(i=0; i < nby; i++){

          ypos = *(gstat1->ygrid + i);

          yinit = ypos;
       
          for(j=0; j < nbx; j++){

             pst->ys = ypos;
             pst->xs = *(gstat1->xgrid + j);
             xinit = pst->xs;

             if(gr->prty){

                switch(gr->prgr){

                   case 0:

                     (*prj)(&pst->xs, 'x', 0);
                     (*prj)(&pst->ys, 'y', 0);
                     xaa = pst->xs;
                     yaa = pst->ys;
                     pst->apos = i * nbx + j;
                     ++pst; ++ioff;
                     break;

                   case 1:

                     if(azimuthal(&pst->ys, &pst->xs, gr1->alat, gr1->alng, 0, gr->prty)){
                        xaa = pst->xs;
                        yaa = pst->ys;
                        pst->apos = i * nbx + j;
                        ++pst; ++scc; ++ioff;
                     }
                     break;
		     
		   case 2:
		     xat1 = pst->xs;
		     yat1 = pst->ys;		   
		     conic(&pst->ys, &pst->xs, gr1->alat, gr1->alng, gr1->sp1, gr1->sp2, 0, gr->prty, &ift);
                     xat2 = pst->xs;
		     yat2 = pst->ys;	     
		     conic(&yat2, &xat2, gr1->alat, gr1->alng, gr1->sp1, gr1->sp2, 1, gr->prty, &ift);		     
                     if(fabs(xat1 - xat2) < TOLCONIC && fabs(yat1 - yat2) < TOLCONIC){	
		        xaa = pst->xs;
                        yaa = pst->ys;		        	     
                        pst->apos = i * nbx + j;
                        ++pst; ++scc; ++ioff;
	             }
                     break;
		    case 3:
		     rotated_cyln(&(pst->ys), &(pst->xs), gr1->alat, gr1->alng, 0, &ift);	
		     xaa = pst->xs;
                     yaa = pst->ys;		        	     
                     pst->apos = i * nbx + j;
                     ++pst; ++scc; ++ioff;		     
		     break;
                }

                if(!ii){xa1 = xa2 = xaa; ya1 = ya2 = yaa; ii = 1;}
                else{

                    if(xaa < xa1) xa1 = xaa;
                    else if(xaa > xa2) xa2 = xaa;

                    if(yaa < ya1) ya1 = yaa;
                    else if(yaa > ya2) ya2 = yaa;

                }

             }

             else {pst->apos = i * nbx + j; ++pst; ++ioff;}            

          } 

      }


      printf("****INFORMATION****, upper right point for statistical estimation grid is\r\n"
             "                     X = %f, Y = %f\n\n", xinit, yinit);


      if(gr->prgr){

         trsav->ptnum = scc;

         trsav->ptst = (struct pt_stat * )realloc_n(trsav->ptst, scc * sizeof(struct pt_stat));
         mem_er((trsav->ptst == NULL) ? 0 : 1, scc * sizeof(struct pt_stat));

      }

      if(gr2) {

        xmn = *(gr1->xgrid);
        xmx = *(gr1->xgrid + gr1->ix - 1);
        ymn = *(gr1->ygrid);
        ymx = *(gr1->ygrid + gr1->iy - 1);

        printf("***Assign points for second region if using azimuthal projection***\n\n");

        select_region(&xr1, &xr2, &yr1, &yr2);

        ta = (xr2 - xr1) * (yr2 - yr1);

        printf("total area is %f, how many sample points are required\n", ta);
        scanf("%d", &nb);

        bl = (float)sqrt((double)(ta/nb));

        if(bl > (xr2 - xr1) || bl > (yr2 - yr1)){
           printf("***error***, number of sample points is incompatable  \r\n"
                  "             with area, exiting from statistical analysis\n");
           return;
        }

        nbx = (xr2 - xr1) / bl;
        nby = (yr2 - yr1) / bl;

        xt1 = xr1 + (bl / 2.0);
        yt1 = yr1 + (bl / 2.0);

        trsav->ptnum += (nbx * nby);

        trsav->ptst = (struct pt_stat * )realloc_n(trsav->ptst, trsav->ptnum * sizeof(struct pt_stat));
        mem_er((trsav->ptst == NULL) ? 0 : 1, trsav->ptnum * sizeof(struct pt_stat));

        gstat2 = (GRID * )malloc_initl(sizeof(GRID));
        mem_er((gstat2 == NULL) ? 0 : 1, sizeof(GRID));
        memcpy(gstat2, gr2, sizeof(GRID));
        gstat2->ix = nbx;
        gstat2->iy = nby;
        gstat2->xgrid = (float *)calloc(nbx, sizeof(float));
        mem_er((gstat2->xgrid == NULL) ? 0 :1, nbx * sizeof(float));
        gstat2->ygrid = (float *)calloc(nby, sizeof(float));
        mem_er((gstat2->ygrid == NULL) ? 0 :1, nby * sizeof(float));

        for(i=0; i < nbx; i++) *(gstat2->xgrid + i) = xt1 + i * bl;
        for(i=0; i < nby; i++) *(gstat2->ygrid + i) = yt1 + i * bl;

        pst = trsav->ptst + scc;
  
        for(i=0; i < nby; i++){

            ypos = *(gstat2->ygrid + i);
       
            for(j=0; j < nbx; j++){

               pst->ys = ypos;
               pst->xs = *(gstat2->xgrid + j);

               if(azimuthal(&pst->ys, &pst->xs, gr2->alat, gr2->alng, 0, gr->prty)){
                 if(pst->xs < xa1) xa1 = pst->xs;
                 else if(pst->xs > xa2) xa2 = pst->xs;

                 if(pst->ys < ya1) ya1 = pst->ys;
                 else if(pst->ys > ya2) ya2 = pst->ys;
                 pst->apos = i * nbx + j;
                  ++pst; ++scc;

               }

            }

         }

         trsav->ptnum = scc;

         trsav->ptst = (struct pt_stat * )realloc_n(trsav->ptst, trsav->ptnum * sizeof(struct pt_stat));
         mem_er((trsav->ptst == NULL) ? 0 : 1, trsav->ptnum * sizeof(struct pt_stat));

      }

      if(gr->prty) proj_group(0, 0);

      xr1 = xa1;
      xr2 = xa2;
      yr1 = ya1;
      yr2 = ya2;

      xmn = *(gr->xgrid + x1u - 1);
      xmx = *(gr->xgrid + x2u - 1);
      ymn = *(gr->ygrid + y1u - 1);
      ymx = *(gr->ygrid + y2u - 1);         



      trsav->xa1 = xr1;
      trsav->xa2 = xr2;
      trsav->ya1 = yr1;
      trsav->ya2 = yr2;

      printf("****INFORMATION****, grid size for statistical estimation is\r\n"
             "                     nx = %d, ny = %d\n\n", nbx, nby);

/* currently only statistics on the sphere can be computed */

      printf("***WARNING***, kernal estimation on the sphere alone has been\r\n"
             "               implemented so far, hence all data must be    \r\n"
             "               defined on the sphere, and the                \r\n"
             "               geodesic distance measure used.               \n\n");

      proj_report(gr->prty, gr->prgr);

      if(tom != 'g') {

          tom = 'g';

          printf("***Distance measure is now GEODESIC***\n\n");

      }

      if(gr->prty) {

         printf("***ERROR***, track data must be defined on a lat-long domain\r\n"
                "             re-create data on correct domain               \n\n");

         exit(1);

      }

/* convert estimation grid to Cartesians */

      stat = grid_convert(trsav);

/* convert data latitude-longitudes to Cartesian space */

      convert_track(all_tr, trackn, 0, 0);
      if(phsp) convert_track(phsp, trphn, 0, 0);

/* extract data points for track, genesis and lysis density estimation */

      printf("The current region is X = %f, %f, Y = %f %f\n\n", xr1, xr2, yr1, yr2);
      printf("The current full region is X = %f, %f, Y = %f %f\n\n", xmn, xmx, ymn, ymx);

      printf("What is the actual region required? Input x1, x2, y1, y2\n\n");
      scanf("%f %f %f %f", &xr1, &xr2, &yr1, &yr2);

      if(xr1 < GXMN) xr1 = GXMN;
      if(xr2 > GXMX) xr2 = GXMX;
      if(yr1 < GYMN) yr1 = GYMN;
      if(yr2 > GYMX) yr2 = GYMX;


      trsav->xa1 = xr1 * FP_PI;
      trsav->xa2 = xr2 * FP_PI;
      trsav->ya1 = FP_PI2 - yr1 * FP_PI;
      trsav->ya2 = FP_PI2 - yr2 * FP_PI;

      printf("the search region is X = %f, %f, Y = %f %f\n\n", xr1, xr2, yr1, yr2);

/* Set boundary flags for boundary kernels */

/*      if(fabs(trsav->xa1 - GXMN * FP_PI) > TOLBOUND) {trsav->ibb[0] = 1; trsav->ibound = 1;}
      if(fabs(trsav->xa2 - GXMX * FP_PI) > TOLBOUND) {trsav->ibb[1] = 1; trsav->ibound = 1;}
      if(fabs(trsav->ya1 - (FP_PI2 - GYMN * FP_PI)) > TOLBOUND) {trsav->ibb[2] = 1; trsav->ibound = 1;}
      if(fabs(trsav->ya2 - (FP_PI2 - GYMX * FP_PI)) > TOLBOUND) {trsav->ibb[3] = 1; trsav->ibound = 1;}



      if(trsav->ibound){
         printf("Kernel estimation will in be the presence of actual boundaries, \r\n"
                "do you want to use boundary kernel reflection to prevent,       \r\n"
                "leakage of probability mass outside of region, 'y' or 'n'.      \n\n");
         scanf("\n");
         if(getchar() == 'n') trsav->ibound = 0;

      } */

/* check sample period with tracks */

      fnc = 0;

      for(i=0; i < trackn; i++){

         altr = all_tr + i;

         if(altr->trpt != NULL){

           if(!fnc) fnc = (altr->trpt + altr->num - 1)->fr_id;

           at = altr->trpt;
           fe =  (at + (altr->num - 1))->fr_id;

           if(fe > fnc) fnc =  fe;

         }


      }

      ff1 = 1;
      ff2 = fnc;

      if(fnc < f2) f2 = fnc;

      if(f1 < 1) f1 = 1;

      printf("***INFORMATION***, actual frame range is, %d --> %d frames\n\n", ff1, ff2);
      printf("***INFORMATION***, chosen frame range spans, %d --> %d frames\n\n", f1, f2);   


      for(i=0; i<trackn; i++){

          if((all_tr + i)->awt){awt = (all_tr + i)->awt; break;}

      }


      if(awt){

         printf("****INFORMATION****, additional weights for statistical estimation \r\n"
                "                     have been detected for this ensemble of tracks\n\n");
         printf("Are the additional weights required for the statistical estimation,\r\n"
                "e.g., time dependent weighting, 'y' or 'n'.                        \n\n");
         scanf("\n");
         if(getchar() == 'y'){

            if(awt == 1) {

               printf("Specify a lower tolerance on the weights, this will \r\n"
                      "speed up the statistical estimation when many of the\r\n"
                      "weights are small. Value must be positive.          \n\n");
               scanf("%f", &tolwt);

               if(tolwt < 0.0){

                 printf("****WARNING****, tolerance on weights must be positive\r\n"
                      "                 continuing with a tolerance of 0.0   \n\n");
                  tolwt = 0.0;

              }

            }
            else if(awt == 2){

               printf("****INFORMATION****, weights are raw values, these can only be used \r\n"
                      "                     with the SQT estimator to estimate statistics  \r\n"
                      "                     for a particular index value using 1D kernel's \r\n"
                      "                     in index space. Note index should be centered  \r\n"
                      "                     and normalized.                                \n\n"); 

               tele = (TELE *)malloc_initl(sizeof(TELE));
               mem_er((tele == NULL) ? 0 : 1, sizeof(TELE));
               printf("What index value should statistics be estimated at?\n\n");
               scanf("%lf", &(tele->tst));
               printf("What bandwidth is required for the index kernels?\n\n");
               scanf("%lf", &(tele->h));
               printf("What kind of kernel is required for the index dimension: \r\n"
                      "'0' for Triangular.                                      \r\n"
                      "'1' for Rectangular.                                     \r\n"
                      "'2' for Epanechnikov.                                    \r\n"
                      "'3' for Bi-weights.                                      \n\n");
               scanf("%d", &(tele->ktyp));

               if(tele->ktyp == 2) tele->h /= R5;

            }

         }
         else {

           for(i=0; i<trackn; i++){

              (all_tr + i)->awt = 0;

           }

           awt = 0;

         }

      }

/* If shape and area parameters are available, scale to planetary values
   e.g. area.                                                             */

      if(aniso == 'y'){

        sc = 1.0;

        printf("Do you want to scale the area attribute, 'y' or 'n'\n\n");
        scanf("\n");
        if(getchar() == 'y'){
           printf("What type of value do you want                \r\n"
                  "Input '0' for user value                      \r\n"
                  "Input '1' for Earths radius squared in Km^2   \r\n");
           scanf("%d", &isc);

           if(!isc) {

              printf("input scaling\n");
              scanf("%f", &sc);

            }

            else if (isc == 1) sc = EARTH_RADIUS2;

        }

      }
      

/* assign memory for track lifetimes */

      ltm = (float * )calloc(trackn, sizeof(float));
      mem_er((ltm == NULL) ? 0 : 1, trackn * sizeof(float));

      printf("what is the time step in the units required for lifetimes, i.e. s, hr, day etc ?\n\n");

      scanf("%f", &ts);

      strncpy(charin, INITTRS, MAXCHR);
      if(iext) strcpy(strstr(charin, EXTENSION), fext);

      gend = open_file(charin, "w");

      strncpy(charin, DISPTRS, MAXCHR);
      if(iext) strcpy(strstr(charin, EXTENSION), fext);

      lysd = open_file(charin, "w");

      fprintf(gend, "%-5d\n", 0);
      fprintf(lysd, "%-5d\n", 0);

      for(i=0; i < trackn; i++){

         altr = all_tr + i;

         if(altr->trpt != NULL){

           at = altr->trpt;

/* how many missing frames for current track */

           nmiss = 0;
           if(!imiss) for(j=0; j<altr->num; j++) nmiss += (at+j)->nfm;

           fb =  at->fr_id;
           fe =  (at + (altr->num - 1))->fr_id;

           if(at->fr_id > f1 && at->fr_id <= f2){

              ddd = dgnn;

              dgn = extract_data(at, dgn, &dgnn, xr1, xr2, yr1, yr2, 's', &wt, 0);
              if(altr->awt && dgnn > ddd){
                if(!tele){
                   if(wt > tolwt) wdgn = weights(wdgn, dgnn, wt, 0);
                   else --dgnn;
                }
                else wdgn = weights(wdgn, dgnn, wt, 0);
              }

              if(dgnn > ddd){

                 for(j=0; j < trphn; j++){

                     altr2 = phsp + j; 

                     if(altr2->tr_sp == i){

                        spav = 0;
                        at2 = altr2->trpt;

                        npav = 0;
                        for(k=0; k < altr2->num; k++) {
                           if((at2 + k)->zf < ADD_CHECK){
                              spav += (at2 + k)->zf;
                              ++npav;
                           }
                        }

                        spav /= npav;


                     }                   


                 }

                 fprintf(gend, "%-8d %7.2f %7.2f %e %6.2f %6.2f\n", at->fr_id, at->xf, at->yf, at->zf, (fe - fb + nmiss) * ts, spav);

              }

           }


           d1 = d2 = dd=0;

           if(fb >= ff1 && fe <= ff2){

              if(fe < f1 && fb < f1) *(ltm + i) = 0.0;
              else if(fe > f2 && fb > f2)*(ltm + i) = 0.0;
              else {
                 dd = fe - fb + nmiss;
                 if(fb < f1 && fe >= f1) d1 = f1 - fb;
                 if(fe > f2 && fb <= f2) d2 = fe - f2;

                 if((float)(d2 + d1)/(float)dd < LPERC){
        
                    *(ltm + i) = (fe - fb + nmiss) * ts;
                    ++nlif;
                 }

                 else *(ltm + i) = 0.0;

              }

           }

           else *(ltm + i) = 0.0;

           at += altr->num - 1;

           if(at->fr_id >= f1 && at->fr_id < f2) {

              ddd = dlyn;
     
              dly = extract_data(at, dly, &dlyn, xr1, xr2, yr1, yr2, 's', &wt, 0);
              if(altr->awt && dlyn > ddd){
                if(!tele){
                   if(wt > tolwt) wdly = weights(wdly, dlyn, wt, 0);
                   else --dlyn;
                }
                else wdly = weights(wdly, dlyn, wt, 0);
              }

              if(dlyn > ddd){


                 for(j=0; j < trphn; j++){

                     altr2 = phsp + j;  
                     if(altr2->tr_sp == i){

                        spav = 0;
                        at2 = altr2->trpt;

                        npav = 0;
                        for(k=0; k < altr2->num; k++) {
                            if((at2 + k)->zf < ADD_CHECK){
                               spav += (at2 + k)->zf;
                               ++npav;
                            }
                        }

                        spav /= npav;


                     }                   


                 }


                 fprintf(lysd, "%-8d %7.2f %7.2f %e %6.2f %6.2f\n", at->fr_id, at->xf, at->yf, at->zf, (fe - fb + nmiss) * ts, spav);

              }

           }

/* ****************************************** */

/* code for computing max stats

           istmax = 0;
           strmax = altr->trpt->zf;
           for(j=0; j < altr->num; j++){
               at = (altr->trpt) + j;
               if(at->zf > strmax) {istmax = j; strmax = at->zf;}
           }

           at = (altr->trpt) + istmax;
           if(at->fr_id >= f1 && at->fr_id <= f2){
              ddd = dtrn;
              dtr = extract_data(at, dtr, &dtrn, xr1, xr2, yr1, yr2, 's', &wt, 0);
           }

*/

/* ****************************************** */

           for(j=0; j < altr->num; j++){

               at = (altr->trpt) + j;

               if(at->fr_id >= f1 && at->fr_id <= f2){

                  if(at->zf < ADD_CHECK){

                     ddd = dtrn;

                     dtr = extract_data(at, dtr, &dtrn, xr1, xr2, yr1, yr2, 's', &wt, 0);

                     if(altr->awt && dtrn > ddd){ 
                        if(!tele){
                           if(wt > tolwt) wdtr = weights(wdtr, dtrn, wt, 0);
                           else --dtrn;
                        }
                        else wdtr = weights(wdtr, dtrn, wt, 0);
                     }

                  }

                  if(aniso == 'y'){
                    at->area *= sc; 
                    ddd = dann;
                    dan = extract_data(at, dan, &dann, xr1, xr2, yr1, yr2, 'a', &wt, 0);
                    dar = extract_data(at, dar, &darn, xr1, xr2, yr1, yr2, 'r', &wt, 0);
                    if(altr->awt && dann > ddd){
                      if(!tele){
                         if(wt > tolwt) {
                            wdan = weights(wdan, dann, wt, 0);
                            wdar = weights(wdar, darn, wt, 0);
                         }
                         else {--dann; --darn;}
                      }
                      else {
                         wdan = weights(wdan, dann, wt, 0);
                         wdar = weights(wdar, darn, wt, 0);
                      }

                    }

                  }

               }

            } 

         }

      }

      fseeko(gend, (off_t)0, FSTART);
      fprintf(gend, "%-8d", dgnn);

      fseeko(lysd, (off_t)0, FSTART);
      fprintf(lysd, "%-8d", dlyn);

      close_file(gend, charin);
      close_file(lysd, charin);

      printf("***INFORMATION***, number of tracks satisfying the frame \r\n"
             "                   constraints for track density and mean\r\n"
             "                   lifetime is %d\n\n", nlif);

/* extract propogation speeds */

      if(phsp != NULL){

         printf("For growth/decay rates and tendencies, do you want:            \r\n"
                "mixed positive and negative growths and tendencies , input '0' \r\n"
                "only positive growth and tendencies                , input '1' \r\n"
                "only negative growth and tendencies                , input '2' \n\n");

         scanf("%d", &igrwth);
         if(igrwth < 0 || igrwth > 2){
            printf("****WARNING****, incorrect identifier %d for growth/decay and tendencies,  \r\n"
                   "                 using default, both positive and negative values are used.\n\n", igrwth);
            igrwth = 0;
         }

         for(i=0; i < trphn; i++){

             altr = phsp + i;

/* ****************************************** */
/* code for computing max stats

             istmax = 0;
             strmax = altr->trpt->gwthr;
             for(j=0; j < altr->num; j++){
                 at = (altr->trpt) + j;
                 if(at->gwthr > strmax) {istmax = j; strmax = at->gwthr;}
             }

             at = (altr->trpt) + istmax;
             if(at->fr_id >= f1 && at->fr_id <= f2){
                ddd = dgrn;
                dgr = extract_data(at, dgr, &dgrn, xr1, xr2, yr1, yr2, 'g', &wt, 0);
             }

             istmax = 0;
             strmax = altr->trpt->tend;
             for(j=0; j < altr->num; j++){
                 at = (altr->trpt) + j;
                 if(at->tend > strmax) {istmax = j; strmax = at->tend;}
             }

             at = (altr->trpt) + istmax;
             if(at->fr_id >= f1 && at->fr_id <= f2){
                ddd = dtnn;
                dtn = extract_data(at, dtn, &dtnn, xr1, xr2, yr1, yr2, 't', &wt, 0);
             }

*/

/* ****************************************** */

             for(j=0; j < altr->num; j++) {

                 at = (altr->trpt) + j;

                 ngrwth = 0;

                 switch(igrwth){
                     case 0:
                        ngrwth = 1;
                        break;
                     case 1:
                        ngrwth = (at->gwthr > 0.0) ? 1: 0;
                        break;
                     case 2:
                        ngrwth = (at->gwthr < 0.0) ? 1: 0;
                        break;
                     default:
                        ngrwth = 1;
                 }

                 if(at->zf < ADD_CHECK){
                    ddd = dphn;
                    dph = extract_data(at, dph, &dphn, xr1, xr2, yr1, yr2, 'v', &wt, 0);
                    if(altr->awt && dphn > ddd){
                      if(!tele){
                         if(wt > tolwt) wdph = weights(wdph, dphn, wt, 0);
                         else --dphn;
                      }
                      else wdph = weights(wdph, dphn, wt, 0);
                    }

                 }

                 if(ngrwth && at->gwthr < ADD_CHECK){

                    ddd = dgrn;
                    dgr = extract_data(at, dgr, &dgrn, xr1, xr2, yr1, yr2, 'g', &wt, 0);
                    if(altr->awt && dgrn > ddd){
                      if(!tele){
                         if(wt > tolwt) wdgr = weights(wdgr, dgrn, wt, 0);
                         else --dgrn;
                      }
                      else wdgr = weights(wdgr, dgrn, wt, 0);
                    }

                    ddd = dtnn;
                    dtn = extract_data(at, dtn, &dtnn, xr1, xr2, yr1, yr2, 't', &wt, 0);
                    if(altr->awt && dtnn > ddd){
                      if(!tele){
                         if(wt > tolwt) wdtn = weights(wdtn, dtnn, wt, 0);
                         else --dtnn;
                      }
                      else wdtn = weights(wdtn, dtnn, wt, 0);

                    }

                 } 

             }

          }

       }

/* choose statistics to compute */

   printf("You can choose to compute all statistics or a selection of statistics,\r\n"
          "track density and mean lifetime are selected in a different way so    \r\n"
          "at this stage selecting them has no effect.                           \n\n");

   printf("Do you want all statistics, input 'a', or a selection, input 's'\n\n");
   scanf("\n");
   if(getchar() == 's'){
      available_stats();
      while(1){
         scanf("%d", &styp);
         if(!styp) break;
         else if(styp == 1){free(dtr); free(wdtr); dtrn = 0;}
         else if(styp == 2){free(dph); free(wdph); dphn = 0;}
         else if(styp == 3){free(dgn); free(wdgn); dgnn = 0;} 
         else if(styp == 4){free(dly); free(wdly); dlyn = 0;}   
         else if(styp == 5){free(dgr); free(wdgr); dgrn = 0;} 
         else if(styp == 6){free(dan); free(wdan); dann = 0;}
         else if(styp == 7){free(dtn); free(wdtn); dtnn = 0;}
         else if(styp == 8){free(dar); free(wdar); darn = 0;}
         else {
            printf("Unknown statistic Id.\n\n");
         }
      }
   }

   else {

      printf("****INFORMATION****, computing all available statistics\n\n");

   }


/* Choose to use the simple but computationally expensive extimation or use the
   more efficient estimation based on spherical quad tree partitioning of the data */

   printf("Do you want to use the old estimation procedures, or use the more      \r\n"
          "efficient procedures based on Spherical Quad Tree's (SQT). SQT         \r\n"
          "should be used if computing confidence intervals by re-sampling.       \r\n"
          "                                                                       \r\n"
          "****WARNING****, SQT can only currently be used with isotropic kernels.\r\n"
          "                                                                       \r\n"
          "Input '0' for old estimation and '1' for SQT.                          \n\n");

   scanf("%d", &isamp);

   if(isamp){

      printf("What base tesselation is required, '0' for Octahedron, or '1' for Icosahedron.\r\n"
             "Default is Octahedron.                                                        \n\n");

      scanf("%d", &iplat);

      sqth = (SQT **)malloc_initl(sizeof(SQT *));
      mem_er((sqth == NULL) ? 0 : 1, sizeof(SQT *));

      if(iplat < 0 || iplat > 1) iplat = 0;

      if(!iplat) {
         printf("****INFORMATION****, using Octahedron for base tesselation.\n\n");

         nfc = 8;
         nvec = 6;
         *sqth = oct;
         gvecs = (VEC *)calloc(nvec, sizeof(VEC));
         mem_er((gvecs == NULL) ? 0 : 1, nvec * sizeof(VEC));
         memcpy(gvecs, oct_vecs, nvec * sizeof(VEC));

      }
      else if(iplat == 1) {
         printf("****INFORMATION****, using Icosahedron for base tesselation.\n\n");

         nfc = 20;
         nvec = 12;
         *sqth = icos;
         gvecs = (VEC *)calloc(nvec, sizeof(VEC));
         mem_er((gvecs == NULL) ? 0 : 1, nvec * sizeof(VEC));
         memcpy(gvecs, icos_vecs, nvec * sizeof(VEC));

      }

      rtar = area_sp_triangle(*sqth, gvecs, 0);

/* for(i=0; i < nfc; i++) printf("%f\n", area_sp_triangle(*sqth + i, gvecs, 0)); */


      printf("****INFORMATION*****, area of a root spherical triangle is %f in steradians\n\n", rtar);

      printf("Specify either number of levels or approximate areal size of triangles \r\n" 
             "for SQT decomposition.                                                 \n\n");
      printf("Input number of levels in the hierarchy for each root spherical triangle,\r\n"
             "input '0' if a represenative triangle area will be specified.            \n\n");
      scanf("%d", &nlev); 

      if(nlev > 0){
         ntrang = powi(4, nlev-1);
         printf("****INFORMATION*****, there are %d leaf spherical triangles of  \r\n"
                "                      approximate area %e to each root triangle.\n\n", ntrang, rtar / (float)ntrang);

      }

      else {
         printf("Specify target representative triangle area for tesselation in steradians. \n\n");
         scanf("%f", &tarea);

/* compute required level of hierarchy from chosen area */

         ntrang = rtar / tarea;
         printf("%d\n", ntrang);

         nlev = floor(exp(log((float)ntrang) / 4) + 0.5);

         printf("****INFORMATION****, the number of levels for chosen area is %d        \r\n"
                "                     with %d leaf spherical triangles per root triangle.\n\n", nlev, powi(4, nlev-1));

      }

      if(nlev <= 1) {
         printf("****ERROR****, number of levels in the SQT hierarchy must be greater than 1\n\n");
         exit(1); 
      }

/* need to add read/write to SQT dump */

/* assign rest of memory for SQT and generate tree */

      sqth = (SQT **)realloc_n(sqth, nlev * sizeof(SQT *));
      mem_er((sqth == NULL) ? 0 : 1, nlev * sizeof(SQT *));

/* ***************DEBUG************** */

/*for(i=0; i < nfc; i++){
    printf("%p ", (void *)(sqth[0] + i));
}*/

/* ********************************** */

      nelm = nfc;
      for(i=1; i < nlev-1; i++){
         printf("Generating level %d of SQT\n", i);
         sqth[i] = (SQT *)calloc(4 * nelm, sizeof(SQT));
         mem_er((sqth[i] == NULL) ? 0 : 1, 4 * nelm * sizeof(SQT));
         gvecs = create_sqt(sqth, gvecs, nelm, i, &nvec, nlev);
         nelm *= 4;

      }

      printf("Generating leaf node level of SQT.\n\n");
      sqth[nlev-1] = (SQT *)((LEAF *)calloc(4 * nelm, sizeof(LEAF)));
      mem_er((sqth[nlev-1] == NULL) ? 0 : 1, 4 * nelm * sizeof(LEAF));
      gvecs = create_sqt(sqth, gvecs, nelm, nlev-1, &nvec, nlev);
      lf = (LEAF *)sqth[nlev - 1];
      nelm *= 4;


/* ***************DEBUG************** */

/*for(i=0; i < nelm; i++){
    printf("%p ", (void *)((LEAF *)sqth[nlev-1] + i));
}
printf("\n");*/


/* ********************************** */


/* determine neighbours in the SQT */

      sqt_neighbour_find(sqth, nlev, nelm, iplat, nfc);

/* determine grid point ownership in the SQT */

      data_partition(sqth, gvecs, nfc, nlev, 0, stat, trsav->ptnum);

   }

   

/* optional single density estimation of the smoothing parameter
   or individual density and regression estimations are possible,
   it is the users choice.                                        */

      times(&usage);
      tstartall = usage.tms_utime + usage.tms_stime;

/* assign memory for statistic calculation */


      den = (double ** )calloc(5, sizeof(double *));
      mem_er((den == NULL) ? 0 : 1, 5 * sizeof(double *));

      den[0] = (double * )calloc(trsav->ptnum, sizeof(double));
      mem_er((den[0] == NULL) ? 0 : 1, trsav->ptnum * sizeof(double));

      den[1] = (double * )calloc(trsav->ptnum, sizeof(double));
      mem_er((den[1] == NULL) ? 0 : 1, trsav->ptnum * sizeof(double));

      den[2] = (double * )calloc(trsav->ptnum, sizeof(double));
      mem_er((den[2] == NULL) ? 0 : 1, trsav->ptnum * sizeof(double));


/* compute track densities and mean strength */

      nn = 0.0;

      printf("***INFORMATION***, performing PHENOMENA density estimation\r\n"
             "                   and MEAN STRENGTH\n\n");

/* if using SQT then partition data */

      if(isamp && dtrn) data_partition(sqth, gvecs, nfc, nlev, 1, dtr, dtrn);

      printf("do you want a single density estimation of the smoothing\r\n"
             "parameter, or individual density and regression estimations\r\n"
             "input 's' for single, or 'i' otherwise\n\n");

      scanf("\n");

      if(getchar() == 's'){

         trsav->sm[0] = sp_kernal_estimate(den, stat, dtr, wdtr, dtrn, trsav->ptnum, 2, &(trsav->kern[0]), 0, &nn, lf, gvecs, nelm, 1, tele);
         trsav->sm[4] = trsav->sm[1] = trsav->sm[0];
         trsav->kern[4] = trsav->kern[1] = trsav->kern[0];

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            pst->stat3 = *(den[0] + i);
            (pst->stat1).mean = *(den[1] + i);
            (pst->stat1).var = *(den[2] + i);
            *(den[0] + i) = *(den[1] + i) = *(den[2] + i) = 0.0;
         }

      }

      else {

         printf("***INFORMATION***, performing PHENOMENA density estimation\n\n");

         trsav->sm[4] = sp_kernal_estimate(den, stat, dtr, wdtr, dtrn, trsav->ptnum, 1, &(trsav->kern[4]), 0, &nn, lf, gvecs, nelm, 1, tele);

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            pst->stat3 = *(den[0] + i);

            *(den[0] + i) = 0.0;
         }         

         printf("***INFORMATION***, performing MEAN STRENGTH estimation\n\n");

         trsav->sm[0] = sp_kernal_estimate(den, stat, dtr, wdtr, dtrn, trsav->ptnum, 0, &(trsav->kern[0]), 0, &nn, lf, gvecs, nelm, 1, tele);
         trsav->sm[1] = trsav->sm[0];
         trsav->kern[1] = trsav->kern[0];

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            (pst->stat1).mean = *(den[1] + i);
            (pst->stat1).var = *(den[2] + i);

            *(den[0] + i) = *(den[1] + i) = *(den[2] + i) = 0.0;
         }

      }

      if(wdtr) nn /= sum_wt;

      trsav->datnm[0] = trsav->datnm[1] = trsav->datnm[4] = nn ;
       
      free(dtr);
      free(wdtr);

/* compute mean speeds */

      den[3] = (double * )calloc(trsav->ptnum, sizeof(double));
      mem_er((den[3] == NULL) ? 0 : 1, trsav->ptnum * sizeof(double));

      den[4] = (double * )calloc(trsav->ptnum, sizeof(double));
      mem_er((den[4] == NULL) ? 0 : 1, trsav->ptnum * sizeof(double));

      nn = 0.0;

      printf("***INFORMATION***, performing MEAN SPEED estimation, and\r\n"
             "                   TRACK FLUX estimation\n\n");

      if(isamp && dphn) data_partition(sqth, gvecs, nfc, nlev, 1, dph, dphn);

      printf("do you want a single density estimation of the smoothing\r\n"
             "parameter, or individual density and regression estimations\r\n"
             "input 's' for single, or 'i' otherwise\n\n");

      scanf("\n");
      if(getchar() == 's'){


         trsav->sm[2] = sp_kernal_estimate(den, stat, dph, wdph, dphn, trsav->ptnum, 3, &(trsav->kern[2]), 0, &nn, lf, gvecs, nelm, 1, tele);
         trsav->sm[8] = trsav->sm[3] = trsav->sm[2];
         trsav->kern[8] = trsav->kern[3] = trsav->kern[2];

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            (pst->stat2).mean = *(den[1] + i);
            (pst->stat2).var = *(den[2] + i);
            (pst->stat7).xcomp = *(den[3] + i);
            (pst->stat7).ycomp = *(den[4] + i);

            *(den[0] + i) = 0.0;
            *(den[1] + i) = *(den[2] + i) = *(den[3] + i) = *(den[4] + i) = 0.0;

         }

      }

      else {

         printf("***INFORMATION***, performing MEAN SPEED estimation.\n\n");


         trsav->sm[2] = sp_kernal_estimate(den, stat, dph, wdph, dphn, trsav->ptnum, 0, &(trsav->kern[2]), 0, &nn, lf, gvecs, nelm, 1, tele);
         trsav->sm[3] = trsav->sm[2];
         trsav->kern[3] = trsav->kern[2];

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            (pst->stat2).mean = *(den[1] + i);
            (pst->stat2).var = *(den[2] + i);

            *(den[0] + i) = 0.0;
            *(den[1] + i) = *(den[2] + i) = 0.0;

         }

         printf("***INFORMATION***, performing TRACK FLUX estimation\n\n");

         trsav->sm[8] = sp_kernal_estimate(den, stat, dph, wdph, dphn, trsav->ptnum, -1, &(trsav->kern[8]), 0, &nn, lf, gvecs, nelm, 1, tele);

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            (pst->stat7).xcomp = *(den[3] + i);
            (pst->stat7).ycomp = *(den[4] + i);

            *(den[0] + i) = 0.0;
            *(den[3] + i) = *(den[4] + i) = 0.0;

         }

      }

      
      if(wdph) nn /= sum_wt;

      trsav->datnm[2] = trsav->datnm[3] = trsav->datnm[8] = nn ;
       
      free(dph);
      free(wdph);


      nn = 0.0;

      printf("***INFORMATION***, performing MEAN ANISOTROPY estimation, and\r\n"
             "                   MEAN ORIENTATION estimation\n\n");

      if(isamp && dann) data_partition(sqth, gvecs, nfc, nlev, 1, dan, dann);

      printf("do you want a single density estimation of the smoothing\r\n"
             "parameter, or individual density and regression estimations\r\n"
             "input 's' for single, or 'i' otherwise\n\n");

      scanf("\n");
      if(getchar() == 's'){

         trsav->sm[11] = sp_kernal_estimate(den, stat, dan, wdan, dann, trsav->ptnum, 3, &(trsav->kern[11]), 0, &nn, lf, gvecs, nelm, 1, tele);
         trsav->sm[12] = trsav->sm[11];
         trsav->kern[12] = trsav->kern[11];

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            pst->stat10 = *(den[1] + i);
            (pst->stat11).xcomp = *(den[3] + i);
            (pst->stat11).ycomp = *(den[4] + i);

            *(den[0] + i) = 0.0;
            *(den[1] + i) = *(den[3] + i) = *(den[4] + i) = 0.0;

         }

      }

      else {

         printf("***INFORMATION***, performing MEAN ANISOTROPY estimation.\n\n");

         trsav->sm[11] = sp_kernal_estimate(den, stat, dan, wdan, dann, trsav->ptnum, 0, &(trsav->kern[11]), 0, &nn, lf, gvecs, nelm, 1, tele);

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            pst->stat10 = *(den[1] + i);


            *(den[0] + i) = *(den[1] + i) = 0.0;


         }

         printf("***INFORMATION***, performing MEAN ORIENTATION estimation\n\n");

         trsav->sm[12] = sp_kernal_estimate(den, stat, dan, wdan, dann, trsav->ptnum, -1, &(trsav->kern[12]), 0, &nn, lf, gvecs, nelm, 1, tele);

         for(i=0; i < trsav->ptnum; i++){

            pst = trsav->ptst + i;

            (pst->stat11).xcomp = *(den[3] + i);
            (pst->stat11).ycomp = *(den[4] + i);

            *(den[0] + i) = 0.0;
            *(den[3] + i) = *(den[4] + i) = 0.0;


         }

      }

      
      if(wdan) nn /= sum_wt;

      trsav->datnm[11] = trsav->datnm[12] = nn;
       
      free(dan);
      free(wdan);



      free(den[4]);
      free(den[3]);

/* compute growth/decay rates */

      nn = 0.0;

      printf("***INFORMATION***, performing mean GROWTH/DECAY rates estimation.\n\n");
      

      if(isamp && dgrn) data_partition(sqth, gvecs, nfc, nlev, 1, dgr, dgrn);

      trsav->sm[10] = sp_kernal_estimate(den, stat, dgr, wdgr, dgrn, trsav->ptnum, 0, &(trsav->kern[10]), 0, &nn, lf, gvecs, nelm, 1, tele);

      for(i=0; i < trsav->ptnum; i++){

         pst = trsav->ptst + i;

         pst->stat9 = *(den[1] + i);

         *(den[0] + i) = *(den[1] + i) = 0.0;

      }

      
      if(wdgr) nn /= sum_wt;

      trsav->datnm[10] = nn;

      free(dgr);
      free(wdgr);

/* compute tendencies */

      nn = 0.0;

      printf("***INFORMATION***, performing mean TENDENCY estimation.\n\n");

      if(isamp && dtnn) data_partition(sqth, gvecs, nfc, nlev, 1, dtn, dtnn);

      trsav->sm[13] = sp_kernal_estimate(den, stat, dtn, wdtn, dtnn, trsav->ptnum, 0, &(trsav->kern[13]), 0, &nn, lf, gvecs, nelm, 1, tele);

      for(i=0; i < trsav->ptnum; i++){

         pst = trsav->ptst + i;

         pst->stat12 = *(den[1] + i);

         *(den[0] + i) = *(den[1] + i) = 0.0;

      }

      
      if(wdtn) nn /= sum_wt;

      trsav->datnm[13] = nn;

      free(dtn);
      free(wdtn);
      
/* compute genesis and lysis densities */

      nn = 0.0;

      printf("***INFORMATION***, performing GENESIS density estimation\n\n");

      if(isamp && dgnn) data_partition(sqth, gvecs, nfc, nlev, 1, dgn, dgnn);

      trsav->sm[5] = sp_kernal_estimate(den, stat, dgn, wdgn, dgnn, trsav->ptnum, 1, &(trsav->kern[5]), 0, &nn, lf, gvecs, nelm, 1, tele);


      for(i=0; i < trsav->ptnum; i++){

         (trsav->ptst + i)->stat4 = *(den[0] + i);
         *(den[0] + i) = 0.0;

      }

      
      if(wdgn) nn /= sum_wt;

      trsav->datnm[5] = nn;

      free(dgn);
      free(wdgn);

      nn = 0.0;

      printf("***INFORMATION***, performing LYSIS density estimation\n\n");

      if(isamp && dlyn) data_partition(sqth, gvecs, nfc, nlev, 1, dly, dlyn);

      trsav->sm[6] = sp_kernal_estimate(den, stat, dly, wdly, dlyn, trsav->ptnum, 1, &(trsav->kern[6]), 0, &nn, lf, gvecs, nelm, 1, tele);


      for(i=0; i < trsav->ptnum; i++) {

         (trsav->ptst + i)->stat5 = *(den[0] + i);
         *(den[0] + i) = 0.0;

      }

      
      if(wdly) nn /= sum_wt;

      trsav->datnm[6] = nn;

      free(dly);
      free(wdly);


/* calculate the track density and mean lifetimes */

      times(&usage);
      tstarttd = usage.tms_utime + usage.tms_stime;

      nn = 0.0;

      printf("***INFORMATION***, performing TRACK density and mean LIFETIME\n"
             "                   calculation. For convience use a smoothing\n"
             "                   parameter of similar value to that used   \n"
             "                   for LYSIS and GENESIS.                  \n\n");

      printf("***WARNING***, this can be very expensive do you wish to continue, 'y' or 'n'\n\n");

      scanf("\n");
      if(getchar() == 'y'){

/* resize memory for returned statistic from sp_kernal_estimate */

        den[2] = (double * )realloc_n(den[2], sizeof(double));
        mem_er((den[2] == NULL) ? 0 : 1, sizeof(double));
        den[1] = (double * )realloc_n(den[1], sizeof(double));
        mem_er((den[1] == NULL) ? 0 : 1, sizeof(double));
        den[0] = (double * )realloc_n(den[0], sizeof(double));
        mem_er((den[0] == NULL) ? 0 : 1, sizeof(double));


/* extract data for track density and lifetime calculation */

        trd = 1;

        mfpts = sizeof(struct fet_pt_tr);

        dtd = (struct dpt *)calloc(trackn, sizeof(struct dpt));
        mem_er((dtd == NULL) ? 0 : 1, trackn * sizeof(struct dpt));


        if(awt){
           wdtd = (float *)calloc(trackn, sizeof(float));
           mem_er((dtd == NULL) ? 0 : 1, trackn * sizeof(float));
        }

        for(i=0; i < trsav->ptnum; i++){

          stt = stat + i;
          pst = trsav->ptst + i;

          for(j=0; j < trackn; j++){

              altr = all_tr + j;

              if(altr->trpt != NULL && *(ltm + j) > TTSTEP){

                 for(k=0; k < altr->num; k++){

                     at = altr->trpt + k;

                     arcl = (float) arc(at, stt);

                     if(k == 0) {
                        arcm = arcl;
                        memcpy(&att, at, mfpts);
                     }

                     else if(arcl > arcm) {
                        arcm = arcl;
                        memcpy(&att, at, mfpts);
                     }

                 }

                 att.zf = *(ltm + j);                 

                 ddd = dtdn;

                 extract_data(&att, dtd, &dtdn, xr1, xr2, yr1, yr2, 's', &wt, 1);
                 if(altr->awt && dtdn > ddd){
                    if(!tele){
                       if(wt > tolwt) weights(wdtd, dtdn, wt, 1);
                       else --dtdn;
                    }
                    else weights(wdtd, dtdn, wt, 1);
                 }

              }

            }

            if(dtd) {

              if(isamp && dtdn) {
                 for(j=0; j< nelm; j++) {free((lf + j)->lgrid); (lf + j)->ng = 0;}
                 data_partition(sqth, gvecs, nfc, nlev, 0, stt, 1);
                 data_partition(sqth, gvecs, nfc, nlev, 1, dtd, dtdn);
              }
              trsav->sm[7] = sp_kernal_estimate(den, stt, dtd, wdtd, dtdn, 1, 2, &(trsav->kern[7]), 1, &nn, lf, gvecs, nelm, 0, tele);
              trsav->sm[9] = trsav->sm[7];
              trsav->kern[9] = trsav->kern[7];
          
              if(wdtd) nn /= sum_wt;

              if(nn > trsav->datnm[7])trsav->datnm[7] = trsav->datnm[9] = nn;

              pst->stat6 = *den[0]; *den[0] = 0.0;
              pst->stat8 = *den[1]; *den[1] = 0.0;

              dtdn = 0;

            }

            else{

               printf("****WARNING****, no data available for Track density and Lifetime statistics.\n\n");


            }


/*            trsav->datnm[7] = trsav->datnm[9] = nlif; */

        }

        free(dtd);
        free(wdtd);

        trd = 0;


      }


      if(sqth){
         for(i=0; i< nelm; i++) free((lf + i)->lgrid);
         for(i=1; i<nlev; i++) free(sqth[i]);
         free(sqth);
         free(gvecs);
      }


/* free memory used for statistical calculation */

      free(den[2]);
      free(den[1]);
      free(den[0]);
      free(den);

      free(tele);


      times(&usage);
      tendall = usage.tms_utime + usage.tms_stime;


/* Set region back to degrees */

      trsav->xa1 = xr1;
      trsav->xa2 = xr2;
      trsav->ya1 = yr1;
      trsav->ya2 = yr2;
      
      printf("%f %f %f %f\n", (double)tstartall, (double)tstarttd, (double)tendall, (double)CLKTCK);

      printf("****INFORMATION****, cpu for all is %f; cpu for track density is %f\n\n",
             (double)(tendall - tstartall) / (double)CLKTCK, (double)(tendall - tstarttd) / (double)CLKTCK);

/* write statistic data to file */

      trsav->add_den_sc = 1.0;

      if(!(trsav->scden)){

         ph = open_file(statsf, "w");
         statdmp(ph, trsav);
         close_file(ph, statsf);

/* write statistics to NETCDF file */

         strncpy(stcdff, statsf, MAXCHR);
         strcat(stcdff, "_1.nc");
         write_gstat(trsav, gstat1, ioff, 1, 0);
         netcdf_write_stats(trsav, gstat1, stcdff, trfil, ioff, 0);
         if(gstat2){
            strncpy(stcdff, statsf, MAXCHR);
            strcat(stcdff, "_2.nc");
            write_gstat(trsav, gstat2, trsav->ptnum - ioff, 2, ioff);
            netcdf_write_stats(trsav, gstat2, stcdff, trfil, trsav->ptnum - ioff, ioff);
         }

         printf("Do you want to scale the pdf's to number densities? 'y' or 'n'\n");
         scanf("\n");
         if(getchar() == 'y'){

            trsav->scden = 1;
            trsav->add_den_sc = 1.0;

            printf("Do you want any additional scaling? 'y' or 'n'\n");
            scanf("\n");
            if(getchar() == 'y'){
              printf("What is the value of the additional (multiplicitive) scaling?\n");
              scanf("%f", &(trsav->add_den_sc));

            }

            sscale = trsav->add_den_sc;

            if(iper_num){

               if(sum_per < TOLIPER){

                  printf("****ERROR****, sum of periods is zero, no scaling to number density.\n\n");
                  exit(1);

               }

               for(i=0; i < trsav->ptnum; i++){

                   pst = trsav->ptst + i;
                   pst->stat3 *= ((double)(trsav->datnm[4]) * sscale / (double)sum_per);
                   pst->stat4 *= ((double)(trsav->datnm[5]) * sscale / (double)sum_per);
                   pst->stat5 *= ((double)(trsav->datnm[6]) * sscale / (double)sum_per);
                   pst->stat6 *= ((double)(trsav->datnm[7]) * sscale / (double)sum_per);

               }

            }

            else {

               for(i=0; i < trsav->ptnum; i++){

                   pst = trsav->ptst + i;
                   pst->stat3 *= (double)(trsav->datnm[4]) * sscale;
                   pst->stat4 *= (double)(trsav->datnm[5]) * sscale;
                   pst->stat5 *= (double)(trsav->datnm[6]) * sscale;
                   pst->stat6 *= (double)(trsav->datnm[7]) * sscale;

               }

            }

            ph = open_file(statss, "w");
            statdmp(ph, trsav);
            close_file(ph, statss);

            strncpy(stcdff, statss, MAXCHR);
            strcat(stcdff, "_1.nc");
            write_gstat(trsav, gstat1, ioff, 1, 0);
            netcdf_write_stats(trsav, gstat1, stcdff, trfil, ioff, 0);
            if(gstat2){
              strncpy(stcdff, statss, MAXCHR);
              strcat(stcdff, "_2.nc");
              write_gstat(trsav, gstat2, trsav->ptnum - ioff, 2, ioff);
              netcdf_write_stats(trsav, gstat2, stcdff, trfil, trsav->ptnum - ioff, ioff);
            }

         }

         else{
            printf("****WARNING*****, statistics are already scaled, no scaling  \r\n"
                   "                  is performed and statistics not written to \r\n"
                   "                  file.                                      \n\n");

         }

      }      

   }

   else if(ifmp == 2){

      printf("what is the name of the statistics file to be read\n");

      scanf("%s", charin);

      ph = open_file(charin, "r");
      trsav = read_stats(ph);
      close_file(ph, charin);

      if(!(trsav->scden)){

         printf("Do you want to scale the pdf's to number densities? 'y' or 'n'\n");
         scanf("\n");
         if(getchar() == 'y'){

            trsav->scden = 1;
            trsav->add_den_sc = 1.0;

            printf("Do you want any additional scaling? 'y' or 'n'\n");
            scanf("\n");
            if(getchar() == 'y'){
              printf("What is the value of the additional (multiplicitive) scaling?\n");
              scanf("%f", &(trsav->add_den_sc));

            }

            sscale = trsav->add_den_sc;

            if(iper_num){

               if(sum_per < TOLIPER){

                  printf("****ERROR****, sum of periods is zero, no scaling to number density.\n\n");
                  exit(1);

               }

               for(i=0; i < trsav->ptnum; i++){

                   pst = trsav->ptst + i;
                   pst->stat3 *= ((double)(trsav->datnm[4]) * sscale / (double)sum_per);
                   pst->stat4 *= ((double)(trsav->datnm[5]) * sscale / (double)sum_per);
                   pst->stat5 *= ((double)(trsav->datnm[6]) * sscale / (double)sum_per);
                   pst->stat6 *= ((double)(trsav->datnm[7]) * sscale / (double)sum_per);

               }

            }

            else {

               for(i=0; i < trsav->ptnum; i++){

                   pst = trsav->ptst + i;
                   pst->stat3 *= (double)(trsav->datnm[4]) * sscale;
                   pst->stat4 *= (double)(trsav->datnm[5]) * sscale;
                   pst->stat5 *= (double)(trsav->datnm[6]) * sscale;
                   pst->stat6 *= (double)(trsav->datnm[7]) * sscale;

               }

            }



            printf("Do you want to output the statistics with scaled densities to file? 'y' or 'n'\n");
            scanf("\n");
            if(getchar() == 'y'){
               ph = open_file(statss, "w");
               statdmp(ph, trsav);
               close_file(ph, statss);
            }

         }

      }

      ta = (trsav->xa2 - trsav->xa1) * (trsav->ya2 - trsav->ya1);

      printf("total area is %f, how many sample points are required for plotting if required\n", ta);
      scanf("%d", &nb);

      bl = (float)sqrt((double)(ta/nb));

      if(bl > (trsav->xa2 - trsav->xa1) || bl > (trsav->ya2 - trsav->ya1)){
         printf("***error***, number of sample points is incompatable  \r\n"
                "             with area, exiting from statistical routine \r\n"
                "%s\n\n", __FILE__);
         return;
      }

      nbx = (trsav->xa2 - trsav->xa1) / bl;
      nby = (trsav->ya2 - trsav->ya1) / bl;


   }

   if(ifmp){

      printf("do you wish to display the statistics, 'y' or 'n'\n");

      scanf("\n");
      if(getchar() == 'y'){

         ++dcount;

         printf("do you want individual plots, '0' \r\n"
                "          or a multiple plot, '1' \n\n");

         scanf("%d", &cpty);

         plot_stats(trsav, cpty, nbx, nby);

         printf("do you want to plot the initiation and dispersion points,\r\n"
                "'y' or 'n'\n");
         scanf("\n");
         if(getchar() == 'y') {

            tinit = (struct tot_tr * )malloc_initl(sizeof(struct tot_tr));
            mem_er((tinit == NULL) ? 0 : 1, sizeof(struct tot_tr));

            strncpy(charin, INITTRS, MAXCHR);
            if(iext) strcpy(strstr(charin, EXTENSION), fext);

            if(!fexist(charin, "r")){

               printf("***ERROR***, track initiation file does not exist,\r\n"
                      "             input an alternative filename or continue\n");

               scanf("%s", charin);

            }

            if(fexist(charin, "r")){

               gend = open_file(charin, "r");

               fscanf(gend, "%d", &(tinit->num));

               tinit->trpt = (struct fet_pt_tr * )calloc(tinit->num, sizeof(struct fet_pt_tr));
               mem_er((tinit->trpt == NULL) ? 0 : 1, tinit->num * sizeof(struct fet_pt_tr));

               for(i=0; i< tinit->num; i++){
                  at = tinit->trpt + i;
                  fscanf(gend, "%d %f %f %e %*f %*f", &(at->fr_id), &(at->xf), &(at->yf), &(at->zf));
               }

               close_file(gend, charin);

               splice_plot(tinit, tinit->num, 1, 1, 1, 's');

            }


            strncpy(charin, DISPTRS, MAXCHR);
            if(iext) strcpy(strstr(charin, EXTENSION), fext);

            if(!fexist(charin, "r")){

               printf("***ERROR***, track initiation file does not exist,\r\n"
                      "             input an alternative filename or continue\n");

               scanf("%s", charin);

            }

            if(fexist(charin, "r")){

               lysd = open_file(charin, "r");

               fscanf(gend, "%d", &(tinit->num));

               tinit->trpt = (struct fet_pt_tr * )realloc_n(tinit->trpt, tinit->num * sizeof(struct fet_pt_tr));
               mem_er((tinit->trpt == NULL) ? 0 : 1, tinit->num * sizeof(struct fet_pt_tr));

               for(i=0; i< tinit->num; i++){
                  at = tinit->trpt + i;
                  fscanf(gend, "%d %f %f %e %*f %*f", &(at->fr_id), &(at->xf), &(at->yf), &(at->zf));

               }

               close_file(lysd, charin);

               splice_plot(tinit, tinit->num, 1, 1, 1, 's');

            }

            if(tinit){

               free(tinit->trpt);
               free(tinit);

            }

         }

      }

   }

   if(gstat1 && gstat1 != gr) { free(gstat1->xgrid); free(gstat1->ygrid); free(gstat1);}
   if(gstat2) { free(gstat2->xgrid); free(gstat2->ygrid); free(gstat2);}

   free(stat);

/* free speed memory */

   if(phsp){             

      for(i=0; i < trphn; i++) free((phsp+i)->trpt);

      free(phsp);

   }

/* free statistic memory */

   if(trsav){

     free(trsav->ptst);

     free(trsav);

   }


   return;

}


/* don't like to do this but it saves on code repitition */


/*#############################################################################

                 SUPPORT ROUTINES FOR statistic()

  ############################################################################# */


void select_region(float *xr1, float *xr2, float *yr1, float *yr2)

{

      int ny;

      printf("Specify a region for the calculation of the statistics,\r\n"
             "the current region is defined by\r\n"
             "X= %f to %f and Y= %f to %f\n", xmn, xmx, ymn, ymx);

      ny = 1;

      while(ny){

         printf("Input X-range 'x1', 'x2'\n");
         scanf("%f %f", xr1, xr2);

         if(*xr1 < xmn || *xr2 > xmx) 

           printf("***error***, invalid region on input, try again\n");

         else ny = 0;

      }

      ny = 1;

      while(ny){

         printf("Input Y-range 'y1', 'y2'\n");
         scanf("%f %f", yr1, yr2);

         if(*yr1 < ymn || *yr2 > ymx) 

           printf("***error***, invalid region on input, try again\n");

         else ny = 0;

      }

      printf("Chosen area is                                     \r\n\n"
             "X= %f to %f                                          \r\n"
             "Y= %f to %f                                            \n",
              *xr1, *xr2, *yr1, *yr2);

    return;

}


struct pt_stat *assign_sample_pt(struct tot_stat *ttst, float x, float y, float yy1, float yy2)

{

    int bn=0;

    struct pt_stat *pst; 

    if(ttst->ptst == NULL){

       ttst->ptst = (struct pt_stat * )malloc_initl(sizeof(struct pt_stat));
       mem_er((ttst->ptst == NULL) ? 0 : 1, sizeof(struct pt_stat));

    }

    else {

       bn = (ttst->ptnum + 1)*sizeof(struct pt_stat);

       ttst->ptst = (struct pt_stat * )realloc_n(ttst->ptst, bn);
       mem_er((ttst->ptst == NULL) ? 0 : 1, bn);

    }

    pst = (ttst->ptst + ttst->ptnum);
    ++(ttst->ptnum);
    pst->xs = x;
    pst->ys = y;

    return ttst->ptst;

}

extern int tom;

struct dpt *extract_data(struct fet_pt_tr *at, struct dpt *dt, int *dn, float xr1, float xr2, float yr1, float yr2, int fld, float *wt, int imem)

{


    float px, py;

    struct dpt *dd=NULL;

    px = (xr2 - at->xf) * (at->xf - xr1);
    py = (yr2 - at->yf) * (at->yf - yr1);

    *wt = at->wght;

    if(px >= 0. && py >= 0.) {

       ++(*dn);

       if(!imem){

          if(dt == NULL) {

             dt = (struct dpt * )malloc_initl(sizeof(struct dpt));
             mem_er((dt == NULL) ? 0 : 1, sizeof(struct dpt));

          }
  
          else {

             dt = (struct dpt *) realloc_n(dt, (*dn) * sizeof(struct dpt));
             mem_er((dt == NULL) ? 0 : 1, (*dn) * sizeof(struct dpt));

          }

       }

       dd = dt + (*dn) - 1;

       dd->xdt = at->pp[0];
       dd->ydt = at->pp[1];
       dd->zdt = at->pp[2];

       if(fld == 'g') dd->sdt = at->gwthr;

       else if(fld == 'a') {

          if(at->sh_an >= 0.){
             dd->sdt = at->sh_an;
             dd->vec[0] = at->or_vec[0];
             dd->vec[1] = at->or_vec[1];
          }

       }

       else if(fld == 't') dd->sdt = at->tend;

       else if(fld == 'r') {
          if(at->area >= 0.) dd->sdt = at->area;
       }
       else dd->sdt = at->zf;

       if(fld == 'v'){dd->vec[0] = at->vec[0]; dd->vec[1] = at->vec[1];}

    }

    return dt;

}



double arc(struct fet_pt_tr *at, struct dpt *pst)

{

    double dd;

    dd = at->pp[0] * pst->xdt + at->pp[1] * pst->ydt + at->pp[2] * pst->zdt;
 
    if(fabs(dd) > 1.) dd = (dd < 0.) ? -1.0 : 1.0;

    return dd;

}


float *weights(float *wght, int np, float wt, int imem)

{

    if(!imem){

       if(!wght){
         if(np != 1){
            printf("****ERROR****, inconsistency between weights and track points.\r\n"
                   "               weight pointer is NULL but current number of   \r\n"
                   "               points is not 1.                               \n\n");
            exit(1);

         }
    
         wght = (float *)malloc_initl(sizeof(float));
         mem_er((wght == NULL) ? 0 : 1, sizeof(float));

       }

       else { 

         wght = (float *) realloc_n(wght, np * sizeof(float));
         mem_er((wght == NULL) ? 0 : 1, np * sizeof(float));

       }

    }

    *(wght + np - 1) = wt;

    return wght;

}

void available_stats(void)

{

    printf("Input the Id. of the statistics **NOT** required.        \n\n"
           " '1',  mean and variance of strength and feature density \r\n"
           " '2',  mean and variance of speed and mean velocity      \r\n"
           " '3',  genesis                                           \r\n"
           " '4',  lysis                                             \r\n"
           " '5',  growth/decay rate                                 \r\n"
           " '6',  feature anistropy and orientation                 \r\n"
           " '7',  tendency                                          \r\n"
           " '8',  mean area (currently not active).                 \r\n"
           "Input '0' to finish selection. Statistics not selected   \r\n"
           "will be computed.                                        \n\n");

    return;

}

/* compute area of a spherical triangle using only side information */

double area_sp_triangle(SQT *ss, VEC *gv, int ty)
{

    double ar1, ar2, ar3;
    double sum, area;

    VEC *v1=NULL, *v2=NULL, *v3=NULL;

    if(!ty){

       v1 = gv + ss->ivec[0];
       v2 = gv + ss->ivec[1];
       v3 = gv + ss->ivec[2];

    }

    else {

       v1 = gv + ((LEAF *)ss)->ivec[0];
       v2 = gv + ((LEAF *)ss)->ivec[1];
       v3 = gv + ((LEAF *)ss)->ivec[2];
    }

    if(fabs(dotp(v1, v1) - 1.0) > TOLVEC || 
       fabs(dotp(v2, v2) - 1.0) > TOLVEC || 
       fabs(dotp(v3, v3) - 1.0) > TOLVEC     )
            printf("****WARNING****, vectors not unit length, area calculation not reliable.\n\n");

/* determine arc lengths of triangle sides */


    ar1 = acos(dotp(v1, v2));
    ar2 = acos(dotp(v2, v3));
    ar3 = acos(dotp(v3, v1));

    sum = (ar1 + ar2 + ar3) * 0.5;

    area = 4.0 * atan(sqrt(tan(0.5 * sum) * tan(0.5 * (sum - ar1)) * tan(0.5 * (sum - ar2)) * tan(0.5 * (sum - ar3))));

    return area;

}

void write_gstat(struct tot_stat *stat, GRID *gstat, int snum, int gid, int off)
{

   int i;

   struct pt_stat *pst=NULL;

   static char *gf[] = {"gstat1.dat", "gstat2.dat"};
   char *gff=NULL;

   FILE *fg=NULL;

   if(gid < 1 || gid > 2){
      printf("****ERROR****, statistics grid ID not known.\n\n");
      exit(1);
   }

   if(gid == 1) gff = gf[0];
   else if(gid == 2) gff = gf[1];

   pst = stat->ptst + off;

   fg = open_file(gff, "w");

   fprintf(fg, "%s\n", trfil);
   fprintf(fg, "%d\n", snum);
   fprintf(fg, "%d %d\n", gstat->ix, gstat->iy);
   fprintf(fg, "%d %d\n", gstat->prgr, gstat->prty);
   fprintf(fg, "%f %f\n", gstat->alng, gstat->alat);
   fprintf(fg, "%s\n", gstat->prgr_nm);
   fprintf(fg, "%s\n", gstat->prty_nm);
   fwrite(gstat->xgrid, gstat->ix * sizeof(float), 1, fg);
   fprintf(fg, "\n");
   fwrite(gstat->ygrid, gstat->iy * sizeof(float), 1, fg);
   fprintf(fg, "\n");
   for(i=0; i < snum; i++) fprintf(fg, "%d ", (pst + i)->apos);
   fprintf(fg, "\n");


   close_file(fg, gff);

   return;

}
