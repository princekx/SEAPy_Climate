#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include <sys/types.h>
#include "grid.h"
#include "file_handle.h"
#include "file_cat_out.h"
#include "mem_er.h"
#include "reg_dat.h"
#include "sphery_dat.h"
#include "bisp.h"
#include "netcdf_info.h"
#include "m_values.h"
#include "geo_values.h"

/* function to compute vorticty, speed or EKE */

#ifdef  REGULAR

#ifdef  NOUNDERSCORE

void bisp(double * , double * , double * , int * , double * , int * , 
          double * , int * , double * , double * , int * );
	  
#else

void bisp_(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#endif

#endif


void surfit(double * , int , int , ... );
int query_interp(int );
int read_field(float * ,float * , float , FILE * , int , int , int , int , int );
void write_header(GRID * , FILE * );
int missing(float, float, int );

extern GRID *gr;
extern int form;
extern int tl, gof;
extern int frnum;
extern float *ap;
extern float period;
extern int tom;
extern int x1u, x2u, y1u, y2u;
extern float xmn, ymn, xmx, ymx;

void compute_vorticity(FILE *fst, off_t pl, int iper)
{

   int i, j, ii;
   int ir=0;
   int nfu=0, nfv=0;
   int smty=0;
   int fruser=0, rf=0;
   int ofil=0;
   int itpfc=0;
   int maxe=0;
   int idif=0;
   int dim=gr->ix*gr->iy, dimn=0;
   int nwr=0;
   int ictype=0;
   int imiss='n';
   int icmp=0;
   int addvor='n';

   int interpy=0;

   int fs[2], fr[2], fe[2];

   off_t uplace1=0, vplace1=0, uplace2=0, vplace2=0;
   off_t place1=0, place2=0, chrnum=0;
   off_t uchrnum=0, vchrnum=0;

   float *uap=NULL, *vap=NULL, *vor=NULL;
   float *mgeo=NULL;
   float *nuvp=NULL;
   float *lat2d=NULL;
   float u1, v1, spp;
   float mval=0.0, mdval=0.0;
   float vorcs=1.0;

   double sm=0.0;
   double udd[4], vdd[4];
   double xx, yy, uz, vz, gz;
   double *coslt=NULL, *sinlt=NULL;
   double *tanr=NULL;
   double coslat=0.0, sinlat=0.0, cc=0.0, sc=0.0;
   double fpp=FP_PI*FP_PI, ff=0.0, fft=OMEGA*EARTH_RADIUS_M*EARTH_RADIUS_M;

   FILE *uu=NULL, *vv=NULL;
   FILE *fgeo=NULL;
   FILE *fvor=NULL;
   FILE *latf=NULL;

   GRID *newg=NULL, *gtmp=NULL;

   char ufil[MAXCHR], vfil[MAXCHR], latfil[MAXCHR];
   char gfil[MAXCHR];
   char vorf[MAXCHR];

   struct sp_dat *coknu=NULL, *coknv=NULL;      /* spline data for any of the methods   */
   struct rspline *rsu=NULL, *rsv=NULL;         /* data structure for use with smoopy   */
   struct sspline *ssu=NULL, *ssv=NULL;         /* data structure for use with sphery   */
   struct savedat *sdu=NULL, *sdv=NULL;         /* data structure for use with smoopy   */
   struct savedat_sphy *ssdu=NULL, *ssdv=NULL;  /* data structure for use with sphery   */

   ap = NULL;
   x1u = 1;
   x2u = gr->ix;
   y1u = 1;
   y2u = gr->iy;

   xmn = *(gr->xgrid);
   xmx = *(gr->xgrid + gr->ix - 1);
   ymn = *(gr->ygrid);
   ymx = *(gr->ygrid + gr->iy - 1);

   if(gr->prty) {
      printf("****ERROR****, data should not be transformed to a projection other than LAT-LONG.\n\n");
      exit(1);
   }

   printf("****INFORMATION****, input fields should be U and V,           \r\n"
          "                     unless Geostrophic Vorticty is required   \r\n"
          "                     when only MSLP or GEOPOTENTIAL field are  \r\n"
          "                     needed for input.                         \n\n");

   printf("Do you want speed,     '0'             \r\n"
          "Kinetic Energy,        '1'             \r\n"
          "Vorticty from winds,   '2'             \r\n"
          "Geostrophic Vorticty,  '3'             \r\n"
          "KE from wind speed,    '4' calculation?\n\n");
   scanf("%d", &ictype);

   if(ictype < 0 || ictype > 4) {
      printf("****ERROR****, unknown option %d\n\n", ictype);
      exit(1);
   }


   if(!ictype){
      printf("****INFORMATION****, computing speeds.\n\n");
   }
   else if(ictype == 1){
      printf("****INFORMATION****, computing eddy kinetic energy.\n\n");
   }
   else if(ictype == 2){
      printf("****INFORMATION****, computing vorticity from winds.\n\n");
      idif = 1;
   }
   else if(ictype == 3){
      printf("****INFORMATION****, computing geostrophic vorticity.\n\n");
      if(tom != 'g'){
        printf("****WARNING****, no geostrophic vorticty for non-spherical data.\n\n");
        exit(1);
      }
      idif = 2;
   }
   else if(ictype == 4){
      printf("****INFORMATION****, computing KE from wind speed.\n\n");
   }

   if((ictype == 2 || ictype == 3) && tom == 'g') {
      if((90.0 + *(gr->ygrid)) <= TOLPOLE || (90.0 - *(gr->ygrid + gr->iy -1)) <= TOLPOLE) {
         printf("****ERROR****, data contains north or south poles, remove poles \r\n"
                "               using initial options before calculation.        \n\n");
         exit(1);
      }
      coslt = (double *)calloc(gr->iy, sizeof(double));
      mem_er((coslt == NULL) ? 0 : 1, (gr->iy)*sizeof(double));
      sinlt = (double *)calloc(gr->iy, sizeof(double));
      mem_er((sinlt == NULL) ? 0 : 1, (gr->iy)*sizeof(double));
      for(i=0; i < gr->iy; i++) sincos((double)(*(gr->ygrid + i)) * FP_PI, sinlt + i, coslt + i);
   }

   if(ictype < 3){

      printf("Does current file have both U and V components? 'y' or 'n'\n\n");
      scanf("\n");
      if(getchar() == 'y'){
        ofil = 1;

        if(form != 4) {
           uu = fst;
           vv = fst;
           fseeko(fst, pl, FSTART);
           uplace1 = vplace1 = pl;	  
        }
        else{
           printf("Choose U field.\n\n");
           uu = (FILE *)nc_clone((NETCDF_INFO *)fst, NULL, NC_SAME);
           ((NETCDF_INFO *)uu)->iframe = 0;
           printf("Choose V field.\n\n");
           vv = (FILE *)nc_clone((NETCDF_INFO *)fst, NULL, NC_SAME);
           ((NETCDF_INFO *)vv)->iframe = 0;
        }

      }
      else {

        printf("What file contains the U field?\n\n");
        scanf("%s", ufil);
        printf("What file contains the V field?\n\n");
        scanf("%s", vfil);

        if(form != 4) {
           uu = open_file(ufil, "r");
           fseeko(uu, pl, FSTART);
           vv = open_file(vfil, "r");
           fseeko(vv, pl, FSTART);
           uplace1 = vplace1 = pl;	  
        }
        else{
           printf("Choose U field.\n\n");
           uu = (FILE *)nc_clone((NETCDF_INFO *)fst, ufil, NC_OPEN_MODE);
           ((NETCDF_INFO *)uu)->iframe = 0;
           printf("Choose V field.\n\n");
           vv = (FILE *)nc_clone((NETCDF_INFO *)fst, vfil, NC_OPEN_MODE);
           ((NETCDF_INFO *)vv)->iframe = 0;
        }

      }

      printf("What are the start, rate and end frame Id's for the U field?\n\n");
      scanf("%d %d %d", fs, fr, fe);
      printf("What are the start, rate and end frame Id's for the V field?\n\n");
      scanf("%d %d %d", fs+1, fr+1, fe+1);
      maxe = (fe[0] > fe[1]) ? fe[0]: fe[1];
      if(fr[0] != fr[1]){
         printf("****WARNING****, frame rates are not identical.\b\n\n");
      }
      
      if(ictype == 2 && tom == 'e'){
         printf("What scaling is required?\n\n");
	 scanf("%f", &vorcs);

         printf("Do you want to add the additional vorticity term for limited areas, 'y' for yes and 'n' for no.\n\n");
	 scanf("\n");
	 addvor = getchar();
	 if(addvor == 'y'){
	 
	    lat2d = (float *)calloc(dim, sizeof(float));
            mem_er((lat2d == NULL) ? 0: 1, dim * sizeof(float));
	    tanr = (double *)calloc(dim, sizeof(double));
            mem_er((tanr == NULL) ? 0: 1, dim * sizeof(double));
	 
	    printf("Which file contains the grid of latitudes?\n\n");
	    scanf("%s", latfil);
            if(form != 4) {
               latf = open_file(latfil, "r");
               fseeko(latf, pl, FSTART);	    
	    }
	    else {
               latf = (FILE *)nc_clone((NETCDF_INFO *)fst, latfil, NC_OPEN_MODE);
               ((NETCDF_INFO *)latf)->iframe = 0;	    
	    }    
	 }
	 
	 read_field(lat2d, NULL, 1.0, latf, 1, 'n', 'n', '0', 'n');
	 
         for(i=0; i < dim; i++){
	     *(tanr + i) = tan(*(lat2d + i) * FP_PI)/ EARTH_RADIUS_M;
	 }
	 
         if(form != 4) close_file(latf, latfil);
         else netcdf_close((NETCDF_INFO *)latf); 
	 free(lat2d);
      }

   }
   else {
      if(ictype == 3)
         printf("Does current file contain the MSLP, GEOPOTENTIAL? 'y' or 'n'\n\n");
      else if(ictype == 4)
         printf("Does current file contain the WIND SPEED? 'y' or 'n'\n\n");

      scanf("\n");
      if(getchar() == 'y'){
        ofil = 1;
        if(form != 4) {
           fgeo = fst;
           fseeko(fst, pl, FSTART);
           place1 = pl;	  
        }
        else{
           if(ictype == 3)
              printf("Choose MSLP or Geopotential field.\n\n");
           else if(ictype == 4)
              printf("Choose the WIND SPEED field.\n\n");
           fgeo = (FILE *)nc_clone((NETCDF_INFO *)fst, NULL, NC_SAME);
           ((NETCDF_INFO *)fgeo)->iframe = 0;
        }

      }
      else {
        if(ictype == 3)
           printf("What file contains the MSLP or Geopotential field?\n\n");
        else if(ictype == 4)
           printf("What file contains the WIND SPEED field?\n\n");

        scanf("%s", gfil);

        if(form != 4) {
           fgeo = open_file(gfil, "r");
           fseeko(fgeo, pl, FSTART);
           place1 = pl;	  
        }
        else{
           if(ictype == 3)
              printf("Choose MSLP or Geopotential field.\n\n");
           else if(ictype == 4)
              printf("Choose the WIND SPEED field.\n\n");

           fgeo = (FILE *)nc_clone((NETCDF_INFO *)fst, gfil, NC_OPEN_MODE);
           ((NETCDF_INFO *)fgeo)->iframe = 0;
        }

      }

      if(ictype == 3){
         printf("Is field MSLP or geopotential, input 'p' or 'g'\n\n");
         scanf("\n");
         if(getchar() == 'p') fft *= RHO;
      }

      printf("What are the start, rate and end frame Id's for the MSLP or Geopotential field?\n\n");
      scanf("%d %d %d", fs, fr, fe);
      maxe = fe[0];
      fs[1] = fs[0]; fr[1] = fr[0]; fe[1] = fe[0];
   }

   if(ictype == 0 || ictype == 1 || ictype == 4){
      if(ictype == 0 || ictype == 1){
         if(form == 4){
            if(((NETCDF_INFO *)uu)->imiss || ((NETCDF_INFO *)vv)->imiss) {
               printf("****WARNING****, field has missing data value or fill value.\n\n");
               imiss = 'y';
            }
         }
         else{
            printf("Are there missing data values, delineated by some extreme value. 'y' or 'n'\n\n");
            scanf("\n");
            imiss=getchar();
         }
      }
      else{
         if(form == 4){
            if(((NETCDF_INFO *)fgeo)->imiss) {
               printf("****WARNING****, field has missing data value or fill value.\n\n");
               imiss = 'y';
            }
         }
         else{
            printf("Are there missing data values, delineated by some extreme value. 'y' or 'n'\n\n");
            scanf("\n");
            imiss=getchar();
         }
      }
      if(imiss == 'y'){
         printf("What is the missing data value to test against?\n\n");
         scanf("%f", &mval);

         printf("How do you want to compare with missing value?\r\n"
                "'0' -- equality.                              \r\n"
                "'1' -- less than.                             \r\n"
                "'2' -- greater than.                          \n\n");
         scanf("%d", &icmp);

         printf("What is the default missing data value which will be used as the missing or fill value?\n\n");
         scanf("%f", &mdval);

      }
   }

   else {
      if(ictype == 2 && form == 4){
         if(((NETCDF_INFO *)uu)->imiss || ((NETCDF_INFO *)vv)->imiss) {
            printf("****WARNING****, field has missing data value or fill value, do you want to continue, 'y' or 'n'\n\n");
            scanf("\n");
            if(getchar() == 'n') exit(1);
         }
      }
      else if(ictype == 3 && form == 4){
         if(((NETCDF_INFO *)fgeo)->imiss) {
            printf("****WARNING****, field has missing data value or fill value, do you want to continue, 'y' or 'n'\n\n");
            scanf("\n");
            if(getchar() == 'n') exit(1);
         }
      }
      else {
         printf("****WARNING****, calculation assumes there are no missing data values.\n\n");
      }

   }

   printf("What is the output file name for the speed, vorticty or KE?\n\n");
   scanf("%s", vorf);

   if(ictype < 3){

      uap = (float *)calloc(dim, sizeof(float));
      mem_er((uap == NULL) ? 0: 1, dim * sizeof(float));  

      vap = (float *)calloc(dim, sizeof(float));
      mem_er((vap == NULL) ? 0: 1, dim * sizeof(float)); 

   }

   else {
      mgeo = (float *)calloc(dim, sizeof(float));
      mem_er((mgeo == NULL) ? 0: 1, dim * sizeof(float));
   }

   vor=(float *)calloc(dim, sizeof(float));
   mem_er((vor == NULL) ? 0: 1, dim * sizeof(float));    

   if(ictype == 2 || ictype == 3){

      smty = query_interp(0);

      coknu = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
      mem_er((coknu == NULL) ? 0 : 1, sizeof(struct sp_dat));

      if(ictype == 2){
         coknv = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
         mem_er((coknv == NULL) ? 0 : 1, sizeof(struct sp_dat));
      }

      if(!smty){
         printf("How mauch wrap round is required for fitting 0 <= w <= %d\n\n", gr->ix);
         scanf("%d", &nwr);
         if(iper == 'n'){
            printf("****WARNING****, wrapround may not be valid for this data set.\n\n");
         }
         if(nwr < 0 || nwr > gr->ix){
           printf("****ERROR****, invalid wrapround.\n\n");
           exit(1);
         }
      }

      if(smty){

         ssu = (struct sspline * )malloc_initl(sizeof(struct sspline));
         mem_er((ssu == NULL) ? 0 : 1, sizeof(struct sspline));

         ssdu = (struct savedat_sphy * )malloc_initl(sizeof(struct savedat_sphy));
         mem_er((ssdu == NULL) ? 0 : 1, sizeof(struct savedat_sphy));

         if(ictype == 2){

            ssv = (struct sspline * )malloc_initl(sizeof(struct sspline));
            mem_er((ssv == NULL) ? 0 : 1, sizeof(struct sspline));

            ssdv = (struct savedat_sphy * )malloc_initl(sizeof(struct savedat_sphy));
            mem_er((ssdv == NULL) ? 0 : 1, sizeof(struct savedat_sphy));

         }

      }

      else{

         rsu = (struct rspline * )malloc_initl(sizeof(struct rspline));
         mem_er((rsu == NULL) ? 0 : 1, sizeof(struct rspline));

         sdu = (struct savedat * )malloc_initl(sizeof(struct savedat));
         mem_er((sdu == NULL) ? 0 : 1, sizeof(struct savedat));

         if(ictype == 2){

            sdv = (struct savedat * )malloc_initl(sizeof(struct savedat));
            mem_er((sdv == NULL) ? 0 : 1, sizeof(struct savedat));

            rsv = (struct rspline * )malloc_initl(sizeof(struct rspline));
            mem_er((rsv == NULL) ? 0 : 1, sizeof(struct rspline));

         }

         newg = (GRID *)malloc_initl(sizeof(GRID));
         mem_er((newg == NULL) ? 0 : 1, sizeof(GRID));
         newg->ix = gr->ix + 2 * nwr;
         newg->iy = gr->iy;
         dimn = newg->ix * newg->iy;

         newg->xgrid = (float * )calloc(newg->ix, sizeof(float));
         mem_er((newg->xgrid == NULL) ? 0 :1, newg->ix * sizeof(float));

         memcpy(newg->xgrid + nwr, gr->xgrid, gr->ix * sizeof(float));

         if(iper == 'y'){

            for(i=0; i < nwr; i++) {
               *(newg->xgrid + i) = *(gr->xgrid + gr->ix - 1 - nwr + i) - period;
               *(newg->xgrid + newg->ix - nwr + i) = *(gr->xgrid + 1 + i) + period;
            }

         }

         else{

            for(i=0; i < nwr; i++){
               *(newg->xgrid + i) = *(gr->xgrid + gr->ix - nwr + i) - period;
               *(newg->xgrid + newg->ix - nwr + i) = *(gr->xgrid + i) + period;
            } 

         }

         newg->ygrid = gr->ygrid;

         nuvp=(float *)calloc(dimn, sizeof(float));
         mem_er((nuvp == NULL) ? 0: 1, dimn * sizeof(float));  

         x1u = 1;
         x2u = newg->ix;
         xmn = *(newg->xgrid);
         xmx = *(newg->xgrid + newg->ix - 1);

      }

   }


/* open file for writing the vorticity field */

   fvor = open_file(vorf, "w");

   write_header(gr, fvor);

/* nfu and nfv are the same for ictype=3 */

   while(nfu <= maxe && nfv <= maxe && !ir){


       if(ictype < 3) {
  
          if(!nfu && !nfv){

/* read U */

             if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir= 1;
       
             if(form != 4) {uplace2 = ftello(uu); uchrnum = uplace2 - uplace1;}
       
             if(fs[0] > 1){
               if(form != 4) fseeko(uu,(fs[0]-2)*uchrnum, ORIGIN);
	       else ((NETCDF_INFO *)uu)->iframe = fs[0] - 1;	  
	       if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir = 1;
               if(form != 4) uplace2 = ftello(uu);
             }
	    
             nfu = fs[0];

             if(ofil && form != 4) fseeko(uu, uplace1, FSTART);

/* read V */

             if((rf = read_field(vap, NULL, 1.0, vv, 1, 'n', 'n', '0', 'y'))) ir = 1;
          
             if(form != 4) {vplace2 = ftello(vv); vchrnum = vplace2 - vplace1;}

             if(uchrnum != vchrnum){
                printf("****ERROR****, field block sizes do not match for U and V. File %s, line %d\n\n", __FILE__, __LINE__);
                exit(1);
             }
       
             if(fs[1] > 1){
               if(form != 4) fseeko(vv,(fs[1]-2)*vchrnum, ORIGIN);
	       else ((NETCDF_INFO *)vv)->iframe = fs[1] - 1;	  
	       if((rf = read_field(vap, NULL, 1.0, vv, 1, 'n', 'n', '0', 'y'))) ir = 1;
               if(form != 4) vplace2 = ftello(uu);
             }
	    
             nfv = fs[1];          

             fruser = 1;
	   }
           else {

             if(ofil){
                if(form != 4) fseeko(uu, uplace2 + (fr[0]-1)*uchrnum, FSTART);
                else ((NETCDF_INFO *)uu)->iframe = nfu - 1;
                if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir = 1;
                if(form != 4) uplace2 = ftello(uu);

                if(form != 4) fseeko(vv, vplace2 + (fr[1]-1)*vchrnum, FSTART);
                else ((NETCDF_INFO *)vv)->iframe = nfv - 1;
                if((rf = read_field(vap, NULL, 1.0, vv, 1, 'n', 'n', '0', 'y'))) ir = 1;
                if(form != 4) vplace2 = ftello(vv);

             }

             else {

                if(form != 4) fseeko(uu, (fr[0]-1)*uchrnum, ORIGIN);
                else ((NETCDF_INFO *)uu)->iframe = nfu - 1;
                if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir = 1; 

                if(form != 4) fseeko(vv, (fr[1]-1)*vchrnum, ORIGIN);
                else ((NETCDF_INFO *)vv)->iframe = nfv - 1;
                if((rf = read_field(vap, NULL, 1.0, vv, 1, 'n', 'n', '0', 'y'))) ir = 1;
             }
          
             if(!ir) ++fruser;

           }

        }

        else {

           if(!nfu){

             if((rf = read_field(mgeo, NULL, 1.0, fgeo, 1, 'n', 'n', '0', 'y'))) ir= 1;
       
             if(form != 4) {place2 = ftello(fgeo); chrnum = place2 - place1;}
       
             if(fs[0] > 1){
               if(form != 4) fseeko(fgeo,(fs[0]-2)*chrnum, ORIGIN);
	       else ((NETCDF_INFO *)fgeo)->iframe = fs[0] - 1;	  
	       if((rf = read_field(mgeo, NULL, 1.0, fgeo, 1, 'n', 'n', '0', 'y'))) ir = 1;
               if(form != 4) place2 = ftello(fgeo);
             }
	    
             nfu = nfv = fs[0];
             fruser = 1;

           }

           else {
             if(ofil){
                if(form != 4) fseeko(fgeo, place2 + (fr[0]-1)*chrnum, FSTART);
                else ((NETCDF_INFO *)fgeo)->iframe = nfu - 1;
                if((rf = read_field(mgeo, NULL, 1.0, fgeo, 1, 'n', 'n', '0', 'y'))) ir = 1;
                if(form != 4) place2 = ftello(fgeo);

             }

             else {

                if(form != 4) fseeko(fgeo, (fr[0]-1)*chrnum, ORIGIN);
                else ((NETCDF_INFO *)fgeo)->iframe = nfu - 1;
                if((rf = read_field(mgeo, NULL, 1.0, fgeo, 1, 'n', 'n', '0', 'y'))) ir = 1; 

             }
          
             if(!ir) ++fruser;

           }

        }

        if(ir) continue;

        nfu += fr[0];
        nfv += fr[1];

        if(ictype == 2) {

#ifdef  REGULAR

           interpy = 1;

           if(tom == 'g'){
              for(i=0; i < gr->iy; i++){
                  coslat = *(coslt + i);
                  for(j=0; j < gr->ix; j++){
                      *(uap + i * gr->ix + j) *= coslat;
                  }
              }
           }

/* fit U  and V */

           if(!smty) {

              ap = nuvp;
              if(iper == 'y'){
                 for(i=0; i < gr->iy; i++){
                     memcpy(nuvp + i * newg->ix + nwr, uap + i * gr->ix, gr->ix*sizeof(float));
                     for(j=0; j < nwr; j++){
                        *(nuvp + i * newg->ix + j) = *(uap + i * gr->ix + gr->ix - 1 - nwr + j);
                        *(nuvp + i * newg->ix + newg->ix - nwr + j) = *(uap + i * gr->ix + 1 + j);
                     }
                 }
              }
              else {
                 for(i=0; i < gr->iy; i++){
                     memcpy(nuvp + i * newg->ix + nwr, uap + i * gr->ix, gr->ix*sizeof(float));

                     for(j=0; j < nwr; j++){
                        *(nuvp + i * newg->ix + j) = *(uap + i * gr->ix + gr->ix - nwr + j);
                        *(nuvp + i * newg->ix + newg->ix - nwr + j) = *(uap + i * gr->ix + j);
                     }
                 }
              }

              gtmp = gr;
              gr = newg;

              surfit(&sm, itpfc, smty, coknu, rsu, sdu);

              gr = gtmp;

              if(iper == 'y'){

                 for(i=0; i < gr->iy; i++){
                     memcpy(nuvp + i * newg->ix + nwr, vap + i * gr->ix, gr->ix*sizeof(float));
                     for(j=0; j < nwr; j++){
                        *(nuvp + i * newg->ix + j) = *(vap + i * gr->ix + gr->ix - 1 - nwr + j);
                        *(nuvp + i * newg->ix + newg->ix - nwr + j) = *(vap + i * gr->ix + 1 + j);
                     }

                 }
              }
              else{
                 for(i=0; i < gr->iy; i++){
                     memcpy(nuvp + i * newg->ix + nwr, vap + i * gr->ix, gr->ix*sizeof(float));
                     for(j=0; j < nwr; j++){
                        *(nuvp + i * newg->ix + j) = *(vap + i * gr->ix + gr->ix - nwr + j);
                        *(nuvp + i * newg->ix + newg->ix - nwr + j) = *(vap + i * gr->ix + j);
                     }
                 }
              }

              gtmp = gr;
              gr = newg;

              surfit(&sm, itpfc, smty, coknv, rsv, sdv);

              gr = gtmp;

           }

           else {
              ap = uap;

	      surfit(&sm, itpfc, smty, coknu, ssu, ssdu);

              ap = vap;

              surfit(&sm, itpfc, smty, coknv, ssv, ssdv);

           }


           itpfc = 1;

           for(i=0; i < gr->iy; i++){

               if(smty){
                  xx = *(gr->ygrid + i);
	          xx = FP_PI2 - xx * FP_PI;
                  if(xx < 0.) xx = 0.0;
	          else if(xx > FPI) xx = FPI;
               }

               else yy = *(gr->ygrid + i);

               if(tom == 'g') coslat = *(coslt + i);

               for(j=0; j < gr->ix; j++){

                  if(smty){
		     yy = *(gr->xgrid + j); 
		     yy *= FP_PI;	 
	           }
		   else {
		     xx = *(gr->xgrid + j);
	           }

#ifdef  NOUNDERSCORE

                   bisp(&uz, udd, coknu->tx, &coknu->nx, coknu->ty, &coknu->ny, coknu->c, &coknu->ncof, &xx, &yy, &idif);
                   bisp(&vz, vdd, coknv->tx, &coknv->nx, coknv->ty, &coknv->ny, coknv->c, &coknv->ncof, &xx, &yy, &idif);

#else

                   bisp_(&uz, udd, coknu->tx, &coknu->nx, coknu->ty, &coknu->ny, coknu->c, &coknu->ncof, &xx, &yy, &idif);
                   bisp_(&vz, vdd, coknv->tx, &coknv->nx, coknv->ty, &coknv->ny, coknv->c, &coknv->ncof, &xx, &yy, &idif);

#endif
                   if(smty){
                     if(tom == 'g') {
                       *(vor + i * gr->ix + j) = (vdd[1] + udd[0]) / (coslat * EARTH_RADIUS_M);
		     }
                     else { /* this one does not really make sense for chosen options */
                       *(vor + i * gr->ix + j) = vorcs * (vdd[1] + udd[0]); 
		       if(addvor == 'y') {ii = i * gr->ix + j; *(vor + ii) += *(uap + ii) * *(tanr + ii);}
		     }
                  
                   }
                   else {
                     if(tom == 'g')
                       *(vor + i * gr->ix + j) = (vdd[0] - udd[1]) / (coslat * EARTH_RADIUS_M * FP_PI);
                     else { 
                       *(vor + i * gr->ix + j) = vorcs * (vdd[0] - udd[1]);
		       if(addvor == 'y') {ii = i * gr->ix + j; *(vor + ii) += *(uap + ii) * *(tanr + ii);}
                     }
                   }
		   

               }
	       
           }   

           if(iper == 'y') {
              for(i=0; i < gr->iy; i++) *(vor + i * gr->ix + gr->ix -1) = *(vor + i * gr->ix);
           }
	   

#else

           printf("****WARNING****, no surface fitting or calculation of vorticity performed.\r\n"
                  "                 May use finite difference in the future as alternative.  \n\n");

#endif


        }

        else if (ictype == 3){

#ifdef  REGULAR

           if(!smty) {

              ap = nuvp;
              if(iper == 'y'){
                 for(i=0; i < gr->iy; i++){
                     memcpy(nuvp + i * newg->ix + nwr, mgeo + i * gr->ix, gr->ix*sizeof(float));
                     for(j=0; j < nwr; j++){
                        *(nuvp + i * newg->ix + j) = *(mgeo + i * gr->ix + gr->ix - 1 - nwr + j);
                        *(nuvp + i * newg->ix + newg->ix - nwr + j) = *(mgeo + i * gr->ix + 1 + j);
                     }
                 }
              }
              else {
                 for(i=0; i < gr->iy; i++){
                     memcpy(nuvp + i * newg->ix + nwr, mgeo + i * gr->ix, gr->ix*sizeof(float));
                     for(j=0; j < nwr; j++){
                        *(nuvp + i * newg->ix + j) = *(mgeo + i * gr->ix + gr->ix - nwr + j);
                        *(nuvp + i * newg->ix + newg->ix - nwr + j) = *(mgeo + i * gr->ix + j);
                     }
                 }
              }

              gtmp = gr;
              gr = newg;

              surfit(&sm, itpfc, smty, coknu, rsu, sdu);

              gr = gtmp;

           }

           else {

              ap = mgeo;

	      surfit(&sm, itpfc, smty, coknu, ssu, ssdu);

           }


           itpfc = 1;

           for(i=0; i < gr->iy; i++){

               if(smty){
                  xx = *(gr->ygrid + i);
	          xx = FP_PI2 - xx * FP_PI;
                  if(xx < 0.) xx = 0.0;
	          else if(xx > FPI) xx = FPI;
               }

               else yy = *(gr->ygrid + i);

               if(tom == 'g') {
                  coslat = *(coslt + i); 
                  sinlat = *(sinlt + i); 
                  cc = coslat * coslat;
                  sc = sinlat * coslat;
                  ff = fft * sinlat;
               }

               for(j=0; j < gr->ix; j++){

                  if(smty){
		     yy = *(gr->xgrid + j); 
		     yy *= FP_PI;	 
	           }
		   else {
		     xx = *(gr->xgrid + j);
	           }

#ifdef  NOUNDERSCORE

                   bisp(&gz, udd, coknu->tx, &coknu->nx, coknu->ty, &coknu->ny, coknu->c, &coknu->ncof, &xx, &yy, &idif);
#else

                   bisp_(&gz, udd, coknu->tx, &coknu->nx, coknu->ty, &coknu->ny, coknu->c, &coknu->ncof, &xx, &yy, &idif);

#endif
                   if(smty){
                       *(vor + i * gr->ix + j) = (udd[3] + sc * udd[0] + cc * udd[2]) / (cc * ff); 
                   }
                   else {
                       *(vor + i * gr->ix + j) = (udd[2] - FP_PI * sc * udd[1] + cc * udd[3]) / (cc * ff * fpp); 
                   }

               }


           }

           if(iper == 'y') {
              for(i=0; i < gr->iy; i++) *(vor + i * gr->ix + gr->ix - 1) = *(vor + i * gr->ix);
           }

#else

           printf("****WARNING****, no surface fitting or calculation of geostrophic vorticity performed.\r\n");

#endif

        }

        else {

/* computing speeds or kinetic energy */
           if(!ictype){

              if(imiss == 'y'){
                 for(i=0; i < dim; i++){
                     u1 = *(uap + i);
                     v1 = *(vap + i);
                     if(missing(u1, mval, icmp) || missing(v1, mval, icmp)) *(vor + i) = mdval;
                     else *(vor + i) = sqrt(u1 * u1 + v1 * v1);
                 }
              }
              else {
                 for(i=0; i < dim; i++){
                     u1 = *(uap + i);
                     v1 = *(vap + i);
                     *(vor + i) = sqrt(u1 * u1 + v1 * v1);
                 }
              }

           }

           else if(ictype == 1){

              if(imiss == 'y'){
                 for(i=0; i < dim; i++){
                     u1 = *(uap + i);
                     v1 = *(vap + i);
                     if(missing(u1, mval, icmp) || missing(v1, mval, icmp)) *(vor + i) = mdval;
                     else *(vor + i) = 0.5 * (u1 * u1 + v1 * v1);
                 }

              }
              else {              
                 for(i=0; i < dim; i++){
                     u1 = *(uap + i);
                     v1 = *(vap + i);
                     *(vor + i) = 0.5 * (u1 * u1 + v1 * v1);
                 }

              }

           }
           else if(ictype == 4){
              if(imiss == 'y'){
                for(i=0; i < dim; i++) { 
                    spp = *(mgeo + i); 
                    if(missing(spp, mval, icmp)) *(vor + i) = mdval;
                    else *(vor + i) = 0.5 * spp * spp;
                }
              }
              else{
                for(i=0; i < dim; i++) { spp = *(mgeo + i); *(vor + i) = 0.5 * spp * spp;}
              }
           }

        }

        fprintf(fvor, "FRAME %6d\n", fruser);
        fwrite(vor, dim * sizeof(float), 1, fvor);
        fprintf(fvor, "\n");

   }

   fseeko(fvor, (off_t)0, FSTART);
   fprintf(fvor, "%6d %6d %6d\n", gr->ix, gr->iy, fruser);

   close_file(fvor, vorf);

   printf("****INFORMATION****, there are %d frames of speed or vorticity.\n\n", fruser);

   x1u = 1;
   x2u = gr->ix;
   xmn = *(gr->xgrid);
   xmx = *(gr->xgrid + gr->ix - 1);

   if(!ofil) {
      if(ictype < 3){
         if(form != 4) {
            close_file(uu, ufil);
            close_file(vv, vfil);
         }
         else {
            netcdf_close((NETCDF_INFO *)uu);
            netcdf_close((NETCDF_INFO *)vv);
         }
      }
      else{
         if(form != 4) close_file(fgeo, gfil);
         else netcdf_close((NETCDF_INFO *)fgeo);
      }
   }


   if(newg){
      free(newg->xgrid);
      free(newg);
   }

   free(uap); free(vap); free(vor);
   free(nuvp);
   free(mgeo);
   free(coslt);
   free(sinlt);
   free(tanr);

   if(ictype == 2 || ictype == 3){

      if(smty) surfit(&sm, -1, smty, coknu, ssu, ssdu);
      else surfit(&sm, -1, smty, coknu, rsu, sdu);

      if(ictype == 2){

         if(smty) surfit(&sm, -1, smty, coknv, ssv, ssdv);
         else surfit(&sm, -1, smty, coknv, rsv, sdv);
   
      }

      free(coknu); free(coknv);
      free(ssu); free(ssv);
      free(ssdu); free(ssdv);
      free(rsu); free(rsv);
      free(sdu); free(sdv);

   }

   return;

}
