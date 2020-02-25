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

/* function to compute gradient related fields */

#ifdef  REGULAR

#ifdef  NOUNDERSCORE

void bisp(double * , double * , double * , int * , double * , int * , 
          double * , int * , double * , double * , int * );
	  
#else

void bisp_(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#endif

#endif

#define  TOLGRAD  1.0e-10

void surfit(double * , int , int , ... );
int query_interp(int );
int read_field(float * ,float * , float , FILE * , int , int , int , int , int );
void write_header(GRID * , FILE * );

extern GRID *gr;
extern int form;
extern int tl, gof;
extern int frnum;
extern float *ap;
extern float period;
extern int tom;
extern int x1u, x2u, y1u, y2u;
extern float xmn, ymn, xmx, ymx;

double *apd=NULL;

void compute_gradient(FILE *fst, off_t pl, int iper)
{

   int i=0, j=0, ii=0;
   int ofil=0, uvfil=0;
   int smty=0;
   int nwr=0;
   int nff=0;
   int idif=1;
   int itpfc=0;
   int ir=0;
   int fruser=0, rf=0;
   int fs, fr, fe;
/*   int wfs[2], wfr[2], wfe[2]; */
   int dim=gr->ix*gr->iy, dimn=0;
   int itfp=0, igrdo=0;
   int icrf=0, ictaf=0;
   int maxe=1, nfu=0, nfv=0;
   
   int igw=0, igwt=0;
   
   int fuvs[2], fuvr[2], fuve[2];

   off_t fplace1=0, fplace2=0, fchrnum=0;
   off_t uplace1=0, vplace1=0, uplace2=0, vplace2=0;
   off_t uchrnum=0, vchrnum=0;
   

   float *cfld=NULL, *grad=NULL;
   float *uap=NULL, *vap=NULL;
   float *ncfp=NULL;
   float *ggg=NULL;

   double *coslt=NULL, *sinlt=NULL;
   double sm=0.0;
   double coslat=0.0;
   double fdd[4];
   double xx, yy, fz;
   double c1, c2;
   double *gxcmp=NULL, *gycmp=NULL, *gradd=NULL, *gradtfa=NULL;
   double *ncfpd=NULL;
   double sc1=EARTH_RADIUS_M, sc2=EARTH_RADIUS_M*FP_PI;
   
   double cori=0.0;

   FILE *ffld=NULL, *fgrd=NULL, *ftfp=NULL;
   FILE *fgrdx=NULL, *fgrdy=NULL;
   FILE *fcr=NULL, *fta=NULL;
   FILE *uu=NULL, *vv=NULL;

   GRID *newg=NULL, *gtmp=NULL;

   char field[MAXCHR];
   char fgrad[MAXCHR];
   char ftfpf[MAXCHR], fcrf[MAXCHR], ftaf[MAXCHR];
   char fgrdc[MAXCHR], fgrdcx[MAXCHR], fgrdcy[MAXCHR];
   char ufil[MAXCHR], vfil[MAXCHR];

   struct sp_dat *coknf=NULL;        /* spline data for any of the methods   */
   struct rspline *rsf=NULL;         /* data structure for use with smoopy   */
   struct sspline *ssf=NULL;         /* data structure for use with sphery   */
   struct savedat *sdf=NULL;         /* data structure for use with smoopy   */
   struct savedat_sphy *ssdf=NULL;   /* data structure for use with sphery   */

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

   if(tom == 'g'){
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

   printf("Does current file have the required field? 'y' or 'n'\n\n");
   scanf("\n");
   if(getchar() == 'y'){
      ofil = 1;

      if(form != 4) {
         ffld = fst;
         fseeko(fst, pl, FSTART);
         fplace1 = pl;
      }
      else{
         printf("Choose field for which gradient is required.\n\n");
         ffld = (FILE *)nc_clone((NETCDF_INFO *)fst, NULL, NC_SAME);
         ((NETCDF_INFO *)ffld)->iframe = 0;
      }
   }
   else {
      printf("Which file contains required field?\n\n");
      scanf("%s", field);

      if(form != 4) {
         ffld = open_file(field, "r");
         fseeko(ffld, pl, FSTART);
         fplace1 = pl;	  
      }
      else{
         printf("Choose field for which gradient is required.\n\n");
         ffld = (FILE *)nc_clone((NETCDF_INFO *)fst, field, NC_OPEN_MODE);
         ((NETCDF_INFO *)ffld)->iframe = 0;
      }
   }

   printf("If current field is temperature, do you want the thermal front parameter (TFP), 'y' or 'n'?\n\n");
   scanf("\n");
   if(getchar() == 'y') {
      itfp = 1;

      printf("Compute Clarke and Renard frontal parameters, '0' for no or '1' for yes.\n\n");
      scanf("%d", &icrf);
      if(icrf < 0 || icrf > 1){
         printf("****ERROR****, incorrect value for this calculation.\n\n");
         exit(1);
      }

      printf("Compute temperature advection front locator (McCann and Whistler), '0' for no or '1' for yes.\n\n");
      scanf("%d", &ictaf);
      if(ictaf < 0 || ictaf > 1){
         printf("****ERROR****, incorrect value for this calculation.\n\n");
         exit(1);
      }

   }

   printf("Do you want to output the components of the gradient, 'y' or 'n'\n\n");
   scanf("\n");
   if(getchar() == 'y') {
      igrdo=1;
      printf("What is the filename stub required for the output of the gradient components?\n\n");
      scanf("%s", fgrdc);
      strcpy(fgrdcx, fgrdc);
      strcpy(fgrdcy, fgrdc);
      strcat(fgrdcx, "_X.dat");
      strcat(fgrdcy, "_Y.dat");
      
      if(tom == 'g'){
         printf("If field is geopotential, height or pressure convert to geostrophic winds, 'y' or 'n'\n\n");
         scanf("\n");
         if(getchar() == 'y'){
            igw = 1;
	    printf("Is field geopotential, height or pressure, '0' for geopotential, '1' for height, '2' for pressure.\n\n");
	    scanf("%d", &igwt);
	    if(igwt < 0 || igwt > 2){
	       printf("****ERROR****, incorrect Id. for type of field.\n\n");
	       exit(1);
	    }
         }
      }
      
   }

   printf("What are the start, rate and end frame Id's for the chosen field?\n\n");
   scanf("%d %d %d", &fs, &fr, &fe);

   if(ictaf){
   
      printf("Does current file have both U and V components? 'y' or 'n'\n\n");
      scanf("\n");
      if(getchar() == 'y'){
        uvfil = 1;
	
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
   
         printf("****INFORMATION****, Open files for reading wind information.\n\n");
      
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
      scanf("%d %d %d", fuvs, fuvr, fuve);
      printf("What are the start, rate and end frame Id's for the V field?\n\n");
      scanf("%d %d %d", fuvs+1, fuvr+1, fuve+1);
      maxe = (fuve[0] > fuve[1]) ? fuve[0]: fuve[1];
      if(fuvr[0] != fuvr[1]){
         printf("****WARNING****, frame rates are not identical.\b\n\n");
      }

   }

   printf("What is the output file name for the gradient?\n\n");
   scanf("%s", fgrad);

   if(itfp){
      printf("What is the output file name for the TFP data?\n\n");
      scanf("%s", ftfpf);

      if(icrf){
         printf("What is the output file name for the Clarke and Renard data?\n\n");
         scanf("%s", fcrf);
      }
      
      if(ictaf){
         printf("What is the output file name for the McCann and Whistler data?\n\n");
         scanf("%s", ftaf);
      }      
      
   }

/* assign memory for input and gradient field */

   cfld = (float *)calloc(dim, sizeof(float));
   mem_er((cfld == NULL) ? 0: 1, dim * sizeof(float));

   grad = (float *)calloc(dim, sizeof(float));
   mem_er((grad == NULL) ? 0: 1, dim * sizeof(float));
   
   if(ictaf){
   
      uap = (float *)calloc(dim, sizeof(float));
      mem_er((uap == NULL) ? 0: 1, dim * sizeof(float));  

      vap = (float *)calloc(dim, sizeof(float));
      mem_er((vap == NULL) ? 0: 1, dim * sizeof(float));
   
   }

   if(itfp || igrdo || ictaf) {
      gxcmp = (double *)calloc(dim, sizeof(double));
      mem_er((gxcmp == NULL) ? 0: 1, dim * sizeof(double));

      gycmp = (double *)calloc(dim, sizeof(double));
      mem_er((gycmp == NULL) ? 0: 1, dim * sizeof(double));      

      gradd = (double *)calloc(dim, sizeof(double));
      mem_er((gradd == NULL) ? 0: 1, dim * sizeof(double));
      
      if(ictaf){
         gradtfa = (double *)calloc(dim, sizeof(double));
         mem_er((gradtfa == NULL) ? 0: 1, dim * sizeof(double));
      }

   }

/* setup interpolation */

   smty = query_interp(0);

   coknf = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
   mem_er((coknf == NULL) ? 0 : 1, sizeof(struct sp_dat));

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

      ssf = (struct sspline * )malloc_initl(sizeof(struct sspline));
      mem_er((ssf == NULL) ? 0 : 1, sizeof(struct sspline));

      ssdf = (struct savedat_sphy * )malloc_initl(sizeof(struct savedat_sphy));
      mem_er((ssdf == NULL) ? 0 : 1, sizeof(struct savedat_sphy));
   }

   else {

       rsf = (struct rspline * )malloc_initl(sizeof(struct rspline));
       mem_er((rsf == NULL) ? 0 : 1, sizeof(struct rspline));

       sdf = (struct savedat * )malloc_initl(sizeof(struct savedat));
       mem_er((sdf == NULL) ? 0 : 1, sizeof(struct savedat));

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

       ncfp=(float *)calloc(dimn, sizeof(float));
       mem_er((ncfp == NULL) ? 0: 1, dimn * sizeof(float));  

       if(itfp || ictaf){
          ncfpd=(double *)calloc(dimn, sizeof(double));
          mem_er((ncfpd == NULL) ? 0: 1, dimn * sizeof(double));
       }
       

       x1u = 1;
       x2u = newg->ix;
       xmn = *(newg->xgrid);
       xmx = *(newg->xgrid + newg->ix - 1);

   }

/* open file for writing the gradient field */

   fgrd = open_file(fgrad, "w");

   write_header(gr, fgrd);

   if(itfp){
      ftfp = open_file(ftfpf, "w");
      write_header(gr, ftfp);

      if(icrf){
         fcr = open_file(fcrf, "w");
         write_header(gr, fcr);
      }
      
      if(ictaf){
         fta = open_file(ftaf, "w");
         write_header(gr, fta);
      }
   }

   if(igrdo){

      fgrdx = open_file(fgrdcx, "w");
      write_header(gr, fgrdx);
      fgrdy = open_file(fgrdcy, "w");
      write_header(gr, fgrdy);

   }

   while(nff <= fe && nfu <= maxe && nfv <= maxe && !ir){

        if(!nff){

           if((rf = read_field(cfld, NULL, 1.0, ffld, 1, 'n', 'n', '0', 'y'))) ir= 1;
       
           if(form != 4) {fplace2 = ftello(ffld); fchrnum = fplace2 - fplace1;}
       
           if(fs > 1){
             if(form != 4) fseeko(ffld,(fs-2)*fchrnum, ORIGIN);
	     else ((NETCDF_INFO *)ffld)->iframe = fs - 1;	  
	     if((rf = read_field(cfld, NULL, 1.0, ffld, 1, 'n', 'n', '0', 'y'))) ir = 1;
             if(form != 4) fplace2 = ftello(ffld);
           }
	    
           nff = fs;

           if(ofil && form != 4) fseeko(ffld, fplace1, FSTART);

           fruser = 1;

        }

        else {

           if(ofil){
             if(form != 4) fseeko(ffld, fplace2 + (fr-1)*fchrnum, FSTART);
             else ((NETCDF_INFO *)ffld)->iframe = nff - 1;
             if((rf = read_field(cfld, NULL, 1.0, ffld, 1, 'n', 'n', '0', 'y'))) ir = 1;
             if(form != 4) fplace2 = ftello(ffld);
           }

           else {
             if(form != 4) fseeko(ffld, (fr-1)*fchrnum, ORIGIN);
             else ((NETCDF_INFO *)ffld)->iframe = nff - 1;
             if((rf = read_field(cfld, NULL, 1.0, ffld, 1, 'n', 'n', '0', 'y'))) ir = 1; 
           }
          
           if(!ir) ++fruser;

        }
	
	
/* read winds if needed */

        if(ictaf){
	
	   if(!nfu && !nfv){

/* read U */

             if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir= 1;
       
             if(form != 4) {uplace2 = ftello(uu); uchrnum = uplace2 - uplace1;}
       
             if(fuvs[0] > 1){
               if(form != 4) fseeko(uu,(fuvs[0]-2)*uchrnum, ORIGIN);
	       else ((NETCDF_INFO *)uu)->iframe = fuvs[0] - 1;	  
	       if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir = 1;
               if(form != 4) uplace2 = ftello(uu);
             }
	    
             nfu = fuvs[0];

             if(uvfil && form != 4) fseeko(uu, uplace1, FSTART);

/* read V */

             if((rf = read_field(vap, NULL, 1.0, vv, 1, 'n', 'n', '0', 'y'))) ir = 1;
          
             if(form != 4) {vplace2 = ftello(vv); vchrnum = vplace2 - vplace1;}

             if(uchrnum != vchrnum){
                printf("****ERROR****, field block sizes do not match for U and V. File %s, line %d\n\n", __FILE__, __LINE__);
                exit(1);
             }
       
             if(fuvs[1] > 1){
               if(form != 4) fseeko(vv,(fuvs[1]-2)*vchrnum, ORIGIN);
	       else ((NETCDF_INFO *)vv)->iframe = fuvs[1] - 1;	  
	       if((rf = read_field(vap, NULL, 1.0, vv, 1, 'n', 'n', '0', 'y'))) ir = 1;
               if(form != 4) vplace2 = ftello(uu);
             }
	    
             nfv = fuvs[1];          

	   }
           else {

             if(uvfil){
                if(form != 4) fseeko(uu, uplace2 + (fuvr[0]-1)*uchrnum, FSTART);
                else ((NETCDF_INFO *)uu)->iframe = nfu - 1;
                if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir = 1;
                if(form != 4) uplace2 = ftello(uu);

                if(form != 4) fseeko(vv, vplace2 + (fuvr[1]-1)*vchrnum, FSTART);
                else ((NETCDF_INFO *)vv)->iframe = nfv - 1;
                if((rf = read_field(vap, NULL, 1.0, vv, 1, 'n', 'n', '0', 'y'))) ir = 1;
                if(form != 4) vplace2 = ftello(vv);

             }

             else {

                if(form != 4) fseeko(uu, (fuvr[0]-1)*uchrnum, ORIGIN);
                else ((NETCDF_INFO *)uu)->iframe = nfu - 1;
                if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir = 1; 

                if(form != 4) fseeko(vv, (fuvr[1]-1)*vchrnum, ORIGIN);
                else ((NETCDF_INFO *)vv)->iframe = nfv - 1;
                if((rf = read_field(vap, NULL, 1.0, vv, 1, 'n', 'n', '0', 'y'))) ir = 1;
             }
	     
	   }
	    
	}

        if(ir) continue;

        nff += fr;

        if(!smty) {

           ap = ncfp;

           if(iper == 'y'){
              for(i=0; i < gr->iy; i++){
                  memcpy(ncfp + i * newg->ix + nwr, cfld + i * gr->ix, gr->ix*sizeof(float));
                  for(j=0; j < nwr; j++){
                     *(ncfp + i * newg->ix + j) = *(cfld + i * gr->ix + gr->ix - 1 - nwr + j);
                     *(ncfp + i * newg->ix + newg->ix - nwr + j) = *(cfld + i * gr->ix + 1 + j);
                  }
              }
           }
           else {
              for(i=0; i < gr->iy; i++){
                  memcpy(ncfp + i * newg->ix + nwr, cfld + i * gr->ix, gr->ix*sizeof(float));
                  for(j=0; j < nwr; j++){
                     *(ncfp + i * newg->ix + j) = *(cfld + i * gr->ix + gr->ix - nwr + j);
                     *(ncfp + i * newg->ix + newg->ix - nwr + j) = *(cfld + i * gr->ix + j);
                  }
              }
           }

           gtmp = gr;
           gr = newg;

           surfit(&sm, itpfc, smty, coknf, rsf, sdf);

           gr = gtmp;

        }

        else {
           ap = cfld;

           surfit(&sm, itpfc, smty, coknf, ssf, ssdf);

        }

        ap = NULL;

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

#ifdef REGULAR

#ifdef  NOUNDERSCORE

               bisp(&fz, fdd, coknf->tx, &coknf->nx, coknf->ty, &coknf->ny, coknf->c, &coknf->ncof, &xx, &yy, &idif);

#else

               bisp_(&fz, fdd, coknf->tx, &coknf->nx, coknf->ty, &coknf->ny, coknf->c, &coknf->ncof, &xx, &yy, &idif);

#endif

#else
                printf("***ERROR***, surface fitting impossible unless correct libraries\r\n"
                       "             are linked, see compilation options.               \n\n");
                exit(1);

#endif

               if(smty){
                 if(tom == 'g') {
                   c1 = fdd[1] / coslat;
                   c2 = fdd[0];
                 }
                 else {
                   c1 = fdd[1];
                   c2 = fdd[0];
 /* this one does not really make sense for chosen options */
                 }
                  
               }
               else {
                 if(tom == 'g') {
                   c1 = fdd[0] / coslat;
                   c2 = -fdd[1];
                 }
                 else {
                   c1 = fdd[0];
                   c2 = fdd[1];
                 }
               }

               if(itfp || igrdo || ictaf){
                 *(gxcmp + i * gr->ix + j) = c1;
                 *(gycmp + i * gr->ix + j) = c2;
                 *(gradd + i * gr->ix + j) = sqrt(c1*c1 + c2 * c2);
                 *(grad + i * gr->ix + j) = *(gradd + i * gr->ix + j);
               }
               else
                 *(grad + i * gr->ix + j) = sqrt(c1*c1 + c2 * c2);


               if(smty) { 
                  if(tom == 'g') {
                     *(grad + i * gr->ix + j) /= sc1;
                     if(itfp || igrdo){
                        *(gradd + i * gr->ix + j) /= sc1;
                        *(gxcmp + i * gr->ix + j) /= sc1;
                        *(gycmp + i * gr->ix + j) /= sc1;
                     }
                  }
               }
               else {
                  if(tom == 'g') {
                     *(grad + i * gr->ix + j) /= sc2;
                     if(itfp || igrdo){
                        *(gradd + i * gr->ix + j) /= sc2;
                        *(gxcmp + i * gr->ix + j) /= sc2;
                        *(gycmp + i * gr->ix + j) /= sc2;
                     }
                  }
               }
            }

        }

        if(iper == 'y') {
           for(i=0; i < gr->iy; i++) *(grad + i * gr->ix + gr->ix -1) = *(grad + i * gr->ix);
        }

        fprintf(fgrd, "FRAME %6d\n", fruser);
        fwrite(grad, dim * sizeof(float), 1, fgrd);
        fprintf(fgrd, "\n");

        if(igrdo){

	   if(igw && tom == 'g'){
	      for(i=0; i < gr->iy; i++){
	          cori = OMEGA * *(sinlt + i);
		  ii = i * gr->ix;
	          for(j=0; j < gr->ix; j++){
		      ggg = grad + ii + j;
		      *ggg = *(gycmp + ii + j);
		      if(!igwt){
                         *ggg = *ggg / cori ;
		      }
		      else if(igwt == 1){
		          *ggg =  GGRAV * *ggg / cori ;      
		      }
		      else if (igwt == 2){
		          *ggg = *ggg / (cori * RHO);
		      }
		  }
	      }
	   }
	   else{
	      for(i=0; i < dim; i++) *(grad + i) = *(gxcmp + i);
	   }
           fprintf(fgrdx, "FRAME %6d\n", fruser);
           fwrite(grad, dim * sizeof(float), 1, fgrdx);
           fprintf(fgrdx, "\n");
	   if(igw && tom == 'g'){
	      for(i=0; i < gr->iy; i++){
	          coslat = *(coslt + i);
	          cori = OMEGA * *(sinlt + i);
		  ii = i * gr->ix;
	          for(j=0; j < gr->ix; j++){
		      ggg = grad + ii + j;
		      *ggg = *(gxcmp + ii + j);
		      if(!igwt){
                         *ggg = *ggg / cori;
		      }
		      else if(igwt == 1){
		          *ggg = GGRAV * (*ggg) / cori;
		      }
		      else if (igwt == 2){
		          *ggg = *ggg / (cori *  RHO);      
		      }
		  }
	      }
	   }
	   else{
	      for(i=0; i < dim; i++) *(grad + i) = *(gycmp + i);
	   }	   
           fprintf(fgrdy, "FRAME %6d\n", fruser);
           fwrite(grad, dim * sizeof(float), 1, fgrdy);
           fprintf(fgrdy, "\n");
        }

/* thermal front parameter calculation */

        if(itfp){

           if(!smty) {

              apd = ncfpd;

              if(iper == 'y'){
                 for(i=0; i < gr->iy; i++){
                     memcpy(ncfpd + i * newg->ix + nwr, gradd + i * gr->ix, gr->ix*sizeof(double));
                     for(j=0; j < nwr; j++){
                        *(ncfpd + i * newg->ix + j) = *(gradd + i * gr->ix + gr->ix - 1 - nwr + j);
                        *(ncfpd + i * newg->ix + newg->ix - nwr + j) = *(gradd + i * gr->ix + 1 + j);
                     }
                 }
              }
              else {
                 for(i=0; i < gr->iy; i++){
                     memcpy(ncfpd + i * newg->ix + nwr, gradd + i * gr->ix, gr->ix*sizeof(double));
                     for(j=0; j < nwr; j++){
                        *(ncfpd + i * newg->ix + j) = *(gradd + i * gr->ix + gr->ix - nwr + j);
                        *(ncfpd + i * newg->ix + newg->ix - nwr + j) = *(gradd + i * gr->ix + j);
                     }                    
                 }
              }

              gtmp = gr;
              gr = newg;

              surfit(&sm, itpfc, smty, coknf, rsf, sdf);

              gr = gtmp;

           }

           else {
              apd = gradd;

              surfit(&sm, itpfc, smty, coknf, ssf, ssdf);

           }

           apd = NULL;

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

#ifdef REGULAR

#ifdef  NOUNDERSCORE

                  bisp(&fz, fdd, coknf->tx, &coknf->nx, coknf->ty, &coknf->ny, coknf->c, &coknf->ncof, &xx, &yy, &idif);

#else

                  bisp_(&fz, fdd, coknf->tx, &coknf->nx, coknf->ty, &coknf->ny, coknf->c, &coknf->ncof, &xx, &yy, &idif);

#endif

#else
                printf("***ERROR***, surface fitting impossible unless correct libraries\r\n"
                       "             are linked, see compilation options.               \n\n");
                exit(1);
#endif

                  if(smty){
                    if(tom == 'g') {
                      c1 = fdd[1] / coslat;
                      c2 = fdd[0];
                    }
                    else {
                      c1 = fdd[1];
                      c2 = fdd[0];
                    }
                  
                  }
                  else {
                    if(tom == 'g') {
                      c1 = fdd[0] / coslat;
                      c2 = -fdd[1];
                    }
                    else {
                      c1 = fdd[0];
                      c2 = fdd[1];
                    }
                  }

                  if(fabs(*(gradd + i * gr->ix + j)) < TOLGRAD) *(grad + i * gr->ix + j) = 0.0;
                  else {
                     *(gradd + i * gr->ix + j) = - (*(gxcmp + i * gr->ix + j) * c1 + *(gycmp + i * gr->ix + j) * c2) / *(gradd + i * gr->ix + j);
/*                     *(gradd + i * gr->ix + j) = - (*(gxcmp + i * gr->ix + j) * c1 + *(gycmp + i * gr->ix + j) * c2); */
                     if(smty) { if(tom == 'g') *(gradd + i * gr->ix + j) /= sc1;}
                     else {if(tom == 'g') *(gradd + i * gr->ix + j) /= sc2;}
                     *(grad + i * gr->ix + j) = *(gradd + i * gr->ix + j);
                  }

               }

           }

           if(iper == 'y') {
              for(i=0; i < gr->iy; i++) {
                  *(grad + i * gr->ix + gr->ix -1) = *(grad + i * gr->ix);
                  *(gradd + i * gr->ix + gr->ix -1) = *(gradd + i * gr->ix);
              }
           }

           fprintf(ftfp, "FRAME %6d\n", fruser);
           fwrite(grad, dim * sizeof(float), 1, ftfp);
           fprintf(ftfp, "\n");

/* compute Clarke and Renard */

           if(icrf){

              if(!smty) {

                 apd = ncfpd;

                 if(iper == 'y'){
                    for(i=0; i < gr->iy; i++){
                        memcpy(ncfpd + i * newg->ix + nwr, gradd + i * gr->ix, gr->ix*sizeof(double));
                        for(j=0; j < nwr; j++){
                           *(ncfpd + i * newg->ix + j) = *(gradd + i * gr->ix + gr->ix - 1 - nwr + j);
                           *(ncfpd + i * newg->ix + newg->ix - nwr + j) = *(gradd + i * gr->ix + 1 + j);
                        }
                    }
                 }
                 else {
                    for(i=0; i < gr->iy; i++){
                        memcpy(ncfpd + i * newg->ix + nwr, gradd + i * gr->ix, gr->ix*sizeof(double));
                        for(j=0; j < nwr; j++){
                           *(ncfpd + i * newg->ix + j) = *(gradd + i * gr->ix + gr->ix - nwr + j);
                           *(ncfpd + i * newg->ix + newg->ix - nwr + j) = *(gradd + i * gr->ix + j);
                        }                    
                    }
                 }

                 gtmp = gr;
                 gr = newg;

                 surfit(&sm, itpfc, smty, coknf, rsf, sdf);

                 gr = gtmp;

              }

              else {
                 apd = gradd;

                 surfit(&sm, itpfc, smty, coknf, ssf, ssdf);

              }
	      
              apd = NULL;

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

#ifdef REGULAR

#ifdef  NOUNDERSCORE

                     bisp(&fz, fdd, coknf->tx, &coknf->nx, coknf->ty, &coknf->ny, coknf->c, &coknf->ncof, &xx, &yy, &idif);

#else

                     bisp_(&fz, fdd, coknf->tx, &coknf->nx, coknf->ty, &coknf->ny, coknf->c, &coknf->ncof, &xx, &yy, &idif);

#endif

#else
                     printf("***ERROR***, surface fitting impossible unless correct libraries\r\n"
                       "             are linked, see compilation options.               \n\n");
                     exit(1);
#endif

                     if(smty){
                       if(tom == 'g') {
                         c1 = fdd[1] / coslat;
                         c2 = fdd[0];
                       }
                       else {
                         c1 = fdd[1];
                         c2 = fdd[0];
                       }
                  
                     }
                     else {
                       if(tom == 'g') {
                         c1 = fdd[0] / coslat;
                         c2 = -fdd[1];
                       }
                       else {
                         c1 = fdd[0];
                         c2 = fdd[1];
                       }

                     }

                     *(gradd + i * gr->ix + j) = sqrt(c1*c1 + c2 * c2);
                     if(fabs(*(gradd + i * gr->ix + j)) < TOLGRAD) *(grad + i * gr->ix + j) = 0.0;
                     else {
                        *(grad + i * gr->ix + j) = (*(gxcmp + i * gr->ix + j) * c1 + *(gycmp + i * gr->ix + j) * c2) / *(gradd + i * gr->ix + j);
                     }

                  }

              }

              if(iper == 'y') {
                 for(i=0; i < gr->iy; i++) {
                  *(grad + i * gr->ix + gr->ix -1) = *(grad + i * gr->ix);
                 }
              }

              fprintf(fcr, "FRAME %6d\n", fruser);
              fwrite(grad, dim * sizeof(float), 1, fcr);
              fprintf(fcr, "\n");

           }
	   
/* compute McCann and Whistler thermal advection front locator */

           if(ictaf){
	      
	      for(i=0; i < gr->iy; i++){

                  if(tom == 'g'){
                     for(j=0; j < gr->ix; j++){
		  
		       *(gradtfa + i * gr->ix + j) = *(gxcmp + i * gr->ix + j) * *(uap + i * gr->ix + j) - *(gycmp + i * gr->ix + j) * *(vap + i * gr->ix + j);
		  
		     }		  	  
		  }
		  else {
                     for(j=0; j < gr->ix; j++){
		  
		       *(gradtfa + i * gr->ix + j) = *(gxcmp + i * gr->ix + j) * *(uap + i * gr->ix + j) + *(gycmp + i * gr->ix + j) * *(vap + i * gr->ix + j);
		  
		     }
		  
		  }	               
	      
	      }
	            
	      
	      if(!smty) {

                 apd = ncfpd;

                 if(iper == 'y'){
                    for(i=0; i < gr->iy; i++){
                        memcpy(ncfpd + i * newg->ix + nwr, gradtfa + i * gr->ix, gr->ix*sizeof(double));
                        for(j=0; j < nwr; j++){
                           *(ncfpd + i * newg->ix + j) = *(gradtfa + i * gr->ix + gr->ix - 1 - nwr + j);
                           *(ncfpd + i * newg->ix + newg->ix - nwr + j) = *(gradd + i * gr->ix + 1 + j);
                        }
                    }
                 }
                 else {
                    for(i=0; i < gr->iy; i++){
                        memcpy(ncfpd + i * newg->ix + nwr, gradtfa + i * gr->ix, gr->ix*sizeof(double));
                        for(j=0; j < nwr; j++){
                           *(ncfpd + i * newg->ix + j) = *(gradtfa + i * gr->ix + gr->ix - nwr + j);
                           *(ncfpd + i * newg->ix + newg->ix - nwr + j) = *(gradd + i * gr->ix + j);
                        }                    
                    }
                 }

                 gtmp = gr;
                 gr = newg;

                 surfit(&sm, itpfc, smty, coknf, rsf, sdf);

                 gr = gtmp;

              }

              else {
                 apd = gradtfa;

                 surfit(&sm, itpfc, smty, coknf, ssf, ssdf);

              }
	      
	      apd = NULL;
	      
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

#ifdef REGULAR

#ifdef  NOUNDERSCORE

                     bisp(&fz, fdd, coknf->tx, &coknf->nx, coknf->ty, &coknf->ny, coknf->c, &coknf->ncof, &xx, &yy, &idif);

#else

                     bisp_(&fz, fdd, coknf->tx, &coknf->nx, coknf->ty, &coknf->ny, coknf->c, &coknf->ncof, &xx, &yy, &idif);

#endif

#else
                     printf("***ERROR***, surface fitting impossible unless correct libraries\r\n"
                       "             are linked, see compilation options.               \n\n");
                     exit(1);
#endif


                     if(smty){
                       if(tom == 'g') {
                         c1 = fdd[1] / coslat;
                         c2 = fdd[0];
                       }
                       else {
                         c1 = fdd[1];
                         c2 = fdd[0];
                       }
                  
                     }
                     else {
                       if(tom == 'g') {
                         c1 = fdd[0] / coslat;
                         c2 = -fdd[1];
                       }
                       else {
                         c1 = fdd[0];
                         c2 = fdd[1];
                       }

                     }

                     if(tom == 'g') *(gradtfa + i * gr->ix + j) = c1 * *(uap + i * gr->ix + j) - c2 * *(vap + i * gr->ix + j);
		     else *(gradtfa + i * gr->ix + j) = c1 * *(uap + i * gr->ix + j) + c2 * *(vap + i * gr->ix + j);
		     if(smty) { if(tom == 'g') *(gradtfa + i * gr->ix + j) /= sc1;}
                     else {if(tom == 'g') *(gradtfa + i * gr->ix + j) /= sc2;}
		     
		     *(grad + i * gr->ix + j) = *(gradtfa + i * gr->ix + j);

                  }

              }
	      
	      if(iper == 'y') {
                 for(i=0; i < gr->iy; i++) {
                  *(grad + i * gr->ix + gr->ix -1) = *(grad + i * gr->ix);
                 }
              }

              fprintf(fta, "FRAME %6d\n", fruser);
              fwrite(grad, dim * sizeof(float), 1, fta);
              fprintf(fta, "\n");
	   
	   }
	  
        }

   }


   fseeko(fgrd, (off_t)0, FSTART);
   fprintf(fgrd, "%6d %6d %6d\n", gr->ix, gr->iy, fruser);

   close_file(fgrd, fgrad);

   if(igrdo){
      fseeko(fgrdx, (off_t)0, FSTART);
      fprintf(fgrdx, "%6d %6d %6d\n", gr->ix, gr->iy, fruser);
      close_file(fgrdx, fgrdcx);
      fseeko(fgrdy, (off_t)0, FSTART);
      fprintf(fgrdy, "%6d %6d %6d\n", gr->ix, gr->iy, fruser);
      close_file(fgrdy, fgrdcy);
   }

   if(itfp){
      fseeko(ftfp, (off_t)0, FSTART);
      fprintf(ftfp, "%6d %6d %6d\n", gr->ix, gr->iy, fruser);
      close_file(ftfp, ftfpf);

      if(icrf){
         fseeko(fcr, (off_t)0, FSTART);
         fprintf(fcr, "%6d %6d %6d\n", gr->ix, gr->iy, fruser);
         close_file(fcr, fcrf);
      }

      if(ictaf){
         fseeko(fta, (off_t)0, FSTART);
         fprintf(fta, "%6d %6d %6d\n", gr->ix, gr->iy, fruser);
         close_file(fta, ftaf);
      }      
      
   }

   printf("****INFORMATION****, there are %d frames of gradient.\n\n", fruser);

   x1u = 1;
   x2u = gr->ix;
   xmn = *(gr->xgrid);
   xmx = *(gr->xgrid + gr->ix - 1);


   if(!ofil) {
     if(form != 4) close_file(ffld, field);
     else netcdf_close((NETCDF_INFO *)ffld);
   }

   if(itfp){
      free(gradd);
      free(gxcmp);
      free(gycmp);
      free(ncfpd);
   }
   
   if(ictaf){
      free(uap);
      free(vap);
   }

   if(newg){
      free(newg->xgrid);
      free(newg);
   }

   free(coslt);
   free(sinlt);

   if(smty) surfit(&sm, -1, smty, coknf, ssf, ssdf);
   else surfit(&sm, -1, smty, coknf, rsf, sdf);

   free(coknf); 
   free(ssf);
   free(ssdf);
   free(rsf);
   free(sdf);


   free(ncfp);
   free(cfld);
   free(grad);

   return;
}
