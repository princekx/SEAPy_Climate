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

#include <Stdio.h>
#include <Math.h>
#include <sys/types.h>
#include <string.h>
#include "mem_er.h"
#include "reg_dat.h"
#include "sphery_dat.h"
#include "bisp.h"
#include "file_handle.h"
#include "file_cat_out.h"
#include "grid.h"
#include "netcdf_info.h"
#include "m_values.h"

/* function to interpolate between two different grids. Useful for transforming from non-gaussian to 
   gaussian grids.                                                                                     */
   
#ifdef  NOUNDERSCORE

void bisp(double * , double * , double * , int * , double * , int * , 
          double * , int * , double * , double * , int * );
	  
#else

void bisp_(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#endif  

GRID *gauss_grid(int );
void surfit(double * , int , int , ... );
int query_interp(int );
void write_header(GRID * , FILE * );
int read_field(float * ,float * , float , FILE * , int , int , int , int , int );

   
extern GRID *gr;
extern int form;
extern float *ap;
extern int x1u, x2u, y1u, y2u;
extern float xmn, ymn, xmx, ymx;
extern float period;
   
void interp_ng(FILE *fdat, off_t pl, int iper)
{

   int i=0, j=0;
   
   int smty=0;
   int intyp=0;
   int dimn=0, dimt=0;
   int igt=0;
   int ntrunc=0;
   int fr1=0, fri=0, frl=0;
   int nf=0, fruser=0;
   int nwr=0;
   int itpfc=0;
   int dim=gr->ix*gr->iy;
   int idif=0;
   
   off_t place1=0, place2=0, chrnum=0;
   
   char gfile[MAXCHR];

   FILE *fing=NULL;
   FILE *fout=NULL;
   
   GRID *newg=NULL, *tmpg=NULL;
   GRID *gtmp=NULL;

   char fintrp[MAXCHR];

   float *ffld=NULL;
   float *ncfp=NULL;
   float *ffint=NULL;
   float *app=NULL;
   
   double sm=0.0;
   double dd[4];
   double fz=0.0;
   double xx, yy;

   struct sp_dat *cokn=NULL;        /* spline data for any of the methods   */
   struct rspline *rs=NULL;         /* data structure for use with smoopy   */
   struct sspline *ss=NULL;         /* data structure for use with sphery   */
   struct savedat *sd=NULL;         /* data structure for use with smoopy   */
   struct savedat_sphy *ssd=NULL;   /* data structure for use with sphery   */

   ap = NULL;
   
   x1u = 1;
   x2u = gr->ix;
   y1u = 1;
   y2u = gr->iy;

   xmn = *(gr->xgrid);
   xmx = *(gr->xgrid + gr->ix - 1);
   ymn = *(gr->ygrid);
   ymx = *(gr->ygrid + gr->iy - 1);

   printf("Do you want to read new grid information from file or compute a Gaussian grid, \r\n"
          "input '0' to read from file or '1' to compute a Gaussian grid.                 \n\n");
   scanf("%d", &igt);
   
   if(igt < 0 || igt > 1){
      printf("****ERROR****, incorrect value for option, %d\n", igt);
      exit(1);
   }
   if(!igt){
      printf("What is the name of the grid information file?\n\n");
      scanf("%s", gfile);
      
      newg = (GRID *)malloc_initl(sizeof(GRID));
      mem_er((newg == NULL) ? 0 : 1, sizeof(GRID));

      fing = open_file(gfile, "r");

      fscanf(fing, "%d %d", &(newg->ix), &(newg->iy));

      newg->xgrid = (float *)calloc(newg->ix, sizeof(float));
      mem_er((newg->xgrid == NULL) ? 0 : 1, newg->ix * sizeof(float));

      newg->ygrid = (float *)calloc(newg->iy, sizeof(float));
      mem_er((newg->ygrid == NULL) ? 0 : 1, newg->iy * sizeof(float));

      for(i=0; i < newg->ix; i++) fscanf(fing, "%f", newg->xgrid + i);
      for(i=0; i < newg->iy; i++) fscanf(fing, "%f", newg->ygrid + i);

      close_file(fing, gfile);     
   }
   else{
      printf("What is the spectral truncation (triangular) for the output gaussian grid?\n");
      scanf("%d", &ntrunc);

      if(ntrunc < 0){
         printf("****ERROR****, truncation must be positive, exiting.\n\n");
         exit(1);
      }

      newg = gauss_grid(ntrunc);

      --(newg->ix);
   }
   
/* Check new grid is within the old grid */

   if(*(newg->xgrid) < *(gr->xgrid) || *(newg->xgrid + newg->ix - 1) > *(gr->xgrid + gr->ix - 1)) {
      printf("****WARNING****, new grid extends outside the range of the old grid in the X-coordinate (longtitude)  \r\n"
             "                 this maybe accomodated by periodically extending the data in X if the data is global.\n\n");
   }
   if(*(newg->ygrid) < *(gr->ygrid) || *(newg->ygrid + newg->iy - 1) > *(gr->ygrid + gr->iy - 1)) {
      printf("****WARNING****, new grid extends outside the range of the old grid in the Y-coordinate (latitude)       \r\n"
             "                 if this is because the poles have been removed for global data consider keeping them in.\n\n"); 
   
   }


   ffld = (float *)calloc(dim, sizeof(float));
   mem_er((ffld == NULL) ? 0: 1, dim * sizeof(float));
      
   dimn = newg->ix * newg->iy;

   printf("What are the start, interval and final frame numbers?\n\n");
   scanf("%d %d %d", &fr1, &fri, &frl);

/*   printf("What form of interpolation is required,               \r\n"
          "Input '0' for bi-cubic spline, no missing values.     \r\n"
          "Input '1' for bi-linear interpolation, missing values.\r\n"
          "Input '2' for bi-cubic interpolation, missing values. \n\n");
   scanf("%d", &intyp);
   if(intyp < 0 || intyp > 2){
      printf("****ERROR****, interpolation type unknown, exiting.\n\n");
      exit(1);
   } 
   
   printf("****INFORMATION****, using interpolation option %d\n", intyp);
   
*/

   intyp = 0;
   
   if(!intyp){

      smty = query_interp(0);

      cokn = (struct sp_dat * )malloc_initl(sizeof(struct sp_dat));
      mem_er((cokn == NULL) ? 0 : 1, sizeof(struct sp_dat));

      if(smty){

         ss = (struct sspline * )malloc_initl(sizeof(struct sspline));
         mem_er((ss == NULL) ? 0 : 1, sizeof(struct sspline));

         ssd = (struct savedat_sphy * )malloc_initl(sizeof(struct savedat_sphy))
;
         mem_er((ssd == NULL) ? 0 : 1, sizeof(struct savedat_sphy));
	 
	 ap = ffld;

      }

      else{

         rs = (struct rspline * )malloc_initl(sizeof(struct rspline));
         mem_er((rs == NULL) ? 0 : 1, sizeof(struct rspline));

         sd = (struct savedat * )malloc_initl(sizeof(struct savedat));
         mem_er((sd == NULL) ? 0 : 1, sizeof(struct savedat));

         printf("How much wrap round is required for fitting 0 <= w <= %d\n\n", gr->ix);
         scanf("%d", &nwr);
         if(iper == 'n'){
            printf("****WARNING****, wrapround may not be valid for this data set.\n\n");
         }
         if(nwr < 0 || nwr > gr->ix){
            printf("****ERROR****, invalid wrapround.\n\n");
            exit(1);
         }

         tmpg = (GRID *)malloc_initl(sizeof(GRID));
         mem_er((tmpg == NULL) ? 0 : 1, sizeof(GRID));
         tmpg->ix = gr->ix + 2 * nwr;
         tmpg->iy = gr->iy;
         dimt = tmpg->ix * tmpg->iy;

         tmpg->xgrid = (float * )calloc(tmpg->ix, sizeof(float));
         mem_er((tmpg->xgrid == NULL) ? 0 :1, tmpg->ix * sizeof(float));

         memcpy(tmpg->xgrid + nwr, gr->xgrid, gr->ix * sizeof(float));

         if(iper == 'y'){

            for(i=0; i < nwr; i++) {
               *(tmpg->xgrid + i) = *(gr->xgrid + gr->ix - 1 - nwr + i) - period;
               *(tmpg->xgrid + tmpg->ix - nwr + i) = *(gr->xgrid + 1 + i) + period;
            }

         }

         else{

            for(i=0; i < nwr; i++){
               *(tmpg->xgrid + i) = *(gr->xgrid + gr->ix - nwr + i) - period;
               *(tmpg->xgrid + tmpg->ix - nwr + i) = *(gr->xgrid + i) + period;
            }

         }
	 
	 if(*(newg->xgrid) < *(tmpg->xgrid) || *(newg->xgrid + newg->ix - 1) > *(tmpg->xgrid + tmpg->ix - 1)) {
               printf("****WARNING****, new grid still extends outside the range of the extended old grid in the X-coordinate (longtitude).\r\n"
	              "               consider extending old grid further.                                                                 \n\n");
	 }
	 
         tmpg->ygrid = gr->ygrid;

         ncfp=(float *)calloc(dimt, sizeof(float));
         mem_er((ncfp == NULL) ? 0: 1, dimt * sizeof(float));
	 
	 ap = ncfp; 

         x1u = 1;
         x2u = tmpg->ix;
         xmn = *(tmpg->xgrid);
         xmx = *(tmpg->xgrid + tmpg->ix - 1);

      }

   }
   
   printf("What is the output file name for the interpolated data?\n\n");
   scanf("%s", fintrp);
   
   ffint=(float *)calloc(dimn, sizeof(float));
   mem_er((ffint == NULL) ? 0: 1, dimn * sizeof(float));

   fout = open_file(fintrp, "w");

   write_header(newg, fout);

   nf = 0;

   while(nf <= frl){

      if(form != 4){

         if(!nf) {

            place1 = ftello(fdat);

            if(read_field(ffld, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;

            if(fr1 == 1) ++nf;;

            place2 = ftello(fdat);
            chrnum = place2 - place1;

            if(fr1 > 1){

               fseeko(fdat, (fr1-2)*chrnum, ORIGIN);
               if(read_field(ffld, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;

               nf = fr1;

            }

            nf += fri;
	    
	    fruser = 1;

         }

         else {

            fseeko(fdat, (fri-1)*chrnum, ORIGIN);
            if(read_field(ffld, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;

            nf += fri;
            ++fruser;
         }

      }

      else {

        if(!nf) nf = fr1;
        ((NETCDF_INFO *)fdat)->iframe = nf - 1;
        if(read_field(ffld, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;
        nf += fri;
	++fruser;

      }
      


      if(!intyp){
      
#ifdef  REGULAR

         if(smty) surfit(&sm, itpfc, smty, cokn, ss, ssd);
         else {
	
            if(iper == 'y'){
               for(i=0; i < gr->iy; i++){
                   memcpy(ncfp + i * tmpg->ix + nwr, ffld + i * gr->ix, gr->ix*sizeof(float));
                   for(j=0; j < nwr; j++){
                      *(ncfp + i * tmpg->ix + j) = *(ffld + i * gr->ix + gr->ix - 1 - nwr + j);
                      *(ncfp + i * tmpg->ix + tmpg->ix - nwr + j) = *(ffld + i * gr->ix + 1 + j);
                   }
               }
            }
            else {
               for(i=0; i < gr->iy; i++){
                   memcpy(ncfp + i * tmpg->ix + nwr, ffld + i * gr->ix, gr->ix*sizeof(float));

                   for(j=0; j < nwr; j++){
                      *(ncfp + i * tmpg->ix + j) = *(ffld + i * gr->ix + gr->ix - nwr + j);
                      *(ncfp + i * tmpg->ix + tmpg->ix - nwr + j) = *(ffld + i * gr->ix + j);
                   }
               }
            }

            gtmp = gr;
            gr = tmpg;
	 
	    surfit(&sm, itpfc, smty, cokn, rs, sd);
	    
	    gr = gtmp;
	    
	 }
	 
         itpfc = 1;
	 
         for(i=0; i < newg->iy; i++){

            if(smty){
               xx = *(newg->ygrid + i);
	       xx = FP_PI2 - xx * FP_PI;
               if(xx < 0.) xx = 0.0;
	       else if(xx > FPI) xx = FPI;
            }

            else yy = *(newg->ygrid + i);
	    
	    
            for(j=0; j < newg->ix; j++){

               if(smty){
                  yy = *(newg->xgrid + j); 
	          yy *= FP_PI;	 
	       }
               else {
		  xx = *(newg->xgrid + j);
	       }
	       
	       
	       app = ffint + i * newg->ix + j;

#ifdef  NOUNDERSCORE

               bisp(&fz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#else

               bisp_(&fz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#endif


               *app = (float) fz;
	      

           }
	   
	   if(iper == 'y') {
              for(i=0; i < newg->iy; i++) *(ffint + i * newg->ix + newg->ix -1) = *(ffint + i * newg->ix);
           }

         }
	 
#else
         printf("****WARNING****, no surface fitting performed.\n\n");

#endif

      }
      
      fprintf(fout, "FRAME %d\n", fruser);
      fwrite(ffint, dimn * sizeof(float), 1, fout);
      fprintf(fout, "\n");

   }
   
   fseeko(fout, (off_t)0, FSTART);
   fprintf(fout, "%6d %6d %6d\n", newg->ix, newg->iy, fruser);

   close_file(fout, fintrp);

   printf("****INFORMATION****,  %d frames have been interpolated.\n\n", nf);



   if(!intyp){
      if(smty) surfit(&sm, -1, smty, cokn, ss, ssd);
      else surfit(&sm, -1, smty, cokn, rs, sd);

      free(cokn);
      free(ss);
      free(ssd);
      free(rs);
      free(sd);
   }

   free(tmpg->xgrid);
   free(newg->xgrid);
   free(newg->ygrid);
   free(newg);
   free(tmpg);
   free(ncfp);
   free(ffint); 
 

   return;
}

#endif
