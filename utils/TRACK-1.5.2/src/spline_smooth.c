#include <Stdio.h>

#ifndef  REGULAR

void spline_smooth()
{

   printf("***error***, surface fitting impossible unless correct libraries\r\n"
          "             are linked, see compilation options.               \n\n");
   exit(1);

}

#else

#include <Math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "grid.h"
#include "mem_er.h"
#include "reg_dat.h"
#include "sphery_dat.h"
#include "bisp.h"
#include "netcdf_info.h"
#include "file_handle.h"
#include "files_out.h"
#include "m_values.h"

/* function to smooth fields using B-splines */


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
void write_header(GRID * , FILE * );

extern GRID *gr;
extern int form;
extern int tl, gof;
extern float *ap;
extern int frnum;

extern char *fext;
extern int iext;

void spline_smooth(FILE *fin)
{

   int i, j;
   int istr=0, intv=0, ilst=0;
   int smty=0;
   int frc=0, rf=0;
   int itpfc=0;
   int fruser=0;
   int idif=0;

   off_t place1=0, place2=0, chrnum=0;

   char smooth[MAXCHR];
   char newgf[MAXCHR];

   FILE *fout = NULL;
   FILE *fgn = NULL;
   
   GRID *gnew=NULL;

   float scale=0.0;
   float *newf=NULL;

   double sm=0.0;
   double xx, yy, zz, dd[4];
   double xa, ya;

   struct sp_dat *cokn=NULL;        /* spline data for any of the methods   */
   struct rspline *rs=NULL;         /* data structure for use with smoopy   */
   struct sspline *ss=NULL;         /* data structure for use with sphery   */
   struct savedat *sd=NULL;         /* data structure for use with smoopy   */
   struct savedat_sphy *ssd=NULL;   /* data structure for use with sphery   */

   strncpy(smooth, SMOOTH, MAXCHR);
   if(iext) strcpy(strstr(smooth, EXTENSION), fext);

   ap=(float *)calloc(gr->ix * gr->iy, sizeof(float));
   mem_er((ap == NULL) ? 0: 1, gr->ix * gr->iy * sizeof(float));     

   printf("What are the start, interval and final frame numbers required.\n\n");
   scanf("%d %d %d", &istr, &intv, &ilst);

   printf("What scaling is required for chosen fields?\n\n");
   scanf("%f", &scale);
   
   printf("Is a new grid required for writing smoothed field, 'y' or 'n'\n\n");
   scanf("\n");
   
   if(getchar() == 'y'){
      printf("What file contains the new grid information.\n\n");
      scanf("%s", newgf);
      
      fgn = open_file(newgf, "r");
      gnew = (GRID *)malloc_initl(sizeof(GRID));
      mem_er((gnew == NULL) ? 0 : 1, sizeof(GRID));
      fscanf(fgn, "%d %d", &(gnew->ix), &(gnew->iy));
      gnew->xgrid = (float *)calloc(gnew->ix, sizeof(float));
      mem_er((gnew->xgrid == NULL) ? 0 : 1, gnew->ix * sizeof(float));
      gnew->ygrid = (float *)calloc(gnew->iy, sizeof(float));
      mem_er((gnew->ygrid == NULL) ? 0 : 1, gnew->iy * sizeof(float));
      
      for(i=0; i < gnew->ix; i++) fscanf(fgn, "%f", gnew->xgrid + i);
      
      for(i=0; i < gnew->iy; i++) fscanf(fgn, "%f", gnew->ygrid + i);   
      
      close_file(fgn, newgf);
      
/* check new grid lies in the domain of the original grid */

      if((*(gnew->xgrid) < *(gr->xgrid)                             || 
          *(gnew->xgrid + gnew->ix - 1) > *(gr->xgrid + gr->ix - 1))||
	 (*(gnew->ygrid) < *(gr->ygrid)                             ||
	  *(gnew->xgrid + gnew->ix - 1) > *(gr->xgrid + gr->ix - 1))  ){
	  
	  printf("****ERROR****, new grid extends outside the region\r\n"
	         "               of the original grid, no output to \r\n"
		 "               new grid.                          \n\n");
          return;
	  
      }
      
      newf = (float *)calloc(gnew->ix * gnew->iy, sizeof(float));
      mem_er((newf == NULL) ? 0 : 1, gnew->ix * gnew->iy * sizeof(float));
   
   }
   else {gnew = gr; newf = ap;}
  
   
   smty = query_interp(0);

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

   fout = open_file(smooth, "w");

/* write header */

   write_header(gnew, fout);

   frc = 0;

   while(frc <= frnum && !rf){

        if(!frc){
	   rf = read_field(ap, NULL, scale, fin, 1, 'n', 'n', '0', 'y');
       
           if(form != 4) {place2 = ftello(fin); chrnum = place2 - place1;}
       
           if(istr > 1){
              if(form != 4) fseeko(fin,(istr-2)*chrnum, ORIGIN);
	      else ((NETCDF_INFO *)fin)->iframe = istr - 1;	  
	      rf = read_field(ap, NULL, scale, fin, 1, 'n', 'n', '0', 'y');
           }
	    
	   frc = istr;
           fruser = 1;

         }
	 else {

           if(form != 4) fseeko(fin, (intv-1)*chrnum, ORIGIN); 
           else ((NETCDF_INFO *)fin)->iframe = frc - 1;
	  
	   rf = read_field(ap, NULL, scale, fin, 1, 'n', 'n', '0', 'y');

           ++fruser;
	  
	 }

	 if(smty) surfit(&sm, itpfc, smty, cokn, ss, ssd);
         else surfit(&sm, itpfc, smty, cokn, rs, sd);
	 itpfc = 1;

         for(i=0; i < gnew->iy; i++){

             for(j=0; j < gnew->ix; j++){

                 yy = *(gnew->ygrid + i);
                 xx = *(gnew->xgrid + j);

                 if(smty){

                    xx *= FP_PI;
                    yy = FP_PI2 - yy * FP_PI;
                    if(yy < 0.) yy = 0.0;
                    else if(yy > FPI) yy = FPI;

                    xa = yy;
                    ya = xx;

                 }

                 else {xa = xx; ya= yy;}
		 

#ifdef  NOUNDERSCORE

                 bisp(&zz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xa, &ya, &idif);

#else

                 bisp_(&zz, dd, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xa, &ya, &idif);

#endif
                 *(newf + i * gnew->ix + j) = zz;

             }

         }

         fprintf(fout, "FRAME %6d\n", fruser);

         fwrite(newf, gnew->ix * gnew->iy * sizeof(float), 1, fout);
         fprintf(fout, "\n");


	 frc += intv;

         if(frc > ilst) break;

   }

   fseeko(fout, (off_t)0, FSTART);
   fprintf(fout, "%6d %6d %6d\n", gnew->ix, gnew->iy, fruser);

   close_file(fout, smooth);

   if(smty) surfit(&sm, -1, smty, cokn, ss, ssd);
   else surfit(&sm, -1, smty, cokn, rs, sd);
   free(cokn);
   free(ss);
   free(ssd);
   free(rs);
   free(sd);
   free(ap);
   
   if(gr != gnew){
      free(newf);
      free(gnew->xgrid);
      free(gnew->ygrid);
      free(gnew);
   }
      
   return;

}

#endif

