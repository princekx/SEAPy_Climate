#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include <sys/types.h>
#include "grid.h"
#include "file_handle.h"
#include "file_cat_out.h"
#include "mem_er.h"
#include "netcdf_info.h"
#include "m_values.h"
#include "geo_values.h"

/* function to compute vorticty, speed or EKE */


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

void compute_vorticity_fd(FILE *fst, off_t pl, int iper)
{

   int i, j, ii;
   int i1=0, i2=0, j1=0, j2=0;
   int ir=0;
   int nfu=0, nfv=0;
   int fruser=0, rf=0;
   int ofil=0;
   int itpfc=0;
   int maxe=0;
   int dim=gr->ix*gr->iy;
   int nwr=0;
   int imiss='n';
   int icmp=0;
   int addvor='n';
   int imm=0;

   int fs[2], fr[2], fe[2];

   off_t uplace1=0, vplace1=0, uplace2=0, vplace2=0;
   off_t uchrnum=0, vchrnum=0;

   float *uap=NULL, *vap=NULL, *vor=NULL;
   float *lat2d=NULL;
   float mval=0.0, mdval=0.0;
   float vorcs=1.0;
   
   float dx=0.0;

   double *coslt=NULL;
   double *tanr=NULL;
   double coslat=0.0;
   double dy=0.0, dy1=0.0, dy2=0.0;
   double dx1=0.0, dx2=0.0;

   FILE *uu=NULL, *vv=NULL;
   FILE *fvor=NULL;
   FILE *latf=NULL;

   char ufil[MAXCHR], vfil[MAXCHR], latfil[MAXCHR];
   char vorf[MAXCHR];

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

   printf("****INFORMATION****, input fields should be U and V.           \n\n");

   coslt = (double *)calloc(gr->iy, sizeof(double));
   mem_er((coslt == NULL) ? 0 : 1, (gr->iy)*sizeof(double));
   for(i=0; i < gr->iy; i++) *(coslt + i) = cos((double)(*(gr->ygrid + i)) * FP_PI);
   
   dy1 = *(gr->ygrid + 1) - *(gr->ygrid);
   dy2 = *(gr->ygrid + gr->iy - 1) - *(gr->ygrid + gr->iy - 2);

   dx1 = *(gr->xgrid + 1) - *(gr->xgrid);
   dx2 = *(gr->xgrid + gr->ix - 1) - *(gr->xgrid + gr->ix - 2);

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
      
   if(tom == 'e'){
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

   if(form == 4){
      if(((NETCDF_INFO *)uu)->imiss || ((NETCDF_INFO *)vv)->imiss) {
         printf("****WARNING****, field has missing data value or fill value, do you want to continue, 'y' or 'n'\n\n");
         scanf("\n");
         if(getchar() == 'n') exit(1);
	 imiss = 'y';
      }
   }
   
   else{
      printf("Are there missing data values, delineated by some extreme value. 'y' or 'n'\n\n");
      scanf("\n");
      imiss=getchar();
   }

   if(imiss == 'y'){
     printf("What is the missing data value to test against?\n\n");
     scanf("%f", &mval);

     printf("How do you want to compare with missing value?\r\n"
            "'0' -- equality.                              \r\n"
            "'1' -- less than.                             \r\n"
            "'2' -- greater than.                          \n\n");
     scanf("%d", &icmp);
     
     printf("Whats missing value should be used for output?\n\n");
     scanf("%f", &mdval);
   }

   printf("What is the output file name for the vorticty?\n\n");
   scanf("%s", vorf);

   uap = (float *)calloc(dim, sizeof(float));
   mem_er((uap == NULL) ? 0: 1, dim * sizeof(float));  

   vap = (float *)calloc(dim, sizeof(float));
   mem_er((vap == NULL) ? 0: 1, dim * sizeof(float)); 

   vor=(float *)calloc(dim, sizeof(float));
   mem_er((vor == NULL) ? 0: 1, dim * sizeof(float));    

   printf("Use wrapround in longitude, 'y' or 'n'\n\n");
   scanf("\n");
   if(getchar() == 'y') nwr = 1;
   else nwr = 0;
   
   if(iper == 'n'){
       printf("****WARNING****, wrapround may not be valid for this data set.\n\n");
   }


/* open file for writing the vorticity field */

   fvor = open_file(vorf, "w");

   write_header(gr, fvor);

/* nfu and nfv are the same for ictype=3 */

   while(nfu <= maxe && nfv <= maxe && !ir){
  
       if(!nfu && !nfv){

/* read U */

          if((rf = read_field(uap, NULL, 1.0, uu, 1, 'n', 'n', '0', 'y'))) ir = 1;
       
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

        if(ir) continue;

        nfu += fr[0];
        nfv += fr[1];

        if(tom == 'g'){
           for(i=0; i < gr->iy; i++){
               coslat = *(coslt + i);
               for(j=0; j < gr->ix; j++){
                   *(uap + i * gr->ix + j) *= coslat;
               }
           }
        }

        itpfc = 1;

        for(i=0; i < gr->iy; i++){
	   
           if(tom == 'g') coslat = *(coslt + i);
	   
	   if(!i) dy = dy1;
	   else if (i == gr->iy - 1) dy = dy2;
	   else dy = *(gr->ygrid + i + 1) - *(gr->ygrid + i - 1);
	   
           for(j=0; j < gr->ix; j++){
	   
	      ii = i * gr->ix + j;
	      j1 = i * gr->ix + j - 1;
	      j2 = i * gr->ix + j + 1;
	      i1 = (i - 1) * gr->ix + j;
	      i2 = (i + 1) * gr->ix + j;
	      
/* printf("%d %d %d %d %d\n", ii, i1, i2, j1, j2); */
	      
	      dx = *(gr->xgrid + j + 1) - *(gr->xgrid + j - 1);
	      
	      if(!j){
	         if(!nwr){
		    j1 = ii;
		    dx = dx1;
		 }
		 else{
		    if(iper == 'y'){
		       j1 = i * gr->ix + gr->ix - 2;
		       dx = *(gr->xgrid + 1) - *(gr->xgrid) + *(gr->xgrid + gr->ix - 1) - *(gr->xgrid + gr->ix - 2);
		    }
		    else{
		       j1 = i * gr->ix + gr->ix - 1;
		       dx = *(gr->xgrid + 1) - *(gr->xgrid) + period - *(gr->xgrid + gr->ix - 1);
		       if(*(gr->xgrid) > TOLGRID) dx += *(gr->xgrid);      
		    }
		 }
		 
	      }
	      else if(j == gr->ix - 1){
	         if(!nwr){
		    j2 = ii;
		    dx = dx2;
		 }
		 else{
		    if(iper == 'y'){
		       j2 = i * gr->ix + 1;
		       dx = *(gr->xgrid + 1) - *(gr->xgrid) + *(gr->xgrid + gr->ix - 1) - *(gr->xgrid + gr->ix - 2);
		    }
		    else{
		       j2 = i * gr->ix;
		       dx = period - *(gr->xgrid + gr->ix - 1) + *(gr->xgrid + gr->ix - 1) - *(gr->xgrid + gr->ix - 2);
		       if(*(gr->xgrid) > TOLGRID) dx += *(gr->xgrid);	       
		    }	 
		 }
		 
	      }
	      
	      if(!i) i1 = ii;
	      else if(i == gr->iy - 1) i2 = ii; 
	      
	      imm = 0;
	      
	      if(imiss == 'y'){
	         if(missing(*(vap + j2), mval, icmp) || missing(*(vap + j1), mval, icmp) || 
		    missing(*(uap + i2), mval, icmp) || missing(*(uap + i1), mval, icmp)             ) imm = 1;
	      }
	      
	      if(imm) *(vor + ii) = mdval;
	      
	      else{

	         *(vor + ii) = ((*(vap + j2) - *(vap + j1)) / dx) - ((*(uap + i2) - *(uap + i1)) / dy);

                 if(tom == 'g') {
                    *(vor + ii) /= (coslat * EARTH_RADIUS_M * FP_PI);
                 }
                 else {
                    *(vor + ii) *= vorcs;
                    if(addvor == 'y') {*(vor + ii) += *(uap + ii) * *(tanr + ii);}
                 }    
	       
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

      if(form != 4) {
         close_file(uu, ufil);
         close_file(vv, vfil);
      }
      else {
         netcdf_close((NETCDF_INFO *)uu);
         netcdf_close((NETCDF_INFO *)vv);
      }

   }


   free(uap); free(vap); free(vor);
   free(coslt);
   free(tanr);

   return;

}
