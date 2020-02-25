#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include "mem_er.h"
#include "grid.h"
#include "proj.h"

/* function to produce grid for application to a conic projection */

extern GRID *gr1, *gr2;
extern CNTRY *cm1, *cm2;

extern GRID *gt;

extern CNTRY *ct;

extern FILE *finit;
extern int init;

int search_area_check(int , int , int );
void conic_grid_assign(GRID * , CNTRY * , float , float , float , float , int , int , int , int , int , int , int );


void conic_grid(int ptp)

{

   int prty;
   int iaz=0;
   int ift=1;

   static char *prname[] = {"Albers Equal Area", "Lambert Conformal",
                            "Equidistant Conic"};
			    
   printf("Is data already on a conic grid, '0' for no or '1' for yes.\n\n");
   
   
   if(init)fscanf(finit, "%d", &iaz);
   else {
     scanf("%d", &iaz);
     if(finit) fprintf(finit, "%d\n", iaz);
   }
  

   if(iaz <0 || iaz > 1){
      printf("****ERROR****, incorrect input.\n\n");
      exit(1);
   }
   
   gt->iaz = iaz;

   if(ptp < 0){

      printf("which conic projection is required,   \r\n"
             "available projections are:          \r\n\n"
             "Albers Equal Area   input '1'         \r\n"
             "Lambert Conformal   input '2'         \r\n"
             "Equidistant Conic   input '3'         \n\n");

      if(init){
         fscanf(finit, "%d", &prty);
      }
      else{
         scanf("%d", &prty);
         if(finit)fprintf(finit, "%d\n", prty);
      }

   }

   else prty = ptp;

   if(gr1){free(gr1->xgrid); free(gr1->ygrid); free(gr1);}
   if(cm1){free(cm1->cmxg); free(cm1->cmyg); free(cm1->cmi); free(cm1);}

   gr1 = (GRID * )malloc_initl(sizeof(GRID));
   mem_er((gr1 == NULL) ? 0 : 1, sizeof(GRID));

   gr1->prgr = gt->prgr;
   gt->prty = prty;

   cm1 = (CNTRY * )malloc_initl(sizeof(CNTRY));
   mem_er((cm1 == NULL) ? 0 : 1, sizeof(CNTRY));

   printf("what is the origin for the projection, input a (lat, long) pair.\n\n");

   if(ptp < 0){

      if(init){
         fscanf(finit, "%f %f", &gr1->alat, &gr1->alng);
      }
      else{
         scanf("%f %f", &(gr1->alat), &(gr1->alng));
         if(finit)fprintf(finit, "%f\n%f\n", gr1->alat, gr1->alng);
      }
      
      
      printf("What are the two standard parallels, sp1, sp2?\r\n"
	     "If only one standard parallel is required use the same value for both.\n\n");
      
      if(init){
         fscanf(finit, "%f %f", &(gr1->sp1), &(gr1->sp2));
      }
      else{
         scanf("%f %f", &(gr1->sp1), &(gr1->sp2));
         if(finit)fprintf(finit, "%f %f\n", gr1->sp1, gr1->sp2);
      }


   }

   else {

      gr1->alat = gt->alat;
      gr1->alng = gt->alng;
      gr1->sp1 = gt->sp1;
      gr1->sp2 = gt->sp2;

   }

   gr1->prty = prty;

   strncpy(gt->prty_nm, prname[prty-1], MXPRCH);
   strncpy(gr1->prty_nm, prname[prty-1], MXPRCH);
   strncpy(gr1->prgr_nm, gt->prgr_nm, MXPRCH);

   if(gr1->alng > 360.) gr1->alng -= 360.;
   else if(gr1->alng < 0.) gr1->alng += 360.;

   printf("input the region required, in grid points                    \n"
          "xng1, xng2 (xng1 < xng2), ylt1, ylt2 (ylt1 < ylt2).          \n"
          "The current region is %d (X) by %d (Y)\n\n", gt->ix, gt->iy);

   if(init){
      fscanf(finit, "%d %d %d %d", &gr1->ox1u, &gr1->ox2u, &gr1->oy1u, &gr1->oy2u);
   }
   else {
      scanf("%d %d %d %d", &gr1->ox1u, &gr1->ox2u, &gr1->oy1u, &gr1->oy2u);
      if(finit)fprintf(finit, "%d\n%d\n%d\n%d\n", gr1->ox1u, gr1->ox2u, gr1->oy1u, gr1->oy2u);
   }

   search_area_check(gr1->ox1u, gr1->ox2u, gt->ix);
   search_area_check(gr1->oy1u, gr1->oy2u, gt->iy);

/* find region extent in projected space */

   conic_grid_assign(gr1, cm1, gr1->alat, gr1->alng, gr1->sp1, gr1->sp2, gr1->ox1u, gr1->ox2u, gr1->oy1u, gr1->oy2u, prty, iaz, ift);
      
   proj_report(gr1->prty, gr1->prgr);

   return;

}


void conic_grid_assign(GRID *ggr, CNTRY *cmm, float alat, float alng, float sp1, float sp2, int xng1, int xng2, int ylt1, int ylt2, int prty, int iaz, int ift)

{

      int i, j;
      int ifp = 0;
      int igtyp=0;

      float xx, yy, yyy;
      float lngmn=0., lngmx=0., latmn=0., latmx=0.;
      float xlen, ylen;
      float xst, yst;
     
      
      if(!iaz){

         for(i=ylt1-1; i < ylt2; i++){
             yyy = *(gt->ygrid + i);
             for(j=xng1-1; j < xng2; j++){

                 yy = yyy;
                 xx = *(gt->xgrid + j);
		 
                 if(conic(&yy, &xx, alat, alng, sp1, sp2, 1, prty, &ift)){

                    if(ifp){

                       if(xx > lngmx) lngmx = xx;
                       else if (xx < lngmn) lngmn = xx;

                       if(yy > latmx) latmx = yy;
                       else if (yy < latmn) latmn = yy;

                    }

                    else{lngmn = lngmx = xx; latmn = latmx = yy; ifp = 1;}
		    
                 }

             }


         }

         xlen = lngmx-lngmn;
         ylen = latmx-latmn;
	 
	 if(xlen <= 0.0 || ylen <= 0.0){
	    printf("****ERROR****, domain range invalid in either X or Y, %f %f\n\n", xlen, ylen);
	    exit(1);	 
	 }
	 
	 printf("Define number of points in X and Y or the grid spacing, '0' for number, '1' for spacing.\n\n");
	 if(init){	 
	    fscanf(finit, "%d", &igtyp);
	 }
	 else {
	    scanf("%d", &igtyp);
	    if(finit)fprintf(finit, "%d\n", igtyp);
	 }

	 if(igtyp < 0 || igtyp > 1){
	    printf("****WARNING****, icorrect identifier for grid, defaulting to choosing number of points.\n\n");
	    igtyp = 0;	 
	 }
	 
	 if(!igtyp){

            printf("this region has an extent of %f starting at %f (X), and                    \r\n"
                   "                             %f starting at %f (Y) in the projected space, \r\n"
                   "how many new grid points are required in the X-direction > 1 and the       \r\n"
                   "Y-direction > 1\n\n", xlen, lngmn, ylen, latmn);

            if(init){
               fscanf(finit, "%d %d", &(ggr->ix), &(ggr->iy));
            }
            else{
               scanf("%d %d", &(ggr->ix), &(ggr->iy));
               if(finit)fprintf(finit, "%d\n%d\n", ggr->ix, ggr->iy);
            }

            xst = xlen / (ggr->ix - 1);
            yst = ylen / (ggr->iy - 1);
	 }
	 
	 else {
	    printf("What is the grid spacing in X and Y?\n\n");
	    scanf("%f %f", &xst, &yst);
	    ggr->ix = (xlen / xst) + 1;
	    ggr->iy = (ylen / yst) + 1;	 
	 }

/* assign memory, compute new grid points */

         ggr->xgrid = (float * )calloc(ggr->ix, sizeof(float));
         mem_er((ggr->xgrid == NULL) ? 0 : 1, (ggr->ix) * sizeof(float));

         ggr->ygrid = (float * )calloc(ggr->iy, sizeof(float));
         mem_er((ggr->ygrid == NULL) ? 0 : 1, (ggr->iy) * sizeof(float));

         for(i=0; i < ggr->ix; i++) *(ggr->xgrid + i) = lngmn + xst * i;

         for(i=0; i < ggr->iy; i++) *(ggr->ygrid + i) = latmn + yst * i;
      
      }
      
      else {
      
         ggr->ix = gt->ix;
	 ggr->iy = gt->iy;
         ggr->xgrid = (float * )calloc(ggr->ix, sizeof(float));
         mem_er((ggr->xgrid == NULL) ? 0 : 1, (ggr->ix) * sizeof(float));

         ggr->ygrid = (float * )calloc(ggr->iy, sizeof(float));
         mem_er((ggr->ygrid == NULL) ? 0 : 1, (ggr->iy) * sizeof(float));
	 
	 memcpy(ggr->xgrid, gt->xgrid, ggr->ix * sizeof(float));
	 memcpy(ggr->ygrid, gt->ygrid, ggr->iy * sizeof(float));
	 
	 lngmn = lngmx = *(ggr->xgrid);
	 for(i=0; i < ggr->ix; i++){
	     if(*(ggr->xgrid + i) < lngmn) lngmn = *(ggr->xgrid + i);
	     else if (*(ggr->xgrid + i) > lngmx) lngmx = *(ggr->xgrid + i);
	 }
	 
	 latmn = latmx = *(ggr->ygrid);	 
	 for(i=0; i < ggr->iy; i++){
	     if(*(ggr->ygrid + i) < latmn) latmn = *(ggr->ygrid + i);
	     else if (*(ggr->ygrid + i) > latmx) latmx = *(ggr->ygrid + i);
	 }	 
      
      }

      if(ggr->ix <= 0 || ggr->iy <= 0) {ggr->xgrid = NULL; ggr->ygrid = NULL;}


/* assign memory and new country map */

      if(ct->dcm){

         cmm->dcm = ct->dcm;

         cmm->cmxg = (float * )calloc(cmm->dcm, sizeof(float));
         mem_er((cmm->cmxg == NULL) ? 0 : 1, (cmm->dcm) * sizeof(float));

         cmm->cmyg = (float * )calloc(cmm->dcm, sizeof(float));
         mem_er((cmm->cmyg == NULL) ? 0 : 1, (cmm->dcm) * sizeof(float));

         cmm->cmi = (int * )calloc(cmm->dcm, sizeof(int));
         mem_er((cmm->cmi == NULL) ? 0 : 1, (cmm->dcm) * sizeof(int));

         for(i=0; i < ct->dcm; i++){

             xx = *(ct->cmxg + i);
             yy = *(ct->cmyg + i);

             if(conic(&yy, &xx, alat, alng, sp1, sp2, 1, prty, &ift)){

                if((xx-lngmn)*(lngmx-xx) > 0. && (yy-latmn)*(latmx-yy) > 0.){

                   *(cmm->cmxg + i) = xx;
                   *(cmm->cmyg + i) = yy;
                   *(cmm->cmi + i) = *(ct->cmi + i);

                }

                else *(cmm->cmi + i) = - *(ct->cmi + i);
            

             }

             else *(cmm->cmi + i) = - *(ct->cmi + i);


         }


      }

      return;

}
