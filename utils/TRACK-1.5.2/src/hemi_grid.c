#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include "mem_er.h"
#include "grid.h"
#include "proj.h"

/* function to produce grid for application to a hemisphere */

extern GRID *gr1, *gr2;
extern CNTRY *cm1, *cm2;

extern GRID *gt;

extern CNTRY *ct;

extern FILE *finit;
extern int init;

int search_area_check(int , int , int );
void azimuthal_grid_assign(GRID * , CNTRY * , float , float , int , int , int , int , int , int );


void hemi_grid(int ptp)

{

   int prty;
   int wh=0;
   int iaz=0;

   static char *prname[] = {"Orthographic", "Stereographic",
                            "Gnomnic", "Lambert Equal Area", 
                            "Azimuthal Equal Distance"};
			    
   printf("Is data already on an azimuthal grid, '0' for no or '1' for yes.\n\n");
   
   
   if(init)fscanf(finit, "%d", &iaz);
   else {
     scanf("%d", &iaz);
     if(finit) fprintf(finit, "%d\n", iaz);
   }

   if(iaz <0 || iaz > 1){
      printf("****ERROR****, incorrect input.\n\n");
      exit(1);
   }

   if(!iaz) {
      printf("are both hemispheres required or just one,\n"
             " 'b' (for both) or 's' (for single).      \n\n");


      if(init){
         fscanf(finit, "\n");
         wh = getc(finit);
      }
      else {

        scanf("\n");
        wh = getchar();
        if(finit)fprintf(finit, "%c\n", wh);

      }
   }
   else {wh = 's'; gt->iaz = iaz;}

   if(wh == 'b' || wh == 's'){

      if(ptp < 0){

         printf("which azimuthal projection is required, \r\n"
                " available projections are:             \r\n\n"
                "Orthographic    input '1'               \r\n"
                "Stereographic   input '2'               \r\n"
                "Gnomnic         input '3'               \r\n"
                "Lambert EA      input '4'               \r\n"
                "Azimuthal ED    input '5'               \n\n");

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

      printf("what is the origin for the 1st projection, input a (lat, long) pair.\n\n");

      if(ptp < 0){

         if(init){
            fscanf(finit, "%f %f", &gr1->alat, &gr1->alng);
         }
         else{
            scanf("%f %f", &gr1->alat, &gr1->alng);
            if(finit)fprintf(finit, "%f\n%f\n", gr1->alat, gr1->alng);
         }

      }

      else {

         gr1->alat = gt->alat;
         gr1->alng = gt->alng;

      }

      gr1->prty = prty;

      strncpy(gt->prty_nm, prname[prty-1], MXPRCH);
      strncpy(gr1->prty_nm, prname[prty-1], MXPRCH);
      strncpy(gr1->prgr_nm, gt->prgr_nm, MXPRCH);

      if(gr1->alng > 360.) gr1->alng -= 360.;
      else if(gr1->alng < 0.) gr1->alng += 360.;

      printf("input the regions of the hemisphere required, in grid points \n"
             "xng1, xng2 (xng1 < xng2), ylt1, ylt2 (ylt1 < ylt2).          \n"
             "The current region is %d (X) by %d (Y)\n\n", gt->ix, gt->iy);

      printf("1st hemisphere\n");

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

      azimuthal_grid_assign(gr1, cm1, gr1->alat, gr1->alng, gr1->ox1u, gr1->ox2u, gr1->oy1u, gr1->oy2u, prty, iaz);


      if(wh == 'b'){

         if(gr2){free(gr2->xgrid); free(gr2->ygrid); free(gr2);}
         if(cm2){free(cm2->cmxg); free(cm2->cmyg); free(cm2->cmi); free(cm2);}

         gr2 = (GRID * )malloc_initl(sizeof(GRID));
         mem_er((gr2 == NULL) ? 0 : 1, sizeof(GRID));

         cm2 = (CNTRY * )malloc_initl(sizeof(CNTRY));
         mem_er((cm2 == NULL) ? 0 : 1, sizeof(CNTRY));

         gr2->prgr = gt->prgr;

         printf("2nd hemisphere\n\n");


         printf("what is the origin for the 2nd projection, input a (lat, long) pair.\n\n");


         if(init){
            fscanf(finit, "%f %f", &gr2->alat, &gr2->alng);
         }
         else{
            scanf("%f %f", &gr2->alat, &gr2->alng);
            if(finit)fprintf(finit, "%f\n%f\n", gr2->alat, gr2->alng);
         }


         if(gr2->alng > 360.) gr2->alng -= 360.;
         else if(gr2->alng < 0.) gr2->alng += 360.;

         gr2->prty = gr1->prty;
         strncpy(gr2->prty_nm, gr1->prty_nm, MXPRCH);
         strncpy(gr2->prgr_nm, gr1->prgr_nm, MXPRCH);

         printf("input new region, in grid points\n"
                "xng1, xng2 (xng1 < xng2), ylt1, ylt2 (ylt1 < ylt2)\n\n");

         if(init){
            fscanf(finit, "%d %d %d %d", &gr2->ox1u, &gr2->ox2u, &gr2->oy1u, &gr2->oy2u);
         }
         else {
            scanf("%d %d %d %d", &gr2->ox1u, &gr2->ox2u, &gr2->oy1u, &gr2->oy2u);
            if(finit)fprintf(finit, "%d\n%d\n%d\n%d\n", gr2->ox1u, gr2->ox2u, gr2->oy1u, gr2->oy2u);
         }


         search_area_check(gr1->ox1u, gr1->ox2u, gt->ix);
         search_area_check(gr1->oy1u, gr1->oy2u, gt->iy);

         azimuthal_grid_assign(gr2, cm2, gr2->alat, gr2->alng, gr2->ox1u, gr2->ox2u, gr2->oy1u, gr2->oy2u, prty, iaz);

      }


   }

   else {

      printf("***ERROR***, incorrect specifier for hemi-spherical projection\n");
      exit(1);

   }

   proj_report(gr1->prty, gr1->prgr);

   return;

}


void azimuthal_grid_assign(GRID *ggr, CNTRY *cmm, float alat, float alng, int xng1, int xng2, int ylt1, int ylt2, int prty, int iaz)

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

                 if(azimuthal(&yy, &xx, alat, alng, 1, prty)){

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

             if(azimuthal(&yy, &xx, alat, alng, 1, prty)){

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
