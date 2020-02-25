#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"
#include "grid.h"

/* function to store the subsets of data defined by the user */

extern int x1u, x2u, y1u, y2u;         /* search area grid numbers */

extern GRID *gr;

void arrayd(struct image *ip, float *ap, float *tav, FILE *fthro, float th, int frnum, int xdim, int ttav)

{

     int i, j;
     int plus;

     float apt;
     float thresh, atv;

     struct image *pp;

     thresh = th;

/*     fprintf(fthro, "FRAME %6d \n", frnum); */

     for(j=y1u-1; j<=y2u-1; j++){
    
         for(i=x1u-1; i<= x2u-1; i++){

             apt = *(ap + j * gr->ix + i);
             if(ttav == 't') {
               atv = *(tav + j * gr->ix + i); 
               thresh = (atv >= th) ? atv : th;
             }

             plus=(j-y1u+1)*xdim+i-x1u+1;

             pp=ip+plus;

             if(apt - thresh >= 0.0) pp->pval = 1;

             else  pp->pval = 0;

             pp->nw=pp->ne=pp->sw=pp->se=NULL;

/*             fprintf(fthro,"%e ", *apt); */

 
         }

/*         fprintf(fthro,"\n"); */

     }

/*     fprintf(fthro, "\n"); */

     return;

}
