#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"


void hprint(struct image *hierarchy[], int level, int xxdim, int yydim)

{

     int i, j, k, n, xdim, ydim, plus;

     xdim=xxdim;
     ydim=yydim;

     for(k=0; k<level; k++) {
         for(j=0; j<ydim; j++) {

             for(i=0; i<xdim; i++) {

                 plus=j*xdim+i;

                 printf("%3d ", (hierarchy[k]+plus)->label->name);

             }
             printf("\n");
          }
          xdim=xdim/2;
          ydim=ydim/2;
          n=xdim*ydim;
          if(xdim == 1 && n != 1) xdim=2;
          if(ydim == 1 && n != 1) ydim=2; 
          printf("\n\n");
     }

     return;

}
