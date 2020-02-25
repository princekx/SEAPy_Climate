#include <Stdio.h>
#include <stdlib.h>

#include "grid.h"
#include "splice.h"

/* function to perform bi-linear interpolation taking account of missing values. */

int qsearch(float * , float , int , int * , int * );
int missing(float, float, int );

extern GRID *gr;

float bilinear_intrp(float *ap, float xx, float yy, float mval, int icmp)
{

   int ix1=0, ix2=0, iy1=0, iy2=0;

   float val=0.0;
   float p1=0.0, p2=0.0, p3=0.0, p4=0.0;
   float xd=0.0, yd=0.0;
   float pa=0.0, pb=0.0;

/* find neighbours */

   if(xx > *(gr->xgrid + gr->ix - 1) || xx < *(gr->xgrid)) return ADD_UNDEF;

   if(yy > *(gr->ygrid + gr->iy - 1) || yy < *(gr->ygrid)) return ADD_UNDEF;

   qsearch(gr->xgrid, xx, gr->ix, &ix1, &ix2);
   qsearch(gr->ygrid, yy, gr->iy, &iy1, &iy2);

   p1 = *(ap + iy1 * gr->ix + ix1);
   p2 = *(ap + iy1 * gr->ix + ix2);
   p3 = *(ap + iy2 * gr->ix + ix1);
   p4 = *(ap + iy2 * gr->ix + ix2);

   if(missing(p1, mval, icmp) || missing(p2, mval, icmp) || 
      missing(p3, mval, icmp) || missing(p4, mval, icmp)    ) return ADD_UNDEF;

   xd = (xx - *(gr->xgrid + ix1)) / (*(gr->xgrid + ix2) - *(gr->xgrid + ix1));
   yd = (yy - *(gr->ygrid + iy1)) / (*(gr->ygrid + iy2) - *(gr->ygrid + iy1));

   pa = p1 + xd * (p2 - p1);
   pb = p3 + xd * (p4 - p3);

   val = pa + yd * (pb - pa);
    
   return val;
}
