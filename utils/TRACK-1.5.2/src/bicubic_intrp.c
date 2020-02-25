#include <Stdio.h>
#include <stdlib.h>

#include "grid.h"
#include "splice.h"

/* function to perform bi-cubic interpolation taking account of missing values. */

/* macros for cubic splice function evaluations */

#define  f2(X) ((X * X - 1.0) * (X - 2.0) * 0.5)
#define  f3(X) (-(X + 1.0) * (X - 2.0) * X * 0.5)
#define  f4(X) ((X * X - 1.0) * X  / 6.0)

int qsearch(float * , float , int , int * , int * );
int missing(float, float, int );

extern GRID *gr;
extern int pb;

float bicubic_intrp(float *ap, float xx, float yy, float mval, int icmp)
{

   int ix1=0, ix2=0, iy1=0, iy2=0;

   float val=0.0;
   float p0[4], p1[4], p2[4], p3[4];
   float xd, y3[4], yd[3];
   float pp[4];

/* find neighbours */

   if(xx > *(gr->xgrid + gr->ix - 1) || xx < *(gr->xgrid)) return ADD_UNDEF;

   if(yy > *(gr->ygrid + gr->iy - 1) || yy < *(gr->ygrid)) return ADD_UNDEF;

   qsearch(gr->xgrid, xx, gr->ix, &ix1, &ix2);
   qsearch(gr->ygrid, yy, gr->iy, &iy1, &iy2);

   if(pb == 'n'){
     if(ix1 - 1 < 0) return ADD_UNDEF;
     else if(ix2 + 1 >= gr->ix) return ADD_UNDEF;
   }

   if(iy1 - 1 < 0) return ADD_UNDEF;
   else if(iy2 + 1 >= gr->iy) return ADD_UNDEF;

   if(pb == 'y'){
     if(ix1 - 1 < 0) {
       p0[0] = *(ap + (iy1 - 1) * gr->ix + gr->ix - 2);
       p1[0] = *(ap + iy1 * gr->ix + gr->ix - 2);
       p2[0] = *(ap + iy2 * gr->ix + gr->ix - 2);
       p3[0] = *(ap + (iy2 + 1) * gr->ix + gr->ix - 2);
     }
     else{
       p0[0] = *(ap + (iy1 - 1) * gr->ix + ix1 - 1);
       p1[0] = *(ap + iy1 * gr->ix + ix1 - 1);
       p2[0] = *(ap + iy2 * gr->ix + ix1 - 1);
       p3[0] = *(ap + (iy2 + 1) * gr->ix + ix1 - 1);
     }
     if(ix2 + 1 >= gr->ix){
       p0[3] = *(ap + (iy1 - 1) * gr->ix + 1);
       p1[3] = *(ap + iy1 * gr->ix + 1);
       p2[3] = *(ap + iy2 * gr->ix + 1);
       p3[3] = *(ap + (iy2 + 1) * gr->ix + 1);
     }
     else {
       p0[3] = *(ap + (iy1 - 1) * gr->ix + ix2 + 1);
       p1[3] = *(ap + iy1 * gr->ix + ix2 + 1);
       p2[3] = *(ap + iy2 * gr->ix + ix2 + 1);
       p3[3] = *(ap + (iy2 + 1) * gr->ix + ix2 + 1);
     }
   }
   else{
       p0[0] = *(ap + (iy1 - 1) * gr->ix + ix1 - 1);
       p1[0] = *(ap + iy1 * gr->ix + ix1 - 1);
       p2[0] = *(ap + iy2 * gr->ix + ix1 - 1);
       p3[0] = *(ap + (iy2 + 1) * gr->ix + ix1 - 1);
       p0[3] = *(ap + (iy1 - 1) * gr->ix + ix2 + 1);
       p1[3] = *(ap + iy1 * gr->ix + ix2 + 1);
       p2[3] = *(ap + iy2 * gr->ix + ix2 + 1);
       p3[3] = *(ap + (iy2 + 1) * gr->ix + ix2 + 1);
   }

   p0[1] = *(ap + (iy1 - 1) * gr->ix + ix1);
   p0[2] = *(ap + (iy1 - 1) * gr->ix + ix2);
   p1[1] = *(ap + iy1 * gr->ix + ix1);
   p1[2] = *(ap + iy1 * gr->ix + ix2);
   p2[1] = *(ap + iy2 * gr->ix + ix1);
   p2[2] = *(ap + iy2 * gr->ix + ix2);
   p3[1] = *(ap + (iy2 + 1) * gr->ix + ix1);
   p3[2] = *(ap + (iy2 + 1) * gr->ix + ix2);

   if(missing(p0[0], mval, icmp) || missing(p0[1], mval, icmp) || 
      missing(p0[2], mval, icmp) || missing(p0[3], mval, icmp)    ) return ADD_UNDEF;
   else if(missing(p1[0], mval, icmp) || missing(p1[1], mval, icmp) || 
      missing(p1[2], mval, icmp) || missing(p1[3], mval, icmp)    ) return ADD_UNDEF;
   else if(missing(p2[0], mval, icmp) || missing(p2[1], mval, icmp) || 
      missing(p2[2], mval, icmp) || missing(p2[3], mval, icmp)    ) return ADD_UNDEF;
   else if(missing(p3[0], mval, icmp) || missing(p3[1], mval, icmp) || 
      missing(p3[2], mval, icmp) || missing(p3[3], mval, icmp)    ) return ADD_UNDEF;

   xd = (xx - *(gr->xgrid + ix1)) / (*(gr->xgrid + ix2) - *(gr->xgrid + ix1));

   y3[0] = *(gr->ygrid + iy1 - 1);
   y3[1] = *(gr->ygrid + iy1);
   y3[2] = *(gr->ygrid + iy2);
   y3[3] = *(gr->ygrid + iy2 + 1);

   yd[0] = (yy - y3[0]) * (yy - y3[2]) * (yy - y3[3]) / ((y3[1] - y3[0]) * (y3[1] - y3[2]) * (y3[1] - y3[3]));
   yd[1] = (yy - y3[0]) * (yy - y3[1]) * (yy - y3[3]) / ((y3[2] - y3[0]) * (y3[2] - y3[1]) * (y3[2] - y3[3]));
   yd[2] = (yy - y3[0]) * (yy - y3[1]) * (yy - y3[2]) / ((y3[3] - y3[0]) * (y3[3] - y3[1]) * (y3[3] - y3[2]));   


   pp[0] = p0[1] + xd * (p0[2] - p0[1]);
   pp[1] = p1[0] + f2(xd) * (p1[1] - p1[0]) + f3(xd) * (p1[2] - p1[0]) + f4(xd) * (p1[3] - p1[0]);
   pp[2] = p2[0] + f2(xd) * (p2[1] - p2[0]) + f3(xd) * (p2[2] - p2[0]) + f4(xd) * (p2[3] - p2[0]);
   pp[3] = p3[1] + xd * (p3[2] - p3[1]);

   val = pp[0] + yd[0] * (pp[1] - pp[0]) + yd[1] * (pp[2] - pp[0]) + yd[2] * (pp[3] - pp[0]);
    
   return val;
}
