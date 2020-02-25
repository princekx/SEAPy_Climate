#include <Stdio.h>
#include <stdlib.h>
#include <reg_dat.h>
#include <mem_er.h>
#include "grid.h"

void error_interp(int , char * , int , int );


/* function to perform the setup for regular interpolation.
   converted from the smoopy.f code of DIERCKX.                   
   P. Dierckx, SIAM J. Numer. Analy. 19, 1286--1304, 1982         */

struct savedat *sd;
struct rspline *rs;

extern float xmn, ymn, xmx, ymx;
extern int x1u, x2u, y1u, y2u;

extern int delb, bs;

extern GRID *gr;

void smoopy_setup(int mem)

{

   int i, nmaxx, nmaxy;
   int xdm=x2u-x1u+1, ydm=y2u-y1u+1;
   int xd1=x1u-delb, xd2=x2u+delb, yd1=y1u-delb, yd2=y2u+delb;
   int err;

/* free previously assigned memory */

   if(mem){

      free(sd->nrdatx); 
      free(sd->nrdaty);
      free(rs->x); 
      free(rs->y);
      free(rs->z);
      return;

   }

/* assign spline degree data */

   rs->kx = KX;
   rs->ky = KY;
   rs->kx1 = KX + 1;
   rs->ky1 = KY + 1;
   rs->kx2 = KX + 2;
   rs->ky2 = KY + 2;
   rs->kmax = (rs->kx > rs->ky) ? rs->kx : rs->ky;
   rs->kmax2 = rs->kmax + 2;
   rs->k2max2 = 2* rs->kmax + 2;

/* define region of interest in the data */

   if(bs == 'y'){
      rs->sx = (xd1 > 0) ? xd1 : 1;
      rs->sy = (yd1 > 0) ? yd1 : 1;
      rs->ex = (xd2 < gr->ix) ? xd2 : gr->ix;
      rs->ey = (yd2 < gr->iy) ? yd2 : gr->iy;
      rs->xb = (double)*(gr->xgrid + rs->sx - 1);
      rs->xe = (double)*(gr->xgrid + rs->ex - 1);
      rs->yb = (double)*(gr->ygrid + rs->sy - 1);
      rs->ye = (double)*(gr->ygrid + rs->ey - 1);
      rs->mx = rs->ex - rs->sx + 1;
      rs->my = rs->ey - rs->sy + 1;
   }
   else{
      rs->xb = (double)xmn;
      rs->yb = (double)ymn;
      rs->xe = (double)xmx; 
      rs->ye = (double)ymx;
      rs->sx = x1u;
      rs->sy = y1u;
      rs->ex = x2u;
      rs->ey = y2u;
      rs->mx = xdm;
      rs->my = ydm;
    }

   nmaxx = rs->nxest = rs->mx + rs->kx1;
   nmaxy = rs->nyest = rs->my + rs->ky1;

   rs->mxy = rs->mx * rs->my;

   rs->nmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

   rs->ncof = rs->mxy;

/* assign storage to pointers */

   sd->nrdatx = (int * )calloc(rs->nmax, sizeof(int));
   mem_er((sd->nrdatx == NULL) ? 0 : 1, rs->nmax * sizeof(int));

   sd->nrdaty = (int * )calloc(rs->nmax, sizeof(int));
   mem_er((sd->nrdaty == NULL) ? 0 : 1, rs->nmax * sizeof(int));   
 
   rs->x = (double * )calloc(rs->mx, sizeof(double));
   mem_er((rs->x == NULL) ? 0 : 1, rs->mx * sizeof(double));

   rs->y = (double * )calloc(rs->my, sizeof(double));
   mem_er((rs->y == NULL) ? 0 : 1, rs->my * sizeof(double));

   rs->z = (double *)calloc(rs->mxy, sizeof(double));
   mem_er((rs->z == NULL) ? 0 : 1, rs->mxy * sizeof(double));


/* perform initial data check */

   if(             (rs->kx <= 0 || rs->ky <= 0)                                ||
               (rs->mx < rs->kx1 || rs->nxest < 2 * rs->kx1)                   ||
               (rs->my < rs->ky1 || rs->nyest < 2 * rs->ky1)                   ||
     (rs->xb > *(gr->xgrid + rs->sx - 1) || rs->xe < *(gr->xgrid + rs->ex - 1))||
     (rs->yb > *(gr->ygrid + rs->sy - 1) || rs->ye < *(gr->ygrid + rs->ey-1))      ){

       printf("***error***, file = %s, line = %d\n", __FILE__, __LINE__);
       printf("kx = %d <= 0 or ky = %d <= 0\n", rs->kx, rs->ky);
       printf("mx = %d < kx1 = %d or my = %d < ky1 = %d\n", rs->mx, rs->kx1, rs->my, rs->ky1);
       printf("nxest = %d, nyest = %d\n", rs->nxest, rs->nyest);
       printf("xb = %f, xe = %f, yb = %f, ye = %f\n", rs->xb, rs->xe, rs->yb, rs->ye);
       printf("defined region incompatable with data\n\n");
       error_interp(10, __FILE__, __LINE__, 1);



   }

/* initialize grid arrays */

   for(i=0; i < rs->mx; i++) *(rs->x + i) = (double)*(gr->xgrid + rs->sx + i - 1);

/*   for(i=0; i < rs->mx; i++) printf("%f \n", *(rs->x + i)); */

   for(i=0; i < rs->my; i++) *(rs->y + i) = (double)*(gr->ygrid + rs->sy + i - 1);

/*   for(i=0; i < rs->my; i++) printf("%f \n", *(rs->y + i)); */

/* this last check is not strictly neccesary but does'nt hurt */

   err = 0;

   for(i=1; i < rs->mx; i++) 

      if(*(rs->x + i - 1) >= *(rs->x + i)) {err = 1; break;}


   if(err == 1 ) {

      printf("***error***, file = %s, line = %d\n", __FILE__, __LINE__);
      printf("***error*** in x-data %f %f\n\n", *(rs->x + i - 1), *(rs->x + i));
      error_interp(10, __FILE__, __LINE__, 1);

   }

   err = 0;

   for(i=1; i < rs->my; i++)

      if(*(rs->y + i - 1) >= *(rs->y + i)) {err = 1; break;}

   if(err == 1 ){

     printf("***error***, file = %s, line = %d\n", __FILE__, __LINE__);
     printf("***error*** in x-data %f %f\n\n", *(rs->y + i - 1), *(rs->y + i));
     error_interp(10, __FILE__, __LINE__, 1); 

   }  
 
   return;

}
