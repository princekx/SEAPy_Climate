#include <Stdio.h>
#include <Math.h>
#include "grid.h"

/* function to create an orographic mask from the surface geopotential */

#define   PTOL       1.0e-4
#define   TOLIVAL    1.0e-3

int qsearch(float * , float , int , int * , int * );
int mask_data(int * , float , float , int );

const float press[13]={1000.0, 900.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0, 200.0, 100.0, 50.0, 30.0, 10.0};
const float zz[13]={0.111, 0.988, 1.457, 3.012, 4.206, 5.574, 7.185, 9.164, 11.784, 16.180, 20.576, 23.849, 31.055};

extern GRID *gr;

int create_mask(int *imask, float *orog, float mpress, int dim)
{
   int i;
   int iret=0;

   float slp=0.0, icpt=0.0;
   float zt=0.0;

/* determine geopotential height */

   if(mpress > press[0]) zt = zz[0];
   else if(mpress < press[12]) zt = 100.0;
   else {
     for(i=0; i < 12; i++){
         if(fabs(press[i] - mpress) < PTOL) zt = zz[i];
         else if((press[i] - mpress) * (mpress - press[i+1]) >= 0.0){
            slp = (zz[i+1] - zz[i]) / (press[i+1] - press[i]);
            icpt = (zz[i] * press[i+1] - zz[i+1] * press[i]) / (press[i+1] - press[i]);
            zt = slp * mpress + icpt;
            break;
         }
     }
   }

   for(i=0; i < dim; i++){
      *(imask + i) = (*(orog + i) >= zt) ? 1 : 0;
      if(*(imask + i)) iret = 1;
   }

   return iret;
}

int mask_data(int *imsk, float x, float y, int imy)
{

   int ix1=0, ix2=0;
   int iy1=0, iy2=0;
   int ip1=-1, ip2=-1;

   int ii[4];

   int iret=0;

   if(!imy) return iret;

   ip1 = qsearch(gr->xgrid, x, gr->ix, &ix1, &ix2);
   ip2 = qsearch(gr->ygrid, y, gr->iy, &iy1, &iy2);

   if(ip1 >= 0 && ip2 >= 0){
      if(*(imsk + ip2 * gr->ix + ip1)) iret = 1;
   }
   else {
      ii[0] = *(imsk + iy1 * gr->ix + ix1);
      ii[1] = *(imsk + iy1 * gr->ix + ix2);
      ii[2] = *(imsk + iy2 * gr->ix + ix1);
      ii[3] = *(imsk + iy2 * gr->ix + ix2);

      if(ii[0] + ii[1] + ii[2] + ii[3] > 2) iret = 1;
   }

   return iret;
}

int qsearch(float *poss, float pps, int np, int *ip1, int *ip2)
{

   int i1=0, i2=0;
   int ih=0, iht=0, ig=-1;

   float p1, p2, ph;
   float dd=0;

   i1 = 0;
   i2 = np - 1;
   ih = (i1 + i2) / 2;
   p1 = *poss;
   p2 = *(poss + np - 1);
   ph = *(poss + ih);

   while(ih != iht) {

      dd = (p2 - pps) * (pps - ph);
      if(dd >= 0.0) {
         i1 = ih;
         p1 = ph;
      }
      else {
         i2 = ih;
         p2 = ph;
      }
      iht = ih;
      ih = (i1 + i2) / 2;
      ph = *(poss + ih);

   }

   *ip1 = i1;
   *ip2 = i2;

   if(fabs(p1 - pps) < TOLIVAL) ig = i1;
   else if(fabs(p2 - pps) < TOLIVAL) ig = i2;

   return ig;
}
