#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IPT  0

#define   TOLIVAL    1.0e-3

/* test if point is over land */

int qsearch(float * , float , int , int * , int * );

int orog_test(float *glng, float *glat, int lglng, int lglat, int *ilms, float xpos, float ypos)
{
    int ipt=IPT;
    int ix1=0, ix2=0, iy1=0, iy2=0;

    int ip1=0, ip2=0, ip3=0, ip4=0;

    if(xpos > *(glng + lglng - 1) || xpos < *glng) return 0;

    if(ypos > *(glat + lglat - 1) || ypos < *glat) return 0;

    qsearch(glng, xpos, lglng, &ix1, &ix2);
    qsearch(glat, ypos, lglat, &iy1, &iy2);

    ip1 = *(ilms + iy1 * lglng + ix1);
    ip2 = *(ilms + iy1 * lglng + ix2);
    ip3 = *(ilms + iy2 * lglng + ix1);
    ip4 = *(ilms + iy2 * lglng + ix2);

    if(ipt) {
      if(ip1 && ip2 && ip3 && ip4) return 1;   /* all points land points */
      else return 0;
    }
    else {
      if(ip1 || ip2 || ip3 || ip4) return 1;   /* any points are land points */
      else return 0;
    }
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







