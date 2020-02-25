#include <Stdio.h>
#include <Math.h>
#include "statistic.h"

/* function to evaluate the linear density at a point */

double power_density(struct dpt *pt, struct dpt *dt, double sm, int ipow)

{

   double dis=0.0;

   dis = pt->xdt * dt->xdt + pt->ydt * dt->ydt + pt->zdt * dt->zdt;

   if(sm > 0.){

      dis *= sm;
      dis -= 1.0;

   }


   if(dis > 0.){

      switch(ipow){
         case 0:
           return 1.0;
           break;
         case 1:
           return dis;
           break;
         default:
           return pow(dis, (float)ipow);
      }
   }
   else return 0.0;

}
