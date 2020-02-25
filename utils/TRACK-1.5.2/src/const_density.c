#include <Stdio.h>
#include "statistic.h"

/* function to evaluate the linear density at a point */

double const_density(struct dpt *pt, struct dpt *dt, double sm)

{

   double dis;

   dis = pt->xdt * dt->xdt + pt->ydt * dt->ydt + pt->zdt * dt->zdt;

   if(sm > 0.){

      dis *= sm;
      dis -= 1.0;

   }

   if(dis > 0.) return 1.0;

   else return 0.0;


}
