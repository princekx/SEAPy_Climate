#include <Stdio.h>
#include <Math.h>
#include "statistic.h"

/* function to evaluate the fisher density at a point */

double fisher_density(struct dpt *pt, struct dpt *dt, double sm, double lg, int nf)

{

   double dis;

   dis = pt->xdt * dt->xdt + pt->ydt * dt->ydt + pt->zdt * dt->zdt;

   if(nf) dis -= 1.0;

   if(dis > lg){

      dis *= sm;
 
      return exp(dis);

   }

   else return 0.0;


}
