#include <Stdio.h>
#include "statistic.h"

/* function to evaluate the quadratic density at a point */

double quad_density(struct dpt *pt, struct dpt *dt, double sm)

{

   double dis, ds1, ds;

   dis = pt->xdt * dt->xdt + pt->ydt * dt->ydt + pt->zdt * dt->zdt;
   ds = dis;

   if(sm > 0.){

      dis *= sm;

      ds = dis - 1.0;

      ds1 = dis * dis - 1.0;

   }

   else ds1 = ds * ds;
 
   if(ds > 0.)return ds1;

   else return 0.0;


}
