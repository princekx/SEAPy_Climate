#include <Math.h>
#include "statistic.h"
#include "p_vecs.h"

/* function to evaluate non-isotropic kernel, un-normalized */


double non_iso_density(struct dpt *pt, struct dpt *dt, struct cvecs *dtpv, double sm, double beta, int mm, int *reg)

{

     float nn=(float)mm;

     double dis;
     double d1, d2;

     *reg = 0; 

     dis = sm * (pt->xdt * dt->xdt + pt->ydt * dt->ydt + pt->zdt * dt->zdt);

     d1 = pt->xdt * dtpv->p3[0] + pt->ydt * dtpv->p3[1] + pt->zdt * dtpv->p3[2];

     d1 *= d1;

     d2 = pt->xdt * dtpv->p2[0] + pt->ydt * dtpv->p2[1] + pt->zdt * dtpv->p2[2];

     d2 *= d2;

     dis += beta * (d1 - d2) - 1.0;

     if(dis > 0.0){

        *reg = 1;

        if(!mm) return 1.0;

        else if(mm == 1) return dis;
        
        else return pow(dis, nn);

     }

     else return 0.0;

}
