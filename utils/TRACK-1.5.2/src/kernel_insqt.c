#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "sqt.h"

/* function to determine if a kernel is completely enclosed by its spherical
   triangle.                                                                 */


int kernel_insqt(LEAF *lf, VEC *gv, struct dpt *dtt, double sang)
{

    int intr=0;

    double norm=0.0;

    VEC *pa=NULL, *pb=NULL;
    VEC cpv;
    VEC dpp;

    if(sang < 0.){

       printf("***ERROR***, negative smoothing parameter in %s\n", __FILE__);
       exit(1);
    }

    dpp.x = dtt->xdt;
    dpp.y = dtt->ydt;
    dpp.z = dtt->zdt;

    pa = gv + lf->ivec[0];
    pb = gv + lf->ivec[1];
    crosp(pa, pb, &cpv);
    if(lf->ud == 'd') negate(&cpv);
    norm = sqrt(dotp(&cpv, &cpv));
    normv(&cpv, norm);
    if(dotp(&cpv, &dpp) >= sang) ++intr;

    pa = gv + lf->ivec[1];
    pb = gv + lf->ivec[2];
    crosp(pa, pb, &cpv);
    if(lf->ud == 'd') negate(&cpv);
    norm = sqrt(dotp(&cpv, &cpv));
    normv(&cpv, norm);
    if(dotp(&cpv, &dpp) >= sang) ++intr;

    pa = gv + lf->ivec[2];
    pb = gv + lf->ivec[0];
    crosp(pa, pb, &cpv);
    if(lf->ud == 'd') negate(&cpv);
    norm = sqrt(dotp(&cpv, &cpv));
    normv(&cpv, norm);
    if(dotp(&cpv, &dpp) >= sang) ++intr;

    intr = (intr == 3) ? 1 : 0;

    return intr;

}
