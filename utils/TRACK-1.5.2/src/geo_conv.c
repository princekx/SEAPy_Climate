#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "st_fo.h"
#include "st_obj.h"
#include "grid.h"
#include "m_values.h"
#include "mem_er.h"

extern int tf;
extern int x1u, y1u;

extern GRID *gr;


void geo_conv(struct feature_pts *fp)

{

    double xr1, yr1;
    double s1, c1, s2, c2;

    if(tf == 3){

       xr1 = *(gr->xgrid + ((fp->x).ixy)+x1u-2) * FP_PI;
       yr1 = FP_PI2 - *(gr->ygrid + ((fp->y).ixy)+y1u-2) * FP_PI; 

    }    

    else {

       xr1 = (fp->x).xy * FP_PI;
       yr1 = FP_PI2 - (fp->y).xy * FP_PI;

    }

    if(yr1 < 0.) yr1 = 0.;


    sincos(xr1, &s1, &c1);
    sincos(yr1, &s2, &c2);

    fp->gwk[0] = s2 * c1;
    fp->gwk[1] = s2 * s1;
    fp->gwk[2] = c2;

    fp->gwky = 1;


    return;

}
