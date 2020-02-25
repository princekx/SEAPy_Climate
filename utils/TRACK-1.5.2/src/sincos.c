#include <Math.h>

/* function to compute sin and cosine of an angle in radians */

void sincos(double ang, double *sn, double *cn)

{


    *sn = sin(ang);
    *cn = cos(ang);

    return;

}
