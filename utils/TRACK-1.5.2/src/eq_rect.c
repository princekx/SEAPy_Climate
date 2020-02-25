#include <Stdio.h>
#include "m_values.h"

#define STLAT1  0.7933533       /* cosine of the standard lat. of 37.5 deg. */
#define STLAT2  0.7313537       /* cosine of the standard lat. of 40.0 deg. */

/* equirectangular map transformation */


void eq_rect(float *grid, int typ, int itp)

{

    float stlat=STLAT2;

    if(itp){

       *grid *= FP_PI;

       if(typ == 'x') *grid *= stlat;

    }

    else{

       if(typ == 'x') *grid /= stlat;

       *grid /= FP_PI;

    }

    return;

}
      
