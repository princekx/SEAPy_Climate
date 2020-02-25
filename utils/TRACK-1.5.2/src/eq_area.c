#include <Stdio.h>
#include <Math.h>
#include "m_values.h"

#define  STLAT1 1.0        /* Lambert cosine of standard latitude of 0 deg.*/
#define  STLAT2 0.8660254  /* Behrmann cosine of standard latitude of 30 deg. */
#define  STLAT3 0.7071067  /* Peters cosine of standard latitude of 45 deg. */

/* equal area projection */

void eq_area(float *grid, int typ, int itp)

{

    float stlat=STLAT1;

    if(itp){

       *grid *= FP_PI;

       if(typ == 'x'){

          *grid *= stlat;

       }

       else if (typ == 'y'){

          *grid = (float)sin((double)*grid) / stlat;

       }

    }

    else{

       if(typ == 'x') *grid /= stlat;

       else if(typ == 'y') *grid = (float)asin((double)(*grid * stlat));

       *grid /= FP_PI;

    }

    return;

}
