#include <Stdio.h>
#include <Math.h>
#include "m_values.h"

/* cylidrical stereographic projections */

#define  STLAT1   1.0        /* Braun cosine of st. lat. of 0 deg. */
#define  STLAT2   0.8660254  /* BSAM cosine of standard latitude of 30 deg. */
#define  STLAT3   0.7071067  /* Gall cosine of standard latitude of 45 deg. */

void stereog(float *grid, int typ, int itp)

{

    float stlat=STLAT1, st1;

    st1 = 1.0 + stlat;

    if(itp){

       *grid *= FP_PI;

       if(typ == 'x') *grid *= stlat;

       else if (typ == 'y'){

          *grid = st1 * tan((double)(*grid / 2.0));     

       }

    }

    else {

       if(typ == 'x') *grid /= stlat;

       else if(typ == 'y') *grid = 2.0 * (float)atan((double) *grid / st1);

       *grid /= FP_PI;

    }

    return;

}
