#include <Stdio.h>
#include <Math.h>
#include "m_values.h"

/* Mercator conformal projection */

void mercator(float *grid, int typ, int itp)

{

    if(itp){                                  /* forward */

       *grid *= FP_PI;

       if (typ == 'y'){

          *grid += FP_PI2;
          *grid = (float)log(tan((double)*grid / 2.0));       

       }

    }

    else{                                    /* backward */

       if(typ == 'y'){

         *grid = 2.0 * (float)atan(exp((double) *grid));
         *grid -= FP_PI2;

       }

       *grid /= FP_PI;

    }

    return;

}
