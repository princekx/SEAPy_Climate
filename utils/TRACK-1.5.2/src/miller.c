#include <Stdio.h>
#include <Math.h>
#include "m_values.h"

/* Miller projections */

#define  MILLER1  1.25
#define  MILLER2  1.50

void miller(float *grid, int typ, int itp)

{

    float c=MILLER1;

    if(itp){

       *grid *= FP_PI;

       if (typ == 'y'){

          *grid /= (2.0 * c);
          *grid += FP_PI4;
          *grid = c * (float)log(tan((double)*grid));       

       }

    }

    else{

       if(typ == 'y')

          *grid = ((float) atan(exp((double) *grid / c)) - FP_PI4) * (2.0 * c);

       *grid /= FP_PI;

    }

    return;

}
