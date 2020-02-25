#include <Stdio.h>
#include <Math.h>
#include "m_values.h"

#define TYPE  0              /* TYPE=0 for Urmaev, TYPE=1 for Pavlov */

#if TYPE                     /* Pavlov, A0=a0, A2=a2*0.33, A4=a5*0.2 */

#define   A0     1.00000
#define   A2    -0.05103
#define   A4    -0.00534

#else                        /* Urmaev, A0=a0, A2=a2*0.33, A4=a5*0.2 */

#define   A0     0.92810
#define   A2     0.37143
#define   A4     0.00000

#endif

float ufunc(float g)
{

   double pw;

   pw = pow(g, 2.00);
   return A0 * g + g * (A2 * pw + A4 * pw * pw);

}

float secant(float (*ff)(float ), float , float , float );

void  urmaev(float *grid, int typ, int itp)

{

    double grd;

    float p1, p2;

    float (*f)(float );

    f = ufunc;

    if(itp){

      *grid *= FP_PI;

      if (typ == 'y'){

         grd = *grid;
         *grid = ufunc(grd);


      }

    }

    else {

      if(typ == 'y') {

        if(*grid > 0.) {p1 = 0.; p2 = FP_PI2;}
        else{p2 = 0.; p1 = -FP_PI2;}

        *grid = secant(f, p1, p2, *grid);

      }

      *grid /= FP_PI;

    }

    return;

}
