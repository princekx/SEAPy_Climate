#include <Stdio.h>
#include "grid.h"

extern float period;
extern int x1u;
extern int pb;

extern GRID *gr;


/* function to find the real value of a grid point from its integer position.
   This takes account of periodic objects in the X-direction if in force.     */

float obj_xreal(int ptx)

{

   float xx=0.;

   if(pb == 'y') {

     if(ptx < 0) xx = *(gr->xgrid + ptx + gr->ix - 1) - period;
     else if (ptx >= gr->ix) xx = *(gr->xgrid + ptx - gr->ix + 1) + period;
     else xx = *(gr->xgrid + ptx);

     return xx;

   }

   else {

       if(ptx < 0 || ptx >= gr->ix) return GERROR;

       return *(gr->xgrid + ptx);

   }

}
