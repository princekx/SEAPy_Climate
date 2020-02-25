#include <Stdio.h>
#include <proj.h>
#include "grid.h"


/* projection transformation routine */

extern float xmn, ymn, xmx, ymx;
extern float period;

extern GRID *gr;
extern CNTRY *cm;

void trans_proj(PRFP proj, int typ)

{

   int i;

   (*proj)(&xmn, 'x', typ);
   (*proj)(&xmx, 'x', typ);
   (*proj)(&ymn, 'y', typ);
   (*proj)(&ymx, 'y', typ);

   for(i=0; i < gr->ix; i++)(*proj)(gr->xgrid+i, 'x', typ);
   for(i=0; i < gr->iy; i++)(*proj)(gr->ygrid+i, 'y', typ);
   for(i=0; i < cm->dcm; i++){
      (*proj)(cm->cmxg+i, 'x', typ);
      (*proj)(cm->cmyg+i, 'y', typ);
   }

   period = *(gr->xgrid + gr->ix - 1) - *(gr->xgrid);

   return;

}
