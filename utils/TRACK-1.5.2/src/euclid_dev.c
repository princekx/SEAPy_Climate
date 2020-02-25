#include <Stdio.h>
#include <Math.h>
#include "st_fo.h"
#include "grid.h"

/* function to compute the euclidian deviation for three true feature points */

extern int tf;
extern float w1, w2;
extern int x1u, y1u;

extern GRID *gr;

double euclid_dev(struct feature_pts *fp0, struct feature_pts *fp1, struct feature_pts *fp2)

{

    double diff1x, diff1y, diff2x, diff2y, mod1, mod2;
    double mm, phi;

    if(tf == 3){

      diff1x = *(gr->xgrid + ((fp0->x).ixy)+x1u-2) - *(gr->xgrid + ((fp1->x).ixy)+x1u-2);
      diff1y = *(gr->ygrid + ((fp0->y).ixy)+y1u-2) - *(gr->ygrid + ((fp1->y).ixy)+y1u-2);

      diff2x = *(gr->xgrid + ((fp1->x).ixy)+x1u-2) - *(gr->xgrid + ((fp2->x).ixy)+x1u-2);
      diff2y = *(gr->ygrid + ((fp1->y).ixy)+y1u-2) - *(gr->ygrid + ((fp2->y).ixy)+y1u-2);


     }

     else{

       diff1x = (fp0->x).xy - (fp1->x).xy;
       diff1y = (fp0->y).xy - (fp1->y).xy;

       diff2x = (fp1->x).xy - (fp2->x).xy;
       diff2y = (fp1->y).xy - (fp2->y).xy;

     }

     mod1 = sqrt(diff1x*diff1x + diff1y*diff1y);
     mod2 = sqrt(diff2x*diff2x + diff2y*diff2y);


     if(mod1 <= 0.0 && mod2 <= 0.0) phi = 0.0;

     else if(mod1 <= 0.0 || mod2 <= 0.0) phi = w2;

     else{

        mm = mod1*mod2;
        phi = w1*(1.0-((diff1x*diff2x+diff1y*diff2y)/mm))
              + w2*(1.0 - 2.0*(sqrt(mm)/(mod1+mod2)));

     }

     return phi;

}
