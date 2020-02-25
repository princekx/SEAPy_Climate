#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include "st_fo.h"
#include "st_obj.h"
#include "grid.h"

/* function to calculate the euclidian distance between points */

extern int tf;
extern int x1u, y1u;

extern GRID *gr;

double euclid(struct feature_pts *f1, struct feature_pts *f2)

{

     double dd;
     double d1, d2;

     if(tf == 3){

        d1 = obj_xreal(((f1->x).ixy) + x1u - 2) - obj_xreal(((f2->x).ixy) + x1u - 2);
        d2 = *(gr->ygrid + ((f1->y).ixy) + y1u - 2) - *(gr->ygrid + ((f2->y).ixy) + y1u - 2);

     }

     else{

        d1 = (f1->x).xy - (f2->x).xy;
        d2 = (f1->y).xy - (f2->y).xy;

     }

     dd = sqrt((double)(d1*d1 + d2*d2));

     return dd;

}
