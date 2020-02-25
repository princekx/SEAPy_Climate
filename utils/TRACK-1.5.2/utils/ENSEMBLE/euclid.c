#include <stdio.h>
#include <stdlib.h>
#include "st_fo.h"

/* function to calculate the euclidian distance between points */


double euclid(struct feature_pts *f1, struct feature_pts *f2)

{

     double dd;
     double d1, d2;


     d1 = (f1->x).xy - (f2->x).xy;
     d2 = (f1->y).xy - (f2->y).xy;

     dd = (double)(d1*d1 + d2*d2);

     return dd;

}
