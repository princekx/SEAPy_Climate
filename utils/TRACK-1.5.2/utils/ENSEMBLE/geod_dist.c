#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "st_fo.h"
#include "m_values.h"

/* function to calculate the geodisic distance between points */


double geod_dist(struct feature_pts *f1, struct feature_pts *f2)

{

     double dis;
     double thet1, phi1, thet2, phi2, dd;

     phi1 = (f1->x).xy * FP_PI; 
     phi2 = (f2->x).xy * FP_PI;
     thet1 = FP_PI2 - (f1->y).xy * FP_PI;
     thet2 = FP_PI2 - (f2->y).xy * FP_PI;

     if(thet1 < 0.) thet1 = 0.;
     if(thet2 < 0.) thet2 = 0.;

     dd = sin(thet1) * sin(thet2) * cos(phi2 - phi1) + cos(thet1) * cos(thet2);

     if(fabs(dd) > 1.) dd = (dd < 0.) ? -1.0 : 1.0;

     dis = acos(dd);

     return dis;

}
