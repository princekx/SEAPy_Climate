#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "st_fo.h"
#include "st_obj.h"
#include "m_values.h"
#include "grid.h"

/* function to calculate the geodisic distance between points */

void geo_conv(struct feature_pts * );

extern int tf;
extern int x1u, y1u;

extern GRID *gr;
extern int geo_init;

double geod_dist(struct feature_pts *f1, struct feature_pts *f2)

{

     double dis;
     double thet1, phi1, thet2, phi2, dd;

     if(geo_init){

        if(!f1->gwky) geo_conv(f1);
        if(!f2->gwky) geo_conv(f2);

        dd = f1->gwk[0] * f2->gwk[0] + f1->gwk[1] * f2->gwk[1] + f1->gwk[2] * f2->gwk[2];

     }

     else {

        if(tf == 3){

           phi1 = obj_xreal(((f1->x).ixy) + x1u - 2) * FP_PI; 
           phi2 = obj_xreal(((f2->x).ixy) + x1u - 2) * FP_PI;
           thet1 = FP_PI2 - *(gr->ygrid + ((f1->y).ixy) + y1u - 2) * FP_PI; 
           thet2 = FP_PI2 - *(gr->ygrid + ((f2->y).ixy) + y1u - 2) * FP_PI;

        }

        else{

           phi1 = (f1->x).xy * FP_PI; 
           phi2 = (f2->x).xy * FP_PI;
           thet1 = FP_PI2 - (f1->y).xy * FP_PI;
           thet2 = FP_PI2 - (f2->y).xy * FP_PI;

        }

        if(thet1 < 0.) thet1 = 0.;
        if(thet2 < 0.) thet2 = 0.;

        dd = sin(thet1) * sin(thet2) * cos(phi2 - phi1) + cos(thet1) * cos(thet2);


     }
 
     if(fabs(dd) > 1.) dd = (dd < 0.) ? -1.0 : 1.0;

     dis = acos(dd);    

     return dis;

}
