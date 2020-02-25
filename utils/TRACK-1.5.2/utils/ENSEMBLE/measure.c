#include <stdio.h>
#include <stdlib.h>
#include "st_fo.h"

/* function to compute using diferant measures the displacement between 
   points on a track.                                                   */

double euclid(struct feature_pts * , struct feature_pts * );
double geod_dist(struct feature_pts * , struct feature_pts * );

extern int tom;

float measure(struct feature_pts *f1 , struct feature_pts *f2)

{

   switch(tom){
      case 'e':
        return (float) euclid(f1, f2);
      case 'g':
        return (float) geod_dist(f1, f2);
      default:
        printf("***error***, incorrect key used for type of measure in %s\n", __FILE__);
        exit(1);
    }

}

