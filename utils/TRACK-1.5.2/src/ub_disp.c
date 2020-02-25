#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"
#include "zones.h"

/* function to calculate the upper bound displacement for consecutive 
   points on a track.                                                 */

float dmaxx(ZONE * , struct feature_pts *f1, struct feature_pts *f2, float );


float ub_disp(struct frame_objs *f1, struct frame_objs *f2, struct track_points *t1, struct track_points *t2, ZONE *zone, float dmax)

{


    float dd;

    struct feature_pts *fpts1=NULL, *fpts2=NULL;
    struct object *ob;

    if(t1->feature_id != -1){

       ob = (f1->objs) + (t1->object_id) - 1;
       fpts1 = (ob->fet->fpt) + (t1->feature_id) - 1;

    }

    if(t2->feature_id != -1){

       ob = (f2->objs) + (t2->object_id) - 1;
       fpts2 = (ob->fet->fpt) + (t2->feature_id) - 1;

    }

    dd = dmaxx(zone, fpts1, fpts2, dmax);

    return dd;

}
