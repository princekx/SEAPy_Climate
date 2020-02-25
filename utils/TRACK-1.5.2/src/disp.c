#include <Stdio.h>
#include <stdlib.h>
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"

/* function to calculate the displacement between consecutive points on a track. */

float measure(struct feature_pts * , struct feature_pts * );

float disp(struct frame_objs *f1, struct frame_objs *f2, struct track_points *t1, struct track_points *t2, float dmax)

{


    float dd;

    struct feature_pts *fpts1, *fpts2;
    struct object *ob;

    if(t1->feature_id != -1 && t2->feature_id != -1){

        ob = (f1->objs) + (t1->object_id) - 1;
        fpts1 = (ob->fet->fpt) + (t1->feature_id) - 1;

        ob = (f2->objs) + (t2->object_id) - 1;
        fpts2 = (ob->fet->fpt) + (t2->feature_id) - 1;

        dd = measure(fpts1, fpts2);

    }

    else dd = dmax;

    return dd;

}
