#include <Stdio.h>
#include <stdlib.h>
#include "st_fo.h"
#include "st_obj.h"
#include "st_track.h"

/* cost function for the mge algorithm */

double euclid_dev(struct feature_pts * , struct feature_pts * , struct feature_pts * );
double geod_dev(struct feature_pts * , struct feature_pts * , struct feature_pts * );

int noph=0;

extern int tom;

double devn(struct frame_objs *fo, struct track_points *tp0 , struct track_points *tp1 , struct track_points *tp2 , float phimax)

{

    struct feature_pts *fp0, *fp1, *fp2;
    struct object *ob;

    noph = 0;

    if(tp0->feature_id == -1) return 0.0;

    else if(tp0->feature_id > 0 && 
            tp1->feature_id > 0 && 
            tp2->feature_id > 0)  {

               noph = 1;

               ob = ((fo + (tp0->frame_id) -1)->objs) + (tp0->object_id) - 1;
               fp0 = (ob->fet->fpt) + (tp0->feature_id) - 1;

               ob = ((fo + (tp1->frame_id) -1)->objs) + (tp1->object_id) - 1;
               fp1 = (ob->fet->fpt) + (tp1->feature_id) - 1;

               ob = ((fo + (tp2->frame_id) -1)->objs) + (tp2->object_id) - 1;
               fp2 = (ob->fet->fpt) + (tp2->feature_id) - 1;

               switch(tom){
                  case 'e':
                    return euclid_dev(fp0, fp1, fp2);
                  case 'g':
                    return geod_dev(fp0, fp1, fp2);
                  default:
                    printf("***error***, incorrect key used for type of measure in %s\n", __FILE__);
                    exit(1);
               }


    }

    else return phimax;

}
