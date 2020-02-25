#include <Stdio.h>
#include "st_track.h"

/* function to write track data to file */

extern float w1, w2;
extern float dmax, phimax;

void trackd(struct track_ind *tind, int track_num, int frame_num, FILE *tdump, int re_tr_num)

{

    int i, j;
    int numc=0;

    struct track_ind *tr;
    struct track_points *tps;

    fprintf(tdump, "%f %f\n", w1, w2);
    fprintf(tdump, "%f %f\n", dmax, phimax);


    fprintf(tdump, "TRACK_No. %d\n", re_tr_num);

    for(i=0; i < track_num; i++){

        tr = tind + i;

        numc = 0;


        if(tr->num_point_track > 0){

           fprintf(tdump, "TRACK_Id. %d\n", i+1);
           fprintf(tdump, "POINT_No. %d\n", tr->num_point_track);

           for(j=0; j < frame_num; j++){

               tps = (tr->tp) + j;

               if(tps->feature_id != -1) {

                 ++numc;

                 if(tps->nmpt)
                    fprintf(tdump, "%d %d %d %d \n", tps->frame_id, tps->object_id, tps->feature_id, tps->nmpt);
                 else
                    fprintf(tdump, "%d %d %d \n", tps->frame_id, tps->object_id, tps->feature_id);


               }

           }

        }

        if(numc != tr->num_point_track){

           printf("****ERROR****, writing track file, the number of points \r\n"
                  "               on the track %d is not consistent with   \r\n"
                  "               the tagged point number of %d.           \r\n"
                  "               TRACK_ID = %d,                           \r\n"
                  "               continuing writing file.                 \n\n",
                  i+1, numc, tr->num_point_track); 

        }            

     }

     return;

}
