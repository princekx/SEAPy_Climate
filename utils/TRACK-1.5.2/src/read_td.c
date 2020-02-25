#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "mem_er.h"
#include "st_track.h"
#include "file_handle.h"


#define MAXCHAR   11
#define FIELDWD   40
#define MMXCHR   100

extern int track_num;
extern float w1, w2;
extern float dmax, phimax;

/* function to read in existing track data */

struct track_ind *read_td(FILE *tdump, int frame_num)

{

     int k, i, j, tr_st=0;
     int tr_id, fr_id, ob_id, fe_id, fe_nmpt=0;

     char tex1[MAXCHAR];
     char trdt[FIELDWD];
     char tmpch[MMXCHR];

     struct track_ind *tind=NULL;
     struct track_ind *tr=NULL;
     struct track_points *tps=NULL;

/* Check for old format */

     fgets(tmpch, MMXCHR, tdump);

     if(!strstr(tmpch, "TRACK_No.")){

        fseeko(tdump, (off_t)0, FSTART);
        fscanf(tdump, "%f %f", &w1, &w2);
        fscanf(tdump, "%f %f", &dmax, &phimax);

     }

     else {

        fseeko(tdump, (off_t)0, FSTART);
        printf("****WARNING****, old file format, you may have to specify\r\n"
               "                 parameters used in a previous run.      \n\n");

     }

     fscanf(tdump, "%s %d", tex1, &track_num);


     tind = (struct track_ind * )calloc(track_num, sizeof(struct track_ind));
     mem_er((tind == NULL) ? 0 : 1, track_num * sizeof(struct track_ind));

     for(k=0; k < track_num; k++){

          tr = tind + k;

          fscanf(tdump, "%s %d\n", tex1, &tr_id);
          if(!strstr(tex1, "TRACK_Id.")){

             printf("****ERROR****, reading track file \r\n"
                    "               no TRACK_Id. tag,  \r\n"
                    "               exiting.           \n\n");
             exit(1);

          }

          fscanf(tdump, "%s %d\n", tex1, &tr->num_point_track);


          tr->tp = (struct track_points * )calloc(frame_num, sizeof(struct track_points));
          mem_er((tr->tp == NULL) ? 0 : 1, frame_num * sizeof(struct track_points));

          for(i=0; i < tr->num_point_track; i++){

              fgets(trdt, FIELDWD, tdump);

              fe_nmpt = 0;
              sscanf(trdt, "%d %d %d %d", &fr_id, &ob_id, &fe_id, &fe_nmpt);


/*              fscanf(tdump, "%d %d %d", &fr_id, &ob_id, &fe_id); */

              if(i == 0){

                 tr_st = fr_id;

                 for(j=0; j < fr_id -1; j++){

                     tps = tr->tp + j;
                     tps->frame_id = j+1;
                     tps->object_id = 0;
                     tps->feature_id = -1;

                 }
               }

               tps = tr->tp + tr_st + i - 1;
               tps->frame_id = fr_id;
               tps->object_id = ob_id;
               tps->feature_id = fe_id;
               tps->nmpt = fe_nmpt;


          }

          if(tr->num_point_track < frame_num){

             for(i=tr_st+tr->num_point_track-1; i < frame_num; i++){

                 tps = tr->tp + i;
                 tps->frame_id = i+1;
                 tps->object_id = 0;
                 tps->feature_id = -1;

             }

           }

     }

     return tind;

}
