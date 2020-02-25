#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"


extern int track_num;

/* function to split tracks */

struct track_ind *track_split(struct frame_objs *fo, struct track_ind *tind, int fr)

{

     int k, it, trnm= track_num;
     int fst_pt, lst_pt, pt_str;

     struct track_ind *tr, *tr1, *trp;
     struct track_points *tps0, *tps1, *tpsp;
     struct feature_pts *fpts;
     struct object *ob;

     for(it=0; it < trnm; it++){

        tr = tind + it;

        lst_pt = 0;
        pt_str = 0;
        fst_pt = 0;

cont_srch:

        if(tr->num_point_track > 1){

            for(k=lst_pt; k < fr; k++){

                if(((tr->tp)+k)->feature_id > 0){

                   fst_pt = k;
                   break;

                }

             }

             for(k=fst_pt+1; k < fr; k++){

                 if(((tr->tp)+k)->feature_id == -1){

                    lst_pt = k;
                    break;

                 }
 
                 if(k == fr-1) lst_pt = k+1;

             }

             pt_str = lst_pt - fst_pt;

             if(pt_str < tr->num_point_track){

                tr->num_point_track -= pt_str;
                track_num += 2;

                tind = (struct track_ind * )realloc_n(tind, track_num*sizeof(struct track_ind));
                mem_er((tind == NULL) ? 0 : 1, track_num*sizeof(struct track_ind));

                tr1 = tind + track_num -2;
                trp = tind + track_num -1;

                tr1->tp = (struct track_points * )calloc(fr, sizeof(struct track_points));
                mem_er((tr1->tp == NULL) ? 0 : 1, fr * sizeof(struct track_points));

                trp->tp = (struct track_points * )calloc(fr, sizeof(struct track_points));
                mem_er((trp->tp == NULL) ? 0 : 1, fr * sizeof(struct track_points));

                tr = tind + it;

                tr1->num_point_track = pt_str;
                trp->num_point_track = 0;
                   
                for(k=0; k < fr; k++){

                    tps0 = (tr->tp)+k;
                    tps1 = (tr1->tp)+k;
                    tpsp = (trp->tp)+k;

                    if(k >= fst_pt && k < lst_pt){

                      ob = ((fo+k)->objs)+tps0->object_id -1;
                      fpts = (ob->fet->fpt) + tps0->feature_id -1;

                      tps1->frame_id = tps0->frame_id;
                      tps1->object_id = tps0->object_id;
                      tps0->object_id = 0;
                      tps1->feature_id = tps0->feature_id;
                      tps0->feature_id = -1; 
/*                      fpts->track_id = track_num -2; */

                    }

                    else {

                      tps1->frame_id = k+1;
                      tps1->object_id = 0;
                      tps1->feature_id = -1;

                    }
     
                    tpsp->frame_id = k+1;
                    tpsp->object_id = 0;
                    tpsp->feature_id = -1;

               }

               pt_str = 0;

               goto cont_srch;

             }
                    
        }


     }


/*      for(i=0; i < track_num; i++){

         printf("%d\n", i);

         tr = tind + i;
         printf("%d\n", tr->num_point_track);

         for(k=0; k < fr; k++){

             tps0 = (tr->tp)+k;

             printf("%d %d %d\n", tps0->frame_id, tps0->object_id, tps0->feature_id);

         }

         printf("\n");

      } */


     return tind;

}
