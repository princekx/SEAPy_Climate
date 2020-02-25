#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"
#include "zones.h"

#define  DINIT   1.0e+30

/* function to perform the track initialization for the modified greedy 
   exchange optimization algorithm.                                     */

float measure(struct feature_pts * , struct feature_pts * );
float dmaxx(ZONE * , struct feature_pts * , struct feature_pts * , float );

extern int track_num;
extern int izz;

struct track_ind *initialize_mge(ZONE **zone, struct frame_objs *fo, float *ddmax, int frame_num, int npar)

{
      int obn=0, fnum=0, ad;
      int i, j, k, l, m, n;
      int iim=0;

      float dist, distm;
      float ndmax=*ddmax;
      float dmax=*ddmax;

      struct frame_objs *foo=NULL;
      struct object *ob=NULL, *ob1=NULL;
      struct feature_pts *fpts=NULL, *fpts1=NULL, *fptm=NULL;
      struct track_ind *tind=NULL, *tr=NULL;
      struct track_points *tps=NULL;

      track_num = 0;

/* assign initial storage for tracks */

      track_num = 2 * (fo->tot_f_f_num);

      tind = (struct track_ind * )calloc(track_num, sizeof(struct track_ind));
      mem_er((tind == NULL) ? 0 : 1, track_num * sizeof(struct track_ind));


/* initialize for the first frame */

      ad = 0;

      for(i=0; i < fo->obj_num; i++){

          ob = (fo->objs)+i;

          for(j=0; j < ob->fet->feature_num; j++){

              tr = tind + ad +j;

              tr->tp = (struct track_points * )calloc(frame_num, sizeof(struct track_points));
              mem_er((tr->tp == NULL) ? 0 : 1, frame_num * sizeof(struct track_points));

              tr->num_point_track = 1;
              tr->tp->frame_id = 1;
              tr->tp->object_id = i+1;
              tr->tp->feature_id = j+1;
              tr->tp->nmpt = fo->nmiss;
              ((ob->fet->fpt)+j)->track_id = ad+j+1;

              tr = tr + 1;
              tr->num_point_track = 0;
              tr->tp = (struct track_points * )calloc(frame_num, sizeof(struct track_points));
              mem_er((tr->tp == NULL) ? 0 : 1, frame_num * sizeof(struct track_points));

              for(k=0; k < frame_num; k++){

                 tps = (tr->tp) + k;
                 tps->frame_id = k + 1;
                 tps->object_id = 0;
                 tps->feature_id = -1;
                 tps->nmpt = (fo + k)->nmiss;

              }

              ++ad;   

          }

          ad += ob->fet->feature_num;

      }

/* initialize tracks over the remaining frames */

      for(i=0; i < frame_num-1; i++){

          foo = fo+i;
          iim = (foo->nmiss > npar - 1) ? npar - 1 : foo->nmiss;
          dmax = ndmax = *(ddmax + iim);

          for(j=0; j < foo->obj_num; j++){

              ob = (foo->objs)+j;

              for(k=0; k < ob->fet->feature_num; k++){

                  fpts = (ob->fet->fpt)+k;

                  distm = DINIT;

                  for(l=0; l < (foo+1)->obj_num; l++){

                      ob1 = ((foo+1)->objs)+l;

                      for(m=0; m < ob1->fet->feature_num; m++){

                          fpts1 = (ob1->fet->fpt)+m;

                          dist = measure(fpts, fpts1);
                         
                          if(dist <= distm && fpts1->track_id == 0){

                            distm = dist;
                            obn = l+1;
                            fnum = m+1;
                            fptm = fpts1;
                            if(izz) ndmax = dmaxx(zone[iim], fpts, fpts1, dmax);

                          }

                      }

                  }

                  tr = tind + (fpts->track_id) - 1;

                  if(distm <= ndmax){

                     tps = (tr->tp)+i+1;
                     tps->frame_id = i+2;
                     tps->object_id = obn;
                     tps->feature_id = fnum;
                     tps->nmpt = (foo + 1)->nmiss;
                     fptm->track_id = fpts->track_id;
                     ++(tr->num_point_track);

                  }

                  else{

/* assign phantom points to remainder of current tracks */

                     for(m=i+1; m < frame_num; m++){

                         tps = (tr->tp) + m;
                         tps->frame_id = m+1;
                         tps->object_id = 0;
                         tps->feature_id = -1;
                         tps->nmpt = (fo + m)->nmiss;

                     }

                  }

              }

          }


          for(l=0; l < (foo+1)->obj_num; l++){

              ob1 = ((foo+1)->objs)+l;

              for(m=0; m < ob1->fet->feature_num; m++){

                 fpts1 = (ob1->fet->fpt)+m;

                 if(fpts1->track_id == 0){

                    track_num += 2;

                    tind = (struct track_ind * )realloc_n(tind,
track_num*sizeof(struct track_ind));
                    mem_er((tind == NULL) ? 0 : 1, track_num*sizeof(struct track_ind));

                    tr = tind + track_num - 2;

                    tr->tp = (struct track_points * )calloc(frame_num, sizeof(struct track_points));
                    mem_er((tr->tp == NULL) ? 0 : 1, frame_num * sizeof(struct track_points));

                    for(n=0; n < i+1; n++){

                        tps = (tr->tp) + n;
                        tps->frame_id = n+1;
                        tps->object_id = 0;
                        tps->feature_id = -1;
                        tps->nmpt = (fo + n)->nmiss;

                    }


                    tps = (tr->tp)+i+1;
                    tr->num_point_track = 1;
                    tps->frame_id = i+2;
                    tps->object_id = l+1;
                    tps->feature_id = m+1;
                    tps->nmpt = (foo + 1)->nmiss;
                    fpts1->track_id = track_num-1;

                    tr = tr + 1;
                    tr->num_point_track = 0;
                    tr->tp = (struct track_points * )calloc(frame_num, sizeof(struct track_points));
                    mem_er((tr->tp == NULL) ? 0 : 1, frame_num * sizeof(struct track_points));

                    for(k=0; k < frame_num; k++){

                       tps = (tr->tp) + k;
                       tps->frame_id = k + 1;
                       tps->object_id = 0;
                       tps->feature_id = -1;
                       tps->nmpt = (fo + k)->nmiss;
                    }

                 }

              }

          }

      }

      return tind;

}
