#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"
#include "grid.h"

/* function to produce plotting arrays for feature points to be plotted */

#ifdef NOUNDERSCORE

void featpl(float *, float *, int *, int *, float *, float *, float *, 
             int *, int *, int *, float *, float *, int *, int *, float *, 
             float *, float *, float *, int *, float *, float *, float *,
             float *, int *, int *, int *, int *, float *, float *, float *);

#else

void featpl_(float *, float *, int *, int *, float *, float *, float *, 
             int *, int *, int *, float *, float *, int *, int *, float *, 
             float *, float *, float *, int *, float *, float *, float *,
             float *, int *, int *, int *, int *, float *, float *, float *);

#endif

extern float w1, w2, dmax, phimax;
extern int tf, track_num, mnpt, t_id;
extern int x1u, y1u;
extern int izm;
extern float xmn, ymn, xmx, ymx;
extern int delb;

extern GRID *gr;
extern CNTRY *cm;

void featd(struct frame_objs *fo, struct track_ind *tind, int tot_fet, int frame_num, int k)

{

     int i, j, n, l, count = 0;
     int *tid = NULL;
     int arr=0;
     int ntr=0;
     int pl_vec=0;

     struct track_ind *tr;
     struct track_points *tps;
     struct feature_pts *fpts;
     struct frame_objs *ff;
     struct object *ob;

     float *xf=NULL, *yf=NULL, *zf=NULL;

     xf = (float *)calloc(tot_fet, sizeof(float));
     mem_er((xf == NULL) ? 0 : 1, tot_fet * sizeof(float));

     yf = (float *)calloc(tot_fet, sizeof(float));
     mem_er((yf == NULL) ? 0 : 1, tot_fet * sizeof(float));

     zf = (float *)calloc(tot_fet, sizeof(float));
     mem_er((zf == NULL) ? 0 : 1, tot_fet * sizeof(float));

     tid = (int *)calloc(tot_fet, sizeof(int));
     mem_er((tid == NULL) ? 0 : 1, tot_fet * sizeof(int));

     if(t_id == 0){
        n = 0;
        l = track_num;
     }
     else {
        n = t_id - 1;
        l = t_id;
     }

     for(i=n; i < l; i++){

         tr = tind + i;

         if(tr->num_point_track >= mnpt){

            for(j=0; j < frame_num; j++){

               tps = (tr->tp) + j;

               if(tps->feature_id != -1){

                  ff = fo + (tps->frame_id)-1;
                  ob = (ff->objs) + (tps->object_id)-1;
                  fpts = (ob->fet->fpt) + (tps->feature_id)-1;

                  zf[count] = fpts->str;

                  if(tf == 3){
                    xf[count] = *(gr->xgrid + ((fpts->x).ixy)+x1u-2);
                    yf[count] = *(gr->ygrid + ((fpts->y).ixy)+y1u-2);

                  }

                  else{

                    xf[count] = (fpts->x).xy;
                    yf[count] = (fpts->y).xy;

                  }

                  tid[count]=i+1;

                  ++count;

               }


            }

          }


     }


#ifdef   NOUNDERSCORE

     featpl(gr->xgrid, gr->ygrid, &gr->ix, &gr->iy, xf, yf, zf, tid, &tot_fet, 
             cm->cmi, cm->cmxg, cm->cmyg, &cm->dcm , &k, &w1, &w2, &dmax, &phimax,
             &izm, &xmn, &ymn, &xmx, &ymx, &delb, &arr, &ntr, &pl_vec, NULL,
             NULL, NULL);

#else

     featpl_(gr->xgrid, gr->ygrid, &gr->ix, &gr->iy, xf, yf, zf, tid, &tot_fet, 
             cm->cmi, cm->cmxg, cm->cmyg, &cm->dcm , &k, &w1, &w2, &dmax, &phimax,
             &izm, &xmn, &ymn, &xmx, &ymx, &delb, &arr, &ntr, &pl_vec, NULL,
             NULL, NULL);

#endif


     free(xf);
     free(yf);
     free(zf);
     free(tid);

     return;

}
