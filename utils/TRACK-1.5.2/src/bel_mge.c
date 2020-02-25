#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"
#include "zones.h"

#define TOLSM   1.0e-8    /* Tolenence on the max. gain search */

/* function to perform the mge algorithm in the forward direction recursively */

double devn(struct frame_objs * , struct track_points * , struct track_points * , struct track_points * , float );
float disp(struct frame_objs *, struct frame_objs *, struct track_points *, struct track_points *, float );
void track_fail(struct frame_objs *, struct track_ind *, int , int , int , int );
float ub_disp(struct frame_objs *, struct frame_objs *, struct track_points *, struct track_points * , ZONE * , float );

float phi(float , float , ADPT * );

extern int fi_fl, bi_fl;
extern int track_num;
extern int izz, iadpt, pproc;
extern int noph;

void bel_mge(ADPT **add, ZONE **zone, struct frame_objs *fo, struct track_ind *tind, float *dmax, float *phimax, int frame_num, int bi_count, int npar)

{

    int ex_fl;
    int i, j, k, iex=track_num, jex=track_num;
    int im1=0, im2=0;

    float dispij, dispji;
    float dis1, dis2;
    float dm=0.0, phim=0.0;

    float ndmax1=*dmax, ndmax2=*dmax;

    float nphi1, nphi2;

    double *trwsp;
    double gij, gijmax;
    double devii, devjj, devij, devji;
    double totgij;

    struct track_ind *tri, *trj;
    struct track_points *tpsi[3], *tpsj[3], tps;

    ex_fl = 0;

    trwsp = (double * )calloc(track_num, sizeof(double));
    mem_er((trwsp == NULL) ? 0 : 1, track_num * sizeof(double));

    for(k=frame_num-2; k > 0; k--){

        gijmax = 0.0;

        im1 = ((fo + k - 1)->nmiss > (fo + k)->nmiss) ? (fo + k - 1)->nmiss : (fo + k)->nmiss;
        if(im1 > npar - 1) im1 = npar - 1;
        phim = *(phimax + im1);

        for(i=0; i < track_num-1; i++){

            tri = tind + i;

            tpsi[0] = (tri->tp) + k-1;
            tpsi[1] = (tri->tp) + k;
            tpsi[2] = (tri->tp) + k+1;

            devii = (i == 0) ? (*(trwsp + i) = devn(fo, tpsi[2], tpsi[1], tpsi[0], phim)) : *(trwsp + i);

            for(j=i+1; j < track_num; j++){

                trj = tind + j;

                tpsj[0] = (trj->tp) + k-1;
                tpsj[1] = (trj->tp) + k;
                tpsj[2] = (trj->tp) + k+1;

                if(trj->num_point_track > 0 && i == 0)
                  *(trwsp + j) = devn(fo, tpsj[2], tpsj[1], tpsj[0], phim);

                if(pproc && (tri->num_point_track <= 0 || trj->num_point_track <= 0)) continue;

                if(tri->num_point_track <= 1 && trj->num_point_track <= 1) continue;

                if(tpsi[0]->feature_id != -1 || tpsj[0]->feature_id != -1){

/* check feature point pairs are within requisit distance */

                im2 = ((fo + k - 1)->nmiss > npar - 1) ? npar - 1 : (fo + k - 1)->nmiss;
                dm = *(dmax + im2);
                ndmax1 = ndmax2 = dm;

                if(izz){

                   ndmax1 = ub_disp(fo+k, fo+k-1, tpsi[1], tpsj[0], zone[im2], dm);
                   ndmax2 = ub_disp(fo+k-1, fo+k, tpsi[0], tpsj[1], zone[im2], dm);


                }

                dispij = disp(fo+k, fo+k-1, tpsi[1], tpsj[0], ndmax1);

                dispji = disp(fo+k, fo+k-1, tpsj[1], tpsi[0], ndmax2);

/* if feature points are within the requisite distance compute the exchange gain */

                   if(dispij <= ndmax1 && dispji <= ndmax2){

                      devjj = *(trwsp + j);

                      if(iadpt){

                        im2 = ((fo + k)->nmiss > npar - 1) ? npar - 1 : (fo + k)->nmiss;
                        dm = *(dmax + im2);

                        nphi1 = nphi2 = phim;

                        devij = devn(fo, tpsi[2], tpsi[1], tpsj[0], phim);

                        if(noph)

                           nphi1 = phi(dispij, disp(fo+k+1, fo+k, tpsi[2], tpsi[1], dm), add[im1]);

                        devji = devn(fo, tpsj[2], tpsj[1], tpsi[0], phim);

                        if(noph)

                           nphi2 = phi(dispji, disp(fo+k+1, fo+k, tpsj[2], tpsj[1], dm), add[im1]);

                         gij = devii + devjj - devij - devji;

                         if(gij - gijmax > TOLSM && 
                            (devij <= nphi1 && devji <= nphi2)){

                            gijmax = gij;
                            iex = i;
                            jex = j;


                         }


                      }

                      else {

                         devij = devn(fo, tpsi[2], tpsi[1], tpsj[0], phim);
                         devji = devn(fo, tpsj[2], tpsj[1], tpsi[0], phim);

                         gij = devii + devjj - devij - devji;

                         if(gij - gijmax > TOLSM){

                            gijmax = gij;
                            iex = i;
                            jex = j;
  
                         }

                      }

                   }

                }

            }

        }

        for(i=0; i < track_num; i++) *(trwsp+i) = 0.0;

/* exchange points if required */

        if(gijmax > 0.0){

           tri = tind + iex;
           trj = tind + jex;

           tpsi[0] = (tri->tp) + k-1;
           tpsj[0] = (trj->tp) + k-1;

           if(tpsi[0]->feature_id == -1 && tpsj[0]->feature_id != -1){

              ++(tri->num_point_track);
              --(trj->num_point_track);

           }

           else if(tpsi[0]->feature_id != -1 && tpsj[0]->feature_id == -1){

              --(tri->num_point_track);
              ++(trj->num_point_track);

           }

           tps = *tpsj[0];
           *tpsj[0] = *tpsi[0];
           *tpsi[0] = tps;

           if(k-1 > 0){

             im2 = ((fo + k - 2)->nmiss > npar - 1) ? npar - 1 : (fo + k - 2)->nmiss;
             dm = *(dmax + im2);
             ndmax1 = ndmax2 = dm;

             if(izz){

                 ndmax1 = ub_disp(fo+k-1, fo+k-2, tpsi[0], tpsi[0]-1, zone[im2], dm);
                 ndmax2 = ub_disp(fo+k-1, fo+k-2, tpsj[0], tpsj[0]-1, zone[im2], dm);

              }

              dis1 = disp(fo+k-1, fo+k-2, tpsi[0], tpsi[0]-1, ndmax1);
              dis2 = disp(fo+k-1, fo+k-2, tpsj[0], tpsj[0]-1, ndmax2);

              if(dis1 > ndmax1) track_fail(fo, tind, iex, k, frame_num, 'b');
              if(dis2 > ndmax2) track_fail(fo, tind, jex, k, frame_num, 'b');

           }
 
           if(ex_fl == 0) ex_fl = 1;

        }

    }

    free(trwsp);

    totgij = 0.0;

    for(i=0; i < track_num; i++){

       tri = tind + i;

       if(tri->num_point_track != 0){

          for(k=1; k < frame_num-1; k++){

               tpsi[0] = (tri->tp) + k-1;
               tpsi[1] = (tri->tp) + k;
               tpsi[2] = (tri->tp) + k+1;

               im1 = ((fo + k - 1)->nmiss > (fo + k)->nmiss) ? (fo + k - 1)->nmiss : (fo + k)->nmiss;
               if(im1 > npar - 1) im1 = npar - 1;

               totgij += devn(fo, tpsi[0], tpsi[1], tpsi[2], *(phimax + im1));              
          }

       }

    }

    ++bi_count;

    printf("backward iteration %d completed, cost function= %f\n", bi_count, totgij);

    if(ex_fl == 1){

       fi_fl = 1;
       bel_mge(add, zone, fo, tind, dmax, phimax, frame_num, bi_count, npar);

    }

    else bi_fl = 0;

    
    return;

}
