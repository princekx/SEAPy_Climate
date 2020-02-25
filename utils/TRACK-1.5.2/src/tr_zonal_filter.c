#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "files_in.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"
#include "mem_er.h"
#include "file_handle.h"
#include "m_values.h"
#include "zones.h"
#include "files_out.h"

/* function to test adaptive constraints */

double devn(struct frame_objs * , struct track_points * , struct track_points * , struct track_points * , float );
float disp(struct frame_objs *, struct frame_objs *, struct track_points *, struct track_points * ,float );
float ub_disp(struct frame_objs *, struct frame_objs *, struct track_points *, struct track_points * , ZONE * , float );
struct track_ind *add_2track(struct track_ind * , int );
ZONE *read_zones(int );
float phi(float , float , struct adapt * );
ADPT *read_adptp(int );

extern int track_num;
extern int ffirst, fint, frterm;
extern int frcnt;
extern int irmiss;
extern int tom;

struct track_ind *tr_zonal_filter(ADPT **add, ZONE **zone, struct frame_objs *fo, struct track_ind *tind, float **ddmax, float **pphimax, int nfr, int zz, int appt, int fbb, int npar)
{

  int i,j,k;
  int ibt;
  int adtr = 'n';
  int ubreg = 'n';
  int imsg=0;
  int tnump=0;
  int im=0, im1=0, im2=0;
  int nnpar=0;

  float *dmax, *phimax;
  float dm, phim;
  float dptr, dptr2, devt;
  float ndm;


  struct track_ind *tr, *trr;
  struct track_points *tp, *tps, *tpp1;
  struct frame_objs *f0, *f1, *f2;

  char dmphi_file[MAXCHR];

  FILE *fdmphi=NULL;

  dmax = *ddmax;
  phimax = *pphimax;
  ndm = dm = *dmax;
  phim = *phimax;

  if(irmiss && !appt){

     printf("There maybe missing frames do you want the tracking to be aware of them, 'y' or 'n'. \r\n"
            "Note this assumes the missing frames have been detected during the identification.   \n\n");
     scanf("\n");
     if(getchar() == 'y') {
        printf("How many sets of frame tracking parameters are required?     \r\n"
               "If more frames are missing than parameter sets are used,     \r\n"
               "then the last set are used.                                  \n\n");
        scanf("%d", &nnpar);
        if(nnpar < 1) {
           printf("****WARNING****, number of sets of parameters must be at least 1, resetting to 1.\n\n");
           nnpar = 1;
        }

        if(nnpar > npar){

           *ddmax = (float *)realloc_n(*ddmax, npar * sizeof(float));
           mem_er((*ddmax == NULL) ? 0 : 1, npar*sizeof(float));
           *pphimax = (float *)realloc_n(*pphimax, npar * sizeof(float));
           mem_er((*pphimax == NULL) ? 0 : 1, npar*sizeof(float));

           printf("What file contains the dmax,phimax %d pairs.\n\n", npar);
           scanf("%s", dmphi_file); 

           dmax = *ddmax;
           phimax = *pphimax;

           fdmphi = open_file(dmphi_file, "r");
           for(i=0; i < npar; i++) {
               fscanf(fdmphi, "%f %f", dmax + i, phimax + i);
               if(tom == 'g') *(dmax + i) *= FP_PI;
           }
           close_file(fdmphi, dmphi_file);

           ndm = dm = *dmax;
           phim = *phimax;

        }

     }


  }


  if(!zz){

     printf("Do you want regional upper bound displacements, 'y' or 'n'\n\n");
     scanf("\n");
     ubreg = getchar();

  }

  else if (zz < 0) ubreg = 'n';

  else ubreg = 'y';


  if(ubreg == 'y' && !zone) {

     zone = (ZONE **)calloc(npar, sizeof(ZONE *));
     mem_er((zone == NULL) ? 0 : 1, npar * sizeof(ZONE *));

     for(i=0; i < npar; i++) zone[i] = read_zones(i);

  }

  if(!appt){

     printf("Do you want adaptive tracking? 'y' or 'n'\n\n");
     scanf("\n");
     adtr=getchar();

  }

  else if(appt < 0) adtr = 'n';

  else adtr = 'y';

  if(adtr == 'y' && !add){

     add = (ADPT **)calloc(npar, sizeof(ADPT *));
     mem_er((add == NULL) ? 0 : 1, npar * sizeof(ADPT *));

     for(i=0; i < npar; i++) add[i] = read_adptp(i);

   }


/* loop through tracks perform one forward and one backward loop */


  for(k=0; k< track_num; k++){

      tr = tind + k;

/*      if(!(tr->num_point_track)) continue; */

/* perform track check for consistency */

     tnump = 0;

     for(i=0; i < nfr; i++)

         if((tr->tp + i)->feature_id >0) ++tnump;


     if(tnump != tr->num_point_track){


        printf("****WARNING****, inconsistent number of points on Track %d\r\n"
               "                 Resetting point number. ***FIRST***     \n\n", k+1);

        tr->num_point_track = tnump;

      }

/* forward loop */

      if(tr->num_point_track > 0 && (fbb == 'f' || fbb == 'a')){

         if(!imsg)

            printf("***INFORMATION***, forward pass with track filter/checking.\n\n");

         tps = tr->tp;

         ibt = 0;

         while(tps->feature_id == -1) {++tps; ++ibt;}

         tp = tr->tp + ibt;

         for(i=1; i< tr->num_point_track-1; i++){

             tps = tp + i;

             f1 = fo + tps->frame_id - 1;
             f2 = fo + (tps+1)->frame_id - 1;

             im1 = (f1->nmiss > npar - 1) ? npar - 1 : f1->nmiss;
             ndm = dm = *(dmax + im);

             dptr = disp(f1, f2, tps, tps+1, dm);

             f0 = fo + (tps-1)->frame_id - 1;

             im2 = (f0->nmiss > f1->nmiss) ? f0->nmiss : f1->nmiss;
             if(im2 > npar - 1) im2 = npar - 1;

             if(adtr == 'y'){

                im = (f0->nmiss > npar - 1) ? npar - 1 : f0->nmiss;
                dm = *(dmax + im);

                dptr2 = disp(f0, f1, tps-1, tps, dm);
                
                phim = phi(dptr, dptr2, add[im2]);


             }

             else phim = phimax[im2];

             if(ubreg == 'y')

                ndm = ub_disp(f1, f2, tps, tps+1, zone[im1], ndm);


             devt = devn(fo, tps-1, tps, tps+1, phim); 


             if(dptr > ndm || devt > phim){ 

                    tind = add_2track(tind, nfr);

                    tr = tind + k;
                    tps = tr->tp + ibt + i;

                    trr = tind + track_num - 2;

                    trr->num_point_track = tr->num_point_track - i - 1;
                    tr->num_point_track = i + 1; 

                    tpp1 = (trr->tp) + ibt + i + 1;
                    ++tps;

                    for(j=0; j < trr->num_point_track; j++){
                        tpp1->frame_id = tps->frame_id;
                        tpp1->object_id = tps->object_id;
                        tpp1->feature_id = tps->feature_id;
                        tpp1->nmpt = tps->nmpt;
                        tps->feature_id = -1;
                        tps->object_id = 0;

                        ++tpp1; ++tps;
                    }

                    break;

             }


         } 

      } 

/* backward loop */

      if(tr->num_point_track > 0 && (fbb == 'b' || fbb == 'a')){

         if(!imsg)

            printf("***INFORMATION***, backward pass with track filter/checking.\n\n");

         tps = tr->tp + nfr - 1;

         ibt = 0;

         while(tps->feature_id == -1) {--tps; ++ibt;}

         tp = tr->tp + nfr - ibt - 1;

         for(i=1; i < tr->num_point_track-1; i++){

             tps = tp - i;

             f0 = fo + (tps-1)->frame_id - 1;
             f1 = fo + tps->frame_id - 1;

             im1 = (f0->nmiss > npar - 1) ? npar - 1 : f0->nmiss;
             ndm = dm = *(dmax + im);       

             dptr = disp(f1, f0, tps, tps-1, dm);

             f2 = fo + (tps+1)->frame_id - 1;

             im2 = (f0->nmiss > f1->nmiss) ? f0->nmiss : f1->nmiss;
             if(im2 > npar - 1) im2 = npar - 1;

             if(adtr == 'y'){

                im = (f1->nmiss > npar - 1) ? npar - 1 : f1->nmiss;
                dm = *(dmax + im);

                dptr2 = disp(f2, f1, tps+1, tps, dm);
                phim = phi(dptr, dptr2, add[im2]);

             }

             else phim = phimax[im2];


             if(ubreg == 'y')

                ndm = ub_disp(f1, f0, tps, tps-1, zone[im1], ndm);

             devt = devn(fo, tps+1, tps, tps-1, phim);


             if(dptr > ndm || devt > phim){  

                    tind = add_2track(tind, nfr);

                    tr = tind + k;

                    tps = tr->tp + nfr - ibt - 1 - i;

                    trr = tind + track_num - 2;


                    trr->num_point_track = tr->num_point_track - i - 1;
                    tr->num_point_track = i + 1; 

                    tpp1 = (trr->tp) + nfr - 1 - ibt - i - 1;

                    --tps;

                    for(j=0; j < trr->num_point_track; j++){
                        tpp1->frame_id = tps->frame_id;
                        tpp1->object_id = tps->object_id;
                        tpp1->feature_id = tps->feature_id;
                        tpp1->nmpt = tps->nmpt;
                        tps->feature_id = -1;
                        tps->object_id = 0;
                        --tpp1; --tps;
                    }


                    break;

             }

         }

      } 

      if(!imsg) imsg = 1;


/* perform track check for consistency */

     tnump = 0;

     for(i=0; i < nfr; i++)

         if((tr->tp + i)->feature_id > 0) ++tnump;


     if(tnump != tr->num_point_track){


        printf("****WARNING****, inconsistent number of points on Track %d\r\n"
               "                 resetting point number.  ***LAST***     \n\n", k+1);

        tr->num_point_track = tnump;

      }

  }

  if(ubreg == 'y' && !zz){

    for(i=0; i < npar; i++) {
        free(zone[i]->zlat);
        free(zone[i]);
    }
    free(zone);
  }

  if(!appt && add) {
    for(i=0; i < npar; i++) free(add[i]);
    free(add);
  }

  return tind;


}


