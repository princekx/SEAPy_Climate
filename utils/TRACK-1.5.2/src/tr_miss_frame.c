#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "files_in.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"
#include "mem_er.h"
#include "file_handle.h"
#include "m_values.h"

/* function to test track continuity for missing frames. */

double devn(struct frame_objs * , struct track_points * , struct track_points * , struct track_points * , float );
float disp(struct frame_objs *, struct frame_objs *, struct track_points *, struct track_points *, float );
int check_time(struct times *, struct times *, int );
struct track_ind *add_2track(struct track_ind * , int );


extern int track_num;
extern int tom;
extern int ffirst, fint, frterm;
extern int frcnt;

static int days[]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};


struct track_ind *tr_miss_frame(struct frame_objs *fo, struct track_ind *tind, float dmax, float phimax, int nfr)

{

  int nfram;
  int i,j,k;
  int imiss;
  int ibt;

  off_t chrnum, place1, place2;

  float tstep;
  float dm, phim;
  float dptr, devt;

  FILE *tim=NULL;

  char timfil[100];

  struct times *frtim=NULL, *frt, *frn;
  struct track_ind *tr, *trr;
  struct track_points *tp, *tps, *tpp1;

  printf("****WARNING****, this routine, %s, does not have the zonal upper bound \r\n"
         "                 displacment, adaptive tracking or already identified  \r\n"
         "                 missing frame support implemented yet.                \n\n", __FILE__);

  strcpy(timfil, FRTIMES);

  if(!fexist(timfil, "r")){

     printf("***WARNING***, file %s \r\n"
            "               does not exist for read, try a different file\n\n", timfil);

     printf("What is the name of the time file to be read?\n\n");
     scanf("%s", timfil);

     if(!fexist(timfil, "r")){
        printf("***ERROR***, file %s \r\n"
            "               does not exist for read, aborting\n\n", timfil);
        exit(1);
     }

  }

/* read in time data */

  tim = open_file(timfil, "r");

  fscanf(tim, "%d\n", &nfram);

  if(frcnt) frterm = frcnt;

  if(!frterm){
    printf("What where the frame details, START, INTERVAL, END\n\n");
    scanf("%d %d %d", &ffirst, &fint, &frterm);
  }

  if(frterm > nfram){
     printf("***ERROR***, not enough time data for frame range. \r\n"
            "             no missing frame checks performed.    \n\n");
     return tind;
  }

  frtim = (struct times *)calloc(nfr, sizeof(struct times));
  mem_er((frtim == NULL) ? 0 : 1, nfr * sizeof(struct times));

  place1 = ftello(tim);

  fscanf(tim, "%d %d %d %d\n", &frtim->year, &frtim->month, &frtim->day, &frtim->hour);
  place2 = ftello(tim);
  chrnum=place2-place1;
  fseeko(tim, place1, FSTART);

  fseeko(tim, (ffirst-1)*chrnum, ORIGIN);

  for(i=0; i< nfr; i++) {
     frt = frtim+i; 
     fscanf(tim, "%d %d %d %d\n", &frt->year, &frt->month, &frt->day, &frt->hour);
     if(fint > 1)fseeko(tim, (fint-1)*chrnum, ORIGIN);
  }

  close_file(tim, timfil);

  printf("What is the time step between frames?\n\n");
  scanf("%f", &tstep);


/* loop through tracks perform one forward and one backward loop */

retry:

  printf("The old TRACKING values are %f for displacement and %f for smoothness  \r\n"
         "what are the new values which must be more restrictive than old values.\r\n",
         dmax/FP_PI, phimax);

  scanf("%f %f", &dm, &phim);

  if(tom == 'g') dm *= FP_PI;

  if(dm > dmax || phim > phimax){

     printf("***ERROR***, new TRACKING values must be more restrictive, try again!\n");
     goto retry;

  }

  for(k=0; k< track_num; k++){

      tr = tind + k;

/* forward loop */

      if(tr->num_point_track > 0){

         tps = tr->tp;

         ibt = 0;

         while(tps->feature_id == -1) {++tps; ++ibt;}

         tp = tr->tp + ibt;

         for(i=1; i< tr->num_point_track-1; i++){

             tps = tp + i;

             dptr = disp(fo+tps->frame_id-1, fo+(tps+1)->frame_id-1, tps, tps+1, dm);

             devt = devn(fo, tps-1, tps, tps+1, phim); 


             frt = frtim + tps->frame_id - 1;
             frn = frtim + (tps+1)->frame_id - 1;

             imiss = check_time(frt, frn, tstep);

             tps->nmpt = imiss;

             if((dptr > dm || devt > phim) && !imiss){ 

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
                        tps->feature_id = -1;
                        tps->object_id = 0;
                        tpp1->nmpt = 0;

                        ++tpp1; ++tps;
                    }

                    break;

             }


         } 

      } 

/* backward loop */

      if(tr->num_point_track > 0){

         tps = tr->tp + nfr - 1;

         ibt = 0;

         while(tps->feature_id == -1) {--tps; ++ibt;}

         tp = tr->tp + nfr - ibt - 1;

         for(i=1; i < tr->num_point_track-1; i++){

             tps = tp - i;           

             dptr = disp(fo+tps->frame_id-1, fo+(tps-1)->frame_id-1, tps, tps-1, dm);
             devt = devn(fo, tps+1, tps, tps-1, phim);

             frt = frtim + (tps-1)->frame_id - 1;
             frn = frtim + tps->frame_id - 1;

             imiss = check_time(frt, frn, tstep);

             (tps-1)->nmpt = imiss;

             if((dptr > dm || devt > phim) && !imiss){  

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
                        tps->feature_id = -1;
                        tps->object_id = 0;
                        tpp1->nmpt = 0;
                        --tpp1; --tps;
                    }

                    break;

             }

         }

      } 

  }


  free(frtim);

  return tind;


}


int check_time(struct times *frt, struct times *frn, int tstep)

{

    int i;
    int nhr;
    int ydays, nmd, nyd, nstep;

    days[1] = 28;
    ydays = 365;

    if(!(frt->year % 4)) {days[1] = 29; ydays = 366;}

    if(frt->year != frn->year){

        nmd = 0;
        for(i=0; i< frt->month; i++)nmd += days[i];
        nyd = ydays - nmd + days[frt->month - 1] - frt->day;

        for(i=frt->year+1; i < frn->year; i++) nyd += (!(i % 4)) ? 366 : 365;
        days[1] = (!(frn->year % 4)) ? 29 : 28;
        for(i=0; i< frn->month -1; i++) nyd += days[i];
        nyd += frn->day - 1;


        nhr = nyd * 24 + (24 - frt->hour) + frn->hour;
        

    }


    else{

       if(frt->month != frn->month){

          nmd = days[frt->month-1] - frt->day + frn->day - 1;
          for(i=frt->month; i< frn->month-1; i++) nmd += days[i];
          nhr = nmd * 24 + (24 - frt->hour) + frn->hour;

       }

       else

          nhr = 24 * (frn->day - frt->day) - frt->hour + frn->hour;

    
     }


     nstep = nhr / tstep;

     return (nstep == 1) ? 0 : nstep - 1;



}


struct track_ind *add_2track(struct track_ind *tind, int nfr)
{
   int j;

   struct track_ind *trr;
   struct track_points *tpp1, *tpp2;

   track_num += 2;
   tind = (struct track_ind * )realloc_n(tind, track_num*sizeof(struct track_ind));
   mem_er((tind == NULL) ? 0 : 1, track_num*sizeof(struct track_ind));

   trr = tind + track_num - 2;

   trr->tp = (struct track_points * )calloc(nfr, sizeof(struct track_points));
   mem_er((trr->tp == NULL) ? 0 : 1, nfr * sizeof(struct track_points));

   (trr+1)->tp = (struct track_points * )calloc(nfr, sizeof(struct track_points));
   mem_er(((trr+1)->tp == NULL) ? 0 : 1, nfr * sizeof(struct track_points));

    trr->num_point_track = 0;
   (trr+1)->num_point_track = 0;

   tpp1 = trr->tp;
   tpp2 = (trr+1)->tp;
                    
   for(j=0; j < nfr; j++) {
      tpp1->frame_id = tpp2->frame_id = j+1;
      tpp1->object_id = tpp2->object_id = 0;
      tpp1->feature_id = tpp2->feature_id = -1;
      tpp1->nmpt = tpp2->nmpt = 0;
      ++tpp1; ++tpp2;
   }


   return tind;

}
