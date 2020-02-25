#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include <string.h>
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"
#include "files_out.h"
#include "m_values.h"
#include "grid.h"
#include "mem_er.h"
#include "zones.h"
#include "file_handle.h"

#define  iter_term  40
#define  tot_term  3

/* function to perform the tracking optimization using the modified
   greedy exchange algorithm of Sethi and Jain with the modifications
   of Salari and Sethi included to cope with occlusion.               */

void bel_mge(ADPT ** , ZONE ** , struct frame_objs * , struct track_ind * , float * , float * , int , int , int );

double devn(struct frame_objs * , struct track_points * , struct track_points * , struct track_points * , float );

struct track_ind *initialize_mge(ZONE ** , struct frame_objs * , float * , int , int );

void featd(struct frame_objs * , struct track_ind * , int , int , int );

struct frame_objs *feature_pt_filter(ZONE ** , struct frame_objs * , float * , int , int );

void fel_mge(ADPT ** , ZONE ** , struct frame_objs * , struct track_ind * , float * , float * , int , int , int );

void objectd(struct frame_objs *, FILE * ,int , int );

struct track_ind *read_td(FILE * , int );

void trackd(struct track_ind * , int , int , FILE * , int );

struct track_ind *track_split(struct frame_objs * , struct track_ind * , int );

struct track_ind *tr_miss_frame(struct frame_objs * , struct track_ind * , float , float , int );

struct track_ind *tr_zonal_filter(ADPT ** , ZONE ** , struct frame_objs * , struct track_ind * , float ** , float ** , int , int , int , int , int );

ADPT *read_adptp(int );

ZONE *read_zones(int );

int track_num, mnpt, fp_count, t_id;
int fi_fl, bi_fl, izm=0;
float w1=0., w2=0.;
float dmax=0., phimax=0.;
int geo_init=0;
int izz=0, iadpt=0, pproc=0;

/*extern int track_num, mnpt, t_id, tf;*/
extern int tf;
extern int fi_fl, bi_fl;
extern int tom;
extern GRID *gr, *gr1;
extern int aniso;
extern int i_shape;
extern int irmiss;

extern char *fext;

extern int iext;

extern int nf, nfld;
extern int *nfwpos;

struct track_ind *mge_tracks(struct frame_objs *fo, int frame_num)

{

    int i=0, j, k, tot_fet, tot_count;
    int fi_count, bi_count, re_tr_num;
    int ch, rin='n', rc;
    int rd;
    int npar=1;
    int iim=0;

    char charin[MAXCHR], *filpt=NULL;
    char fout1[MAXCHR],fout2[MAXCHR],fout3[MAXCHR];
    char dmphi_file[MAXCHR];

    double totgij;
    float *odmax=NULL;
    float *ddmax=NULL, *pphimax=NULL;

    FILE *tdump=NULL, *idump=NULL, *fobjo=NULL, *fdmphi=NULL;

    struct track_ind *tind=NULL, *tr=NULL;
    struct frame_objs *foo=NULL;
    struct object *ob=NULL;
    struct feature_pts *fpts=NULL;
    struct track_points *tpsi[3]={NULL, NULL, NULL}, *tps=NULL;

    ZONE **zone=NULL;
    ADPT **add=NULL;

    if(i_shape){

       printf("****WARNING****, object data has been re-read, this means that  \r\n"
              "                 thresholds and other settings may have changed.\r\n"
              "                 Do you want to re-do the tracking, 'y' or 'n'. \n\n");
       scanf("\n");
       if(getchar() == 'n') return NULL;
    }

    strncpy(fout1, NEWOBJF, MAXCHR);
    strncpy(fout2, IDUMP, MAXCHR);
    strncpy(fout3, TDUMP, MAXCHR);

    if(iext){
        strcpy(strstr(fout1, EXTENSION), fext);
        strcpy(strstr(fout2, EXTENSION), fext);
        strcpy(strstr(fout3, EXTENSION), fext);
    }

    mnpt = 0;

    t_id = 0;

/* initialize switch if using geodesic measures */

    if(tom == 'g'){

      geo_init = 1;

      for(i=0; i < frame_num; i++){

         foo = fo + i;

         for(j=0; j < foo->obj_num; j++){

             ob = (foo->objs) + j;

             for(k=0; k < ob->fet->feature_num; k++){

                fpts = (ob->fet->fpt) + k;
                fpts->gwky = 0;

             }

          }

      }



    }

/* is an existing data set of tracks required to be plotted */

    printf("do you want to load an existing file of track data, y or n\n");

    scanf("\n");
    rd = getchar();

    if(rd == 'y'){

       printf("do you want to use a differant file to the default 'y' or 'n'\n");
       scanf("\n");
       if(getchar() == 'y'){

          printf("***WARNING***, make sure track file corresponds to the object file already read\n\n");
          printf("input file name\n");
          scanf("%s", charin);
          filpt = charin;

       }

       else filpt = fout3;

       printf("****INFORMATION****, reading track data from file, \r\n"
              "%s\n\n", filpt);

       tdump = open_file(filpt, "r");

       tind = read_td(tdump, frame_num);

       close_file(tdump, fout3);

       printf("Do you want to apply further analysis to the tracks,            \r\n"
              "e.g. missing frame search, zonal displacement bounds, 'y' or 'n'\n\n");

       scanf("\n");
       if(getchar() == 'y'){


          printf("****WARNING****, zonal displacement bounds are only applied on the sphere\n\n");

          if(w1 <= 0. || w2 <= 0.){

             printf("What where the original weights w1, w2?\n\n");

             scanf("%f %f", &w1, &w2);

             printf("What where the original values of the phantom point parameters dmax and phimax?\n\n");

             scanf("%f %f", &dmax, &phimax);

             if(tom == 'g') dmax *= FP_PI;

          }

          npar = 1;

          ddmax = (float *)calloc(npar, sizeof(float));
          mem_er((ddmax == NULL) ? 0 : 1, npar * sizeof(float));

          odmax = (float *)calloc(npar, sizeof(float));
          mem_er((odmax == NULL) ? 0 : 1, npar * sizeof(float));
          for(i=0; i < npar; i++) *(odmax + i) = 0.0;

          pphimax = (float *)calloc(npar, sizeof(float));
          mem_er((pphimax == NULL) ? 0 : 1, npar * sizeof(float));

          *ddmax = dmax;
          *pphimax = phimax;

          goto miss_zone;

       }

       goto replot;

    }

new_param:

    printf("what input parameters are required for the modified\r\n"
           "greedy exchange algorithm\r\n\n"
           "input the weights such that w1+w2=1 and w1, w2 < 1\n");

    scanf("%f %f", &w1, &w2);

    if(fabs(w1 + w2 - 1.0) > 1.0e-4){

       printf("****WARNING****, weights do not sum to unity.\n\n");

    }

    w1 *= 0.5;   /* correction to ensure that local smoothness lies in (0, 1) */

re_initial:

    printf("the current norm is %s \n\n", (tom == 'g') ? "GEODESIC" : " EUCLIDEAN");

irmiss = 1;

    if(irmiss){

       printf("There maybe missing frames do you want the tracking to be aware of them, 'y' or 'n'. \r\n"
              "Note this assumes the missing frames have been detected during the identification.   \n\n");
       scanf("\n");
       if(getchar() == 'y') {
          printf("How many sets of frame tracking parameters are required?     \r\n"
                 "If more frames are missing than parameter sets are used,     \r\n"
                 "then the last set are used.                                  \n\n");
          scanf("%d", &npar);
          if(npar < 1) {
             printf("****WARNING****, number of sets of parameters must be at least 1, resetting to 1.\n\n");
             npar = 1;
          }

       }

    }

    ddmax = (float *)calloc(npar, sizeof(float));
    mem_er((ddmax == NULL) ? 0 : 1, npar * sizeof(float));

    odmax = (float *)calloc(npar, sizeof(float));
    mem_er((odmax == NULL) ? 0 : 1, npar * sizeof(float));
    for(i=0; i < npar; i++) *(odmax + i) = 0.0;

    pphimax = (float *)calloc(npar, sizeof(float));
    mem_er((pphimax == NULL) ? 0 : 1, npar * sizeof(float));


    izz = 0;

    printf("Do you want regional upper bound displacements, 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y') {


       zone = (ZONE **)calloc(npar, sizeof(ZONE *));
       mem_er((zone == NULL) ? 0 : 1, npar * sizeof(ZONE *));

       for(i=0; i < npar; i++) zone[i] = read_zones(i);

       izz = 1;

    }

    iadpt = 0;

    printf("Do you want adaptive tracking, 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y'){

       add = (ADPT **)calloc(npar, sizeof(ADPT *));
       mem_er((add == NULL) ? 0 : 1, npar * sizeof(ADPT *));

       for(i=0; i < npar; i++) add[i] = read_adptp(i);

       iadpt = 1;
       pproc = 1;

    }

    else {

        printf("Do you want to apply pre-proccesing of tracks before\r\n"
               "each optimization loop, this can speed things up,   \r\n"
                "'y' or 'n'.                                        \n\n");
        scanf("\n");

        if(getchar() == 'y') pproc = -1;


    }

    if(npar == 1){
    
       printf("input the phantom point parameters dmax and phimax, \r\n\n"
              "----------------------------------------------------\r\n"
              "dmax is the max. upperbound displacement,           \r\n\n"
              "dmax is the euclidean norm for non-spherical space, \r\n\n"
              "dmax is the angular seperation (in degrees)         \r\n"
              "     on a great circle for the geodeisic norm       \r\n\n"
              "phimax is the penalty value of the cost function for\r\n"
              "a phantom point.                                    \r\n"
              "Value of dmax should be >= ZONEMAX if variable      \r\n"
              "upper bound dispalcements are being used.           \r\n"
              "----------------------------------------------------\n"); 

       scanf("%f %f", &dmax, &phimax);
       if(tom == 'g') dmax *= FP_PI;
       *ddmax = dmax;
       *pphimax = phimax;

    }

    else {

       printf("What file contains the dmax,phimax %d pairs.\n\n", npar);
       scanf("%s", dmphi_file); 

       fdmphi = open_file(dmphi_file, "r");
       for(i=0; i < npar; i++) {
           fscanf(fdmphi, "%f %f", ddmax + i, pphimax + i);
           if(tom == 'g') *(ddmax + i) *= FP_PI;
       }
       close_file(fdmphi, dmphi_file);

       dmax = *ddmax;
       phimax = *pphimax;

    }

    if(iadpt){

       for(i=0; i < npar; i++){

          if(*(pphimax + i) < add[i]->maxad){

             printf("****WARNING****, you may find the value of phimax chosen    \r\n"
                    "                 incompatable with the preprocessing option.\r\n"
                    "                 Resetting phimax to maxadd.                \n\n");
                    *(pphimax + i) = add[i]->maxad;
                    phimax = *pphimax;
          }

       }

     }

     if(izz){

       for(i=0; i < npar; i++){
 
          if(fabs(*(ddmax + i) - zone[i]->zonemax) > 1.0e-4) {

             printf("****WARNING****, chosen value of dmax != ZONEMAX, resetting dmax to zonemax\n\n");
             *(ddmax + i) =  zone[i]->zonemax;
             dmax = *ddmax;

          }

       }

    }

    if(!tind){

       if(rin == 'n'){

          fo = feature_pt_filter(zone, fo, ddmax, frame_num, npar);
       
          fobjo = open_file(file_exist(fout1, "r", APPS), "w");

          if(aniso == 'y') fprintf(fobjo, "%d %d %d\n", tf, 1, irmiss);
          else fprintf(fobjo, "%d %d %d\n", tf, 0, irmiss);
 
          fprintf(fobjo, "PROJ_DETAILS\n");
          fprintf(fobjo, "%d %d\n", gr->prgr, gr->prty);
          if(gr->prgr) fprintf(fobjo, "%f %f\n", gr->alat, gr->alng);
          
          fprintf(fobjo, "ADDITIONAL_FIELDS\n");
          fprintf(fobjo, "%3d %3d &", nf, nfld);
          if(nfwpos){
            for(i=0; i < nf; i++)fprintf(fobjo, "%1d", *(nfwpos + i));
          }
          fprintf(fobjo, "\n");

          for(i=0; i < frame_num; i++) objectd(fo+i, fobjo, i+1, tf);

          close_file(fobjo, fout1);

       }

       tind = initialize_mge(zone, fo, ddmax, frame_num, npar); 


/* how many real tracks are there */

       re_tr_num = 0;
       for(i=0; i < track_num; i++){
 
           tr = tind + i;
           if(tr->num_point_track > 0) ++re_tr_num;

       }


/* write initial tracks to file */

       idump = open_file(file_exist(fout2, "r", APPS), "w");

       trackd(tind, track_num, frame_num, idump, re_tr_num);

       close_file(idump, fout2);

       tot_fet = 0;

       for(i=0; i < track_num; i++){

           tr = tind + i;
           if(tr->num_point_track >= mnpt) tot_fet += tr->num_point_track;

       }

/* plot the initial track data if required by user */

      printf("do you want the initial track data plotted \r\n"
             "type y(not Y) or n(not N)!\n");

      scanf("\n");
      ch = getchar(); 

      if(ch == 'y') {

         printf("do you want area of interest (input 1) or whole area (input 0)\n");
         scanf("%d", &izm);
         featd(fo, tind, tot_fet, frame_num, 1); 

      }

      if(NOTRACK) exit(0);

   }

   else {

      for(i=0; i < npar; i++){

         if(*(ddmax + i) < *(odmax + i)){
            *(ddmax + i) = *(odmax + i);
            printf("***WARNING***, can't use a smaller dmax, resetting to old value\r\n"
                   "               if you want a more restrictive dmax the tracks  \r\n"
                   "               must be re-initialized.\n\n");

         }

      }

   }

   printf("do you want a different initialization or abort before minimization, 'y' or 'n' or 'a' for abort \n");

   scanf("\n");
   rin = getchar();

   if(rin == 'a') exit(0);

   if(rin == 'y') {

      for(i=0; i < track_num; i++) free((tind+i)->tp);

      free(tind);

      for(i=0; i < frame_num; i++){

         foo = fo + i;

         for(j=0; j < foo->obj_num; j++){

             ob = (foo->objs) + j;

             for(k=0; k < ob->fet->feature_num; k++){

                fpts = (ob->fet->fpt) + k;
                fpts->track_id = 0;

             }

          }

      }

      if(izz && zone){

         for(i=0; i < npar; i++) {
             free(zone[i]->zlat);
             free(zone[i]);
         }
         free(zone);

      }

      if(iadpt && add){

         for(i=0; i < npar; i++) free(add[i]);
         free(add);

      }
      free(ddmax);
      free(odmax);
      free(pphimax);

      goto re_initial;

   }

/* calculate the initial value of the cost fuction */

   if(!iadpt){

      totgij = 0.0;

      for(i=0; i < track_num; i++){

          tr = tind + i;
          if(tr->num_point_track != 0){

             for(k=1; k < frame_num-1; k++){ 

                 tpsi[0] = (tr->tp) + k-1;
                 tpsi[1] = (tr->tp) + k;
                 tpsi[2] = (tr->tp) + k+1;

                 iim = ((fo + k - 1)->nmiss > (fo + k)->nmiss) ? (fo + k - 1)->nmiss : (fo + k)->nmiss;
                 if(iim > npar - 1) iim = npar - 1;

                 totgij += devn(fo, tpsi[0], tpsi[1], tpsi[2], *(pphimax + iim));


             }
 
          }

      }

      printf("the initial value of the cost function is %f \n", totgij);

   }

/* begin iteration */

   fi_fl = 1;
   bi_fl = 1;

/* print warning if fewer than three frames */

   if(frame_num <= 3) printf("***WARNING***, insufficient number of frames chosen.\n");

/* begin the iteration */

   tot_count = 1;


   if(frame_num > 3){

      while(fi_fl ==1 || bi_fl == 1){

        printf("%d \n", tot_count);

        fi_count = 0;

        if(fi_fl == 1) {


           if(pproc) {

              tind = track_split(fo, tind, frame_num);

              tind = tr_zonal_filter(add, zone, fo, tind, &ddmax, &pphimax, frame_num, -1, pproc, 'f', npar);

           }


           fel_mge(add, zone, fo, tind, ddmax, pphimax, frame_num, fi_count, npar);

        }

        bi_count = 0;

        if(bi_fl == 1 && tot_count < tot_term) {

           if(pproc) {

              tind = track_split(fo, tind, frame_num);

              tind = tr_zonal_filter(add, zone, fo, tind, &ddmax, &pphimax, frame_num, -1, pproc, 'b', npar);

           }


           bel_mge(add, zone, fo, tind, ddmax, pphimax, frame_num, bi_count, npar);


        }

        ++tot_count;

        if(tot_count > tot_term) fi_fl = bi_fl = 0;

      }

      tind = track_split(fo, tind, frame_num);


      if(izz && zone){

         for(i=0; i < npar; i++) {
             free(zone[i]->zlat);
             free(zone[i]);
         }
         free(zone);

      }

      if(iadpt && add){

         for(i=0; i < npar; i++) free(add[i]);
         free(add);

      }

miss_zone:

/* need to combine the next two options */

      printf("Do you want a missing frame search, 'y' or 'n'\n");
      scanf("\n");

      if(getchar() == 'y') tind = tr_miss_frame(fo, tind, dmax, phimax, frame_num);



     printf("Do you want to apply zonal upper bound displacements, \r\n"
            "and/or adaptive track smoothness, 'y' or 'n'?         \n\n");

     scanf("\n");
     if(getchar() == 'y') tind = tr_zonal_filter(add, zone, fo, tind, &ddmax, &pphimax, frame_num, 0, 0, 'a', npar);



/* calculate the final value of the cost fuction */

      totgij = 0.0;

      for(i=0; i < track_num; i++){

         tr = tind + i;
         if(tr->num_point_track != 0){

           for(k=1; k < frame_num-1; k++){ 

              tpsi[0] = (tr->tp) + k-1;
              tpsi[1] = (tr->tp) + k;
              tpsi[2] = (tr->tp) + k+1;

              iim = ((fo + k - 1)->nmiss > (fo + k)->nmiss) ? (fo + k - 1)->nmiss : (fo + k)->nmiss;
              if(iim > npar - 1) iim = npar - 1;

              totgij += devn(fo, tpsi[0], tpsi[1], tpsi[2], *(pphimax + iim));

           }

         }

      }

      printf("the final value of the cost function is %f \n", totgij);

      re_tr_num = 0;
      for(i=0; i < track_num; i++){

          tr = tind + i; 
          if(tr->num_point_track > 0){

            ++re_tr_num;

            for(j=0; j < frame_num; j++){

                tps = (tr->tp)+j;

                if(tps->feature_id != -1){

                  ob = ((fo+j)->objs) + (tps->object_id) - 1;
                  ((ob->fet->fpt) + (tps->feature_id) - 1)->track_id = i + 1;

                } 

            }

          }   

      }

      tdump = open_file(file_exist(fout3, "r", APPS), "w");

      trackd(tind, track_num, frame_num, tdump, re_tr_num);

      close_file(tdump, fout3);

replot:

      printf("do you want any plots of the tracks\r\n"
             "input the min. number of points in a track for plotting\r\n"
             "input '0' for no further ploting\n");

      scanf("%d", &mnpt);

      if(mnpt > 0){

         ++fp_count;

         tot_fet = 0;

         for(i=0; i < track_num; i++){

             tr = tind + i;
             if(tr->num_point_track >= mnpt) tot_fet += tr->num_point_track;

         }

         printf("do you want area of interest (input 1) or whole area (input 0)\n");
         scanf("%d", &izm);

         featd(fo, tind, tot_fet, frame_num, fp_count);
         goto replot;

      } 

      if(rd == 'y') {

         for(i=0; i< track_num; i++) free((tind+i)->tp);
         free(tind);
         return NULL;

      }

rplottr:

      printf("what track do you want plotted, input '0' for none\n");

      scanf("%d", &t_id);

      if(t_id > track_num){

         printf("track does not exist, track_num = %d \r\n"
                "choose track number again\n", track_num);

         goto rplottr;

      }

      if(t_id > 0){

         ++fp_count;
         tr = tind + t_id - 1;
         tot_fet = tr->num_point_track;

         printf("do you want area of interest (input 1) or whole area (input 0)\n");
         scanf("%d", &izm);

         featd(fo, tind, tot_fet, frame_num, fp_count);
         goto rplottr;

      }

   }


   printf("do you want to repeat the calculation with a different parametization \r\n"
          "using the existing tracks, calculation, y or n. Note, a smaller dmax  \r\n"
          "should not be chosen unless re-initialization occurs.                 \n\n");

   scanf("\n");
   rc = getchar();

   if(rc == 'y') {

     for(i=0; i < track_num; i++) free((tind+i)->tp);

     free(tind);

     for(i=0; i < npar; i++) *(odmax + i) = *(ddmax + i);

     for(i=0; i < frame_num; i++){

        foo = fo + i;

        for(j=0; j < foo->obj_num; j++){

            ob = (foo->objs) + j;

            for(k=0; k < ob->fet->feature_num; k++){

               fpts = (ob->fet->fpt) + k;
               fpts->track_id = 0;

            }

         }

     }

     if(izz && zone){

         for(i=0; i < npar; i++) {
             free(zone[i]->zlat);
             free(zone[i]);
         }
         free(zone);

     }

     if(iadpt && add){

         for(i=0; i < npar; i++) free(add[i]);
         free(add);

     }

     free(ddmax);
     free(odmax);
     free(pphimax);

     goto new_param;

   }

   if(izz && zone){

       for(i=0; i < npar; i++) {
           free(zone[i]->zlat);
           free(zone[i]);
       }
       free(zone);

   }

   if(iadpt && add){

       for(i=0; i < npar; i++) free(add[i]);
       free(add);

   }

   free(ddmax);
   free(odmax);
   free(pphimax);

   geo_init = 0;

   return tind;

}
