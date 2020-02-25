#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "splice.h"
#include "file_handle.h"


#define SCALE   -1.5
#define TOLWT   1.0e-6
#define  FP_PI   0.017453292519943295   
#define  FP_PI2  1.570796326794896619


/* function to write to file the track data obtained by combining
   different sets of track data                                        */

void sincos(double , double * , double * );

extern int aniso;
extern int nfld, nff;
extern int *nfwpos;


int meantrd_tanh(FILE *tsf, struct tot_tr *all_tr, int tr_count, int fld, int *istr, int *iend, float *ival, int num, int iraw, int ireg, float lat1, float lat2, float lng1, float lng2)

{

    int i, j, k;
    int trc=0;
    int npp=0;
    int inreg=1;
    int ifr=0;

    static int ntrack=0;
    static int nottr=0;

    float wt=0.0;
    float xcen, ycen;
    double xx=0.0, yy=0.0, zz=0.0;
    double xp, yp, zp;
    double s1, c1, s2, c2;
    double dist, dmin;
    static float nwt=0.0, nmth=0.0;
    static double swt=0.0;

    struct tot_tr *altr;
    struct fet_pt_tr *atr;

    for(i=0; i < num; i++){

        if(ival[i] > TOLWT) nmth += 1.0;

    }

    if(ireg){
       printf("****WARNING****, note current spatial sampling does not work across periodic boundaries or poles.\n\n");
       xcen = 0.5 * (lng1 + lng2);
       ycen = 0.5 * (lat1 + lat2);
       xcen *= FP_PI;
       ycen = FP_PI2 - ycen * FP_PI;
       sincos(xcen, &s1, &c1);
       sincos(ycen, &s2, &c2);
       xx = s2 * c1;
       yy = s2 * s1;
       zz = c2;
    }


    for(i=0; i < tr_count; i++){

        altr = all_tr + i;

        if(!altr->num) continue;

        if(ireg){
           inreg = 0;
           dmin = 10000000.0;
           for(j=0; j < altr->num; j++){
               atr = (altr->trpt) + j;
               if(((lng2 - atr->xf) * (atr->xf - lng1) >= 0.0) &&
                  ((lat2 - atr->yf) * (atr->yf - lat1) >= 0.0)    ) {

                  inreg = 1;
                  xcen = atr->xf * FP_PI;
                  ycen = FP_PI2 - atr->yf * FP_PI;
                  sincos(xcen, &s1, &c1);
                  sincos(ycen, &s2, &c2);
                  xp = s2 * c1;
                  yp = s2 * s1;
                  zp = c2;
                  dist = xx * xp + yy * yp + zz * zp;
                  if(dist < dmin) {dmin = dist; ifr = atr->fr_id;}
               }
           }

           if(ireg == 2 && inreg){
              wt = 0.0;
              for(j=0; j<num; j++){
                 if((iend[j] - ifr) * (ifr - istr[j]) >= 0) {
                   if(!iraw) wt = ival[j];
                   else {
                      if(ival[j] < TOLWT) wt = 0.0;
                      else wt = ival[j];
                   }
                 }

              }

           }
        }

        if(!inreg) continue;


        npp = 0;
        ++trc;

        fprintf(tsf, "TRACK_ID  %d\n", i+1);
        fprintf(tsf, "POINT_NUM  %d\n", altr->num);

        if(fld == 's'){

           if(aniso == 'y'){

              for(j=0; j < altr->num; j++){

                  atr = (altr->trpt) + j;

                  if(ireg == 0 || ireg == 1){

                     wt = 0.0;
                     for(k=0; k<num; k++){

                         if((iend[k] - atr->fr_id) * (atr->fr_id - istr[k]) >= 0) {
                            if(!iraw) wt = ival[k];
                            else {
                               if(ival[k] < TOLWT) wt = 0.0;
                               else {
                                 wt = ival[k];
                                 swt += wt;
                                 nwt += 1.0;
                                 ++npp;
                               }
                            }
                            break;
                         }

                     }

                  }

                  else if(ireg == 2) {
                     if(wt > 0.0){
                        swt += wt;
                        nwt += 1.0;
                        ++npp;
                     }
                  }

                  fprintf(tsf, "%d %f %f %e ", atr->fr_id, atr->xf, atr->yf, atr->zf);
                  if(nfld){
                     fprintf(tsf, "& ");                  
                     for(k=0; k < nfld; k++) fprintf(tsf, "%e & ", *(atr->add_fld + k));
                  }

                  fprintf(tsf, "%f %f %f %f", atr->sh_an, atr->or_vec[0], atr->or_vec[1], wt);

                  if(atr->nfm)
                     fprintf(tsf, "%d\n", atr->nfm);
                  else
                     fprintf(tsf, "\n");

              }

           }

           else {

              for(j=0; j < altr->num; j++){

                  atr = (altr->trpt) + j;

                  if(ireg == 0 || ireg == 1){

                     for(k=0; k<num; k++){
                         wt = 0.0;
                         if((iend[k] - atr->fr_id) * (atr->fr_id - istr[k]) >= 0) {
                            if(!iraw) wt = ival[k];
                            else {
                               if(ival[k] < TOLWT) wt = 0.0;
                               else {
                                 wt = ival[k];
                                 swt += wt;
                                 nwt += 1.0;
                                 ++npp;
                               }
                            }
                            break;
                         }


                     }

                  }

                  else if(ireg == 2){
                     if(wt > 0.0){
                        swt += wt;
                        nwt += 1.0;
                        ++npp;
                     }
                  }

                  fprintf(tsf, "%d %f %f %e ", atr->fr_id, atr->xf, atr->yf, atr->zf);
                  if(nfld){
                     fprintf(tsf, "& ");                  
                     for(k=0; k < nfld; k++) fprintf(tsf, "%e & ", *(atr->add_fld + k));
                  }                  

                  if(atr->nfm)
                     fprintf(tsf, "%f %d\n", wt, atr->nfm);
                  else
                     fprintf(tsf, "%f\n",  wt);

              }
  
           }

        }

        else if(fld == 'v'){

           for(j=0; j < altr->num; j++){

               atr = (altr->trpt) + j;

               fprintf(tsf, "%d %f %f %e %e", atr->fr_id, atr->xf, atr->yf, atr->zf, atr->gwthr);
               fprintf(tsf, " %e %e\n", atr->vec[0], atr->vec[1]);

           }

        }

        if((float) npp / (float) altr->num > 0.5) ++ntrack;
        else ++nottr;

    } 

    if(iraw){
       printf("Tracks with +ve weights %d, \r\n"
              "Tracks with zero weights %d \r\n"
              "Periods %f\n\n", ntrack, nottr, nmth);
    }

    return trc;

}

/*

void sincos(double ang, double *sn, double *cn)

{


    *sn = sin(ang);
    *cn = cos(ang);

    return;

}

*/
