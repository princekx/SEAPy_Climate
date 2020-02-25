#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "complex.h"
#include "splice.h"
#include "st_fo.h"
#include "m_values.h"

/* function to filter tracks according to their propogation characteristics */

float measure(struct feature_pts * , struct feature_pts * );

extern float period;

void prop_filt(struct tot_tr *all_tr, int trackn, int *totf)

{

    int i, j;
    int tn=0;

    float hp=period/2.;
    float ttint;
    float ang;
    float xd, xdn;
    float xcomp=0.0, ycomp=0.0;
    float a1=0.0, a2=0.0;
    float spd=0.0;

    double xx, yy;

    struct tot_tr *altr;
    struct fet_pt_tr *at, *at1;
    struct feature_pts f1, f2;

    if(trackn <= 0 || all_tr == NULL){

       printf("****WARNING****, no tracks have been loaded, \r\n"
              "               load tracks and try again.    \n\n");

       return;

    }

    printf("Input the angular boundaries in degrees, a1 and a2\r\n"
           "in the range (-180, 180).                         \n\n");
    scanf("%f %f", &a1, &a2);


    *totf = 0;

    for(i=0; i<trackn; i++){

        altr = all_tr + i;
        at = altr->trpt;
        at1 = at + 1;


        xx = 0.0;
        yy = 0.0;

        if(altr->num > 1){

           for(j=1; j<altr->num; j++){

               xd = at1->xf - at->xf;
               xdn = fabs(xd);

               if(xdn > hp){

                  if(at->xf < hp) xd -= period;

                  else xd += period;

               } 

               (f1.x).xy = at->xf;
               (f1.y).xy = at->yf;

               (f2.x).xy = at1->xf;
               (f2.y).xy = at1->yf;
               ttint = (at->nfm + 1);

               spd = measure(&f1, &f2)/ttint;

               (f2.y).xy = at->yf;

               xcomp = measure(&f1, &f2)/ttint;

               if(xd < 0.) xcomp *= -1.;

               (f2.y).xy = at1->yf;
               (f2.x).xy = at->xf;

               ycomp = measure(&f1, &f2)/ttint;
               if(at1->yf - at->yf < 0.) ycomp *= -1.;


               xx += xcomp;
               yy += ycomp;

               ++at;
               ++at1; 


           }

           
           ang = atan2(yy, xx)/FP_PI;

           if((a2 - ang)*(ang - a1) < 0.0){

              for(j=0; j < altr->num; j++) free((altr->trpt + j)->add_fld);
              altr->num = 0;
              free(altr->trpt);
              altr->trpt = NULL;

           }

           else ++tn;


           *totf += altr->num;



         }

         else{

           for(j=0; j < altr->num; j++) {if(altr->trpt) free((altr->trpt + j)->add_fld);}
           altr->num = 0;
           free(altr->trpt);
           altr->trpt = NULL;

         }




    }

    printf("***INFORMATION***, current number of tracks is now %d\n\n", tn);


    return;

}
