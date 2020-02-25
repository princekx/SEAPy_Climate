#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "mem_er.h"
#include "st_fo.h"
#include "geo_values.h"
#include "m_values.h"
#include "grid.h"
#include "splice.h"

/* function to compute the speeds of objects */

float measure(struct feature_pts * , struct feature_pts * );
int tr_overlap(int * , int , int , int , int );

extern int tom;
extern float period;
extern GRID *gr;

struct tot_tr *phase(struct tot_tr *all_tr, int trackn, int *trphn, int fr1, int fr2)

{

    int i, j, isc, trn=0;
    int in=0, nn=0;
    int imiss=0;
    int num=0;
    
    long int is1=0, is2=0;

    float tint, sc=1., xd, xdn;
    float hp=period/2.;
    float xmin, xmax;
    float ttint;

    struct tot_tr *phsp=NULL, *altr, *phtr;
    struct fet_pt_tr *trfp, *at, *at1;
    struct feature_pts f1, f2;

    if(trackn <= 0 || all_tr == NULL){

       printf("****WARNING****, no tracks have been loaded, \r\n"
              "               load tracks and try again.    \n\n");

       return NULL;

    }

    xmin = *(gr->xgrid);
    xmax = *(gr->xgrid+gr->ix-1);

    printf("current distance measure is ");
    if(tom == 'e') printf("euclidean\n");
    else if(tom == 'g') printf("geodesic\n");

    printf("do you wish to change the distance measure, 'y' or 'n'\r\n"
           "(use geodesic 'g' for a spherical domain)\n");

    scanf("\n");
    if(getchar() == 'y') {

      printf("Input 'e' for Euclid or 'g' for geodesic\n");
      scanf("\n");
      tom = getchar();

    }

    printf("what is the time interval\n");
    scanf("%f", &tint);

    printf("do you want to scale the phase speed, 'y' or 'n'\r\n"
           "e.g. speed relative to the Earths surface, scale with radius of the Earth\n");

    scanf("\n");
    if(getchar() == 'y'){

      printf("What type of value do you want      \r\n"
             "Input '0' for user value            \r\n"
             "Input '1' for Earths radius in Km   \r\n");
      scanf("%d", &isc);

      if(isc == 0 ) {

         printf("input scaling\n");
         scanf("%f", &sc);

      }

      else if (isc == 1) sc = EARTH_RADIUS;

    }


    for(i=0; i < trackn; i++){

        altr = all_tr + i;
        at = altr->trpt;

        in = 0;

        if(altr->num > 1) in = tr_overlap(&nn, fr1, fr2, at->fr_id, (at + altr->num - 1)->fr_id);

        if(in){

           ++trn;

/* allocate memory */

           if(trn == 1) {

              phsp = (struct tot_tr * )malloc_initl(sizeof(struct tot_tr));
              mem_er((phsp == NULL) ? 0 : 1, sizeof(struct tot_tr));

           }
           
           else {

              phsp = (struct tot_tr * )realloc_n(phsp, trn*sizeof(struct tot_tr));
              mem_er((phsp == NULL) ? 0 : 1, trn*sizeof(struct tot_tr));

           }

           phtr = phsp + trn - 1;

           phtr->num = nn;
           phtr->awt = altr->awt;
           phtr->tr_sp = i;
           phtr->time = altr->time;
	   phtr->trid = altr->trid;

           phtr->trpt = (struct fet_pt_tr * )calloc(phtr->num, sizeof(struct fet_pt_tr));
           mem_er((phtr == NULL) ? 0 : 1, phtr->num * sizeof(struct fet_pt_tr));

/* compute phase speeds */

           trfp = phtr->trpt;
           at = altr->trpt;
           at1 = at + 1;

           num = 0;

           for(j=1; j < altr->num; j++){

               if(!imiss){

                  if(at->nfm){

                     imiss = 1;
                     printf("***WARNING***, missing frames have been detected whilst \r\n"
                            "               computing speeds and velocities. Speeds  \r\n"
                            "               and velocities are corrected accordingly.\n\n");

                  }


               }

               is1 = (long int)(at->fr_id - fr1) * (long int)(fr2 - at->fr_id);
               is2 = (long int)(at1->fr_id - fr1) * (long int)(fr2 - at1->fr_id);


               if(is1 >= 0 || is2 >= 0){

                 ++num;

                 trfp->fr_id = at->fr_id;
                 trfp->time = at->time;
                 trfp->wght = at->wght;
                 xd = at1->xf - at->xf;
                 xdn = fabs(xd);

                 if(xdn > hp){

                   if(at->xf < hp)

                      xd -= period;

                   else

                      xd += period;

                 } 

                 
                 trfp->xf = at->xf + xd * 0.5;

                 if(! gr->prty){
                   if(trfp->xf < xmin) trfp->xf += period;
                   else if(trfp->xf > xmax) trfp->xf -= period;
                 }
                 trfp->yf = (at->yf + at1->yf) * 0.5;

                 (f1.x).xy = at->xf;
                 (f1.y).xy = at->yf;

                 (f2.x).xy = at1->xf;
                 (f2.y).xy = at1->yf;

                 ttint = tint * (at->nfm + 1);

                 trfp->zf = sc * measure(&f1, &f2)/ttint;

                 if(at1->zf < ADD_CHECK && at->zf < ADD_CHECK){
                    trfp->tend = (at1->zf - at->zf) / ttint;
                    trfp->gwthr = trfp->tend / (0.5 * (at1->zf + at->zf));

/* may need to change this to corect evaluation */

                    (f2.y).xy = at->yf;

                    trfp->vec[0] = sc * measure(&f1, &f2)/ttint;

                    if(xd < 0.) trfp->vec[0] *= -1.;

                    (f2.y).xy = at1->yf;
                    (f2.x).xy = at->xf;

                    trfp->vec[1] = sc * measure(&f1, &f2)/ttint;
                    if(at1->yf - at->yf < 0.) trfp->vec[1] *= -1.;

                 }
                else {
                    trfp->tend = ADD_UNDEF;
                    trfp->gwthr = ADD_UNDEF;
                    trfp->zf = ADD_UNDEF;
                    trfp->vec[0] = ADD_UNDEF;
                    trfp->vec[1] = ADD_UNDEF;
                 }
                            
                 ++trfp;

               }

               ++at;
               ++at1; 

           }

           phtr->num = num;

           in = 0;         

        }

    }

    *trphn = trn;

    return phsp;

}
