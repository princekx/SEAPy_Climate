#include <Stdio.h>
#include <stdlib.h>
#include "st_fo.h"
#include "st_track.h"
#include "complex.h"
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"

/* function to compute the displacement distance along a track */


float measure(struct feature_pts * , struct feature_pts * );

extern int tom;

void disp_filter(struct tot_tr *all_tr, int trackn, int *totf)
{

    int i, j;
    int tn=0;
    int wf;

    float dspfil;
    float dispm=0.;
    float scc=0.;

    struct tot_tr *altr;
    struct feature_pts f1, f2;
    struct fet_pt_tr *at, *at1;

    if(tom == 'e'){

       printf("Current distance measure is Euclidean, you must enter the\r\n"
              "Euclidean filter distance when prompted.                 \n\n");

       scc = 1.0;

    }

    else if(tom == 'g'){

       printf("Current distance measure is Geodesic, you must enter the \r\n"
              "angular filter distance when prompted.                   \n\n");

       scc = FP_PI;

    }

    printf("What filter distance do you require, in the current units\r\n"
           "and for the current distance measure.                    \n\n");

    scanf("%f", &dspfil);

    dspfil *= scc;

    printf("Do you want to filter according to displacement distance along        \r\n"
           "the track or seperation distance between start and finish, 't' or 's' \n\n");

    scanf("\n");
    wf = getchar();

    if(!(wf == 't' || wf == 's')){

       printf("***WARNING***, incorrect filter specifier, no distance filtering performed.\n\n");

       return;

    }

    *totf = 0;

    for(i=0; i<trackn; i++){

        altr = all_tr + i;

        dispm = 0.;

        if(wf == 't' && altr->trpt){

           at = altr->trpt;
           at1 = at + 1;

           for(j=1; j < altr->num; j++){

               (f1.x).xy = at->xf;
               (f1.y).xy = at->yf;

               (f2.x).xy = at1->xf;
               (f2.y).xy = at1->yf;

               dispm += measure(&f1, &f2);

               ++at;
               ++at1;      

           }



        }

        else if(wf == 's' && altr->trpt){

            at = altr->trpt;
            at1 = at + altr->num - 1;
            (f1.x).xy = at->xf;
            (f1.y).xy = at->yf;

            (f2.x).xy = at1->xf;
            (f2.y).xy = at1->yf;
            dispm = measure(&f1, &f2);

        }

        if(dispm < dspfil){

	   for(j=0; j < altr->num; j++) {if(altr->trpt) free((altr->trpt + j)->add_fld);}
           altr->num = 0;
           free(altr->trpt);
           altr->trpt = NULL;

        }

        else ++tn;

        *totf += altr->num;

    }

    printf("***INFORMATION***, current number of tracks is now %d\n\n", tn);

    return;

}
