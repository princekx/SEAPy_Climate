#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file_cat_out.h"
#include "mem_er.h"
#include "st_fo.h"
#include "file_handle.h"
#include "splice.h"
#include "m_values.h"

#define   LARGE_DIST   1.0e+6
#define   LARGE_TRNG   1000000

/* function to merge tracks from a common object based on seperation distance */

void meantrd(FILE * , struct tot_tr * , int , int );
float measure(struct feature_pts * , struct feature_pts * );

extern int nf, nfld;
extern char *fext;
extern int iext;
extern int tom;
extern int tf;

void track_merge(struct tot_tr *alltr, int trackn, int trtyp)
{
    int i, j;
    int ifmn=LARGE_TRNG, ifmx=0;
    int itrm=0;
    int iavg=0;
    int nold=0;

    char tmerge[MAXCHR];

    float mrad=0.0;
    float dd=0.0;
    float dmin=LARGE_DIST;

    FILE *tsf=NULL;

    struct tot_tr *altr1=NULL, *altr2=NULL;
    struct fet_pt_tr *atr1=NULL, *atr2=NULL;
    struct feature_pts fpts1, fpts2;

    tf = 4;

/* check end points */

    printf("Do you want track end point matching, 'y' or 'n'.\n\n");
    scanf("\n");
    if(getchar() == 'y'){

       printf("What is the search radius required for track end point merging?\n\n");
       scanf("%f", &mrad);

       printf("When merging tracks average over common points or keep common points as they are,\r\n"
              "input '0' for keeping points as they are or '1' for averaging.                   \n\n");
       scanf("%d", &iavg);

/* determine frame range */

       for(i=0; i < trackn; i++){
           altr1 = alltr + i;
           if(altr1->num > 0){
              atr1 = altr1->trpt;
              if(atr1->fr_id < ifmn) ifmn = atr1->fr_id;
              atr1 = altr1->trpt + altr1->num - 1;
              if(atr1->fr_id > ifmx) ifmx = atr1->fr_id;
           }
       }

       printf("****INFORMATION****, time step range is %d -- %d\n\n", ifmn, ifmx);

       for(i=0; i < trackn; i++){

           altr1 = alltr + i;

/* match end of track to start of another track */

           if(altr1->num == 0) continue;

           atr1 = altr1->trpt + altr1->num - 1;
           if(atr1->fr_id == ifmx) continue;

           (fpts1.x).xy = atr1->xf;
           (fpts1.y).xy = atr1->yf;

           dmin = LARGE_DIST;

           for(j=0; j < trackn; j++){

              altr2 = alltr + j;

              if(i == j || altr2->num == 0) continue;

              atr2 = altr2->trpt;
              if(atr2->fr_id == ifmn) continue;

              if(atr1->fr_id == atr2->fr_id && atr1->obj_id == atr2->obj_id){

                 (fpts2.x).xy = atr2->xf;
                 (fpts2.y).xy = atr2->yf;
                 dd = measure(&fpts1, &fpts2);
                 if(tom == 'g') dd /= FP_PI;
                 if(dd < dmin && dd <= mrad) {dmin = dd; itrm = j;}
              }

           }

/* merge tracks */

           if(dmin <= mrad) {
              altr1->imt = 1;
              altr2 = alltr + itrm;
              nold = altr1->num;

              if(iavg){
                 altr1->imt = 1;
                 altr1->num += altr2->num - 1;
                 altr1->trpt = (struct fet_pt_tr * )realloc_n(altr1->trpt, (altr1->num)*sizeof(struct fet_pt_tr));
                 mem_er((altr1->trpt == NULL) ? 0 : 1, (altr1->num)*sizeof(struct fet_pt_tr));
                 atr1 = altr1->trpt + nold - 1;
                 atr2 = altr2->trpt;
                 atr1->xf = (atr1->xf + atr2->xf) * 0.5;
                 atr1->yf = (atr1->yf + atr2->yf) * 0.5;
                 atr1->zf = (atr1->zf + atr2->zf) * 0.5;

/* the following averages maybe redundent but are kept for safety */

                 atr1->gwthr = (atr1->gwthr + atr2->gwthr) * 0.5;
                 atr1->tend = (atr1->tend + atr2->tend) * 0.5;
                 atr1->sh_an = (atr1->sh_an + atr2->sh_an) * 0.5;
                 atr1->area = (atr1->area + atr2->area) * 0.5;
                 atr1->wght = (atr1->wght + atr2->wght) * 0.5;
                 atr1->vec[0] = (atr1->vec[0] + atr2->vec[0]) * 0.5;
                 atr1->vec[1] = (atr1->vec[1] + atr2->vec[1]) * 0.5;
                 atr1->or_vec[0] = (atr1->or_vec[0] + atr2->or_vec[0]) * 0.5;
                 atr1->or_vec[1] = (atr1->or_vec[1] + atr2->or_vec[1]) * 0.5;

                 if(nf){
                    for(j=0; j < nfld; j++){
                        *(atr1->add_fld + j) = (*(atr1->add_fld + j) + *(atr2->add_fld + j)) * 0.5;
                    }
                 }

                 for(j=1; j < altr2->num; j++){
                     memcpy(altr1->trpt + nold + j - 1, altr2->trpt + j, sizeof(struct fet_pt_tr));
                 }

              }
              else {
                 altr1->imt = 2;
                 altr1->num += altr2->num;
                 altr1->trpt = (struct fet_pt_tr * )realloc_n(altr1->trpt, (altr1->num)*sizeof(struct fet_pt_tr));
                 mem_er((altr1->trpt == NULL) ? 0 : 1, (altr1->num)*sizeof(struct fet_pt_tr));
                 for(j=0; j < altr2->num; j++){
                     memcpy(altr1->trpt + nold + j, altr2->trpt + j, sizeof(struct fet_pt_tr));
                 }

              }


              if(nf){
                for(j=0; j < altr2->num; j++) free((altr2->trpt + j)->add_fld);
              }
              altr2->num = 0;
              free(altr2->trpt);
              altr2->trpt = NULL;

              --i;

           }

       }

    }

    strncpy(tmerge, FILTRS, MAXCHR);
    if(iext) strcpy(strstr(tmerge, EXTENSION), fext);
    strcat(tmerge, "_merged");

    tsf = open_file(tmerge, "w");

    meantrd(tsf, alltr, trackn, trtyp);  
 
    close_file(tsf, tmerge);


    return;
} 
