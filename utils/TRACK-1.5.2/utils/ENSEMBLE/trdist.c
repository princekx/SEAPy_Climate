#include <stdio.h>
#include <stdlib.h>
#include "splice.h"
#include "st_fo.h"
#include "m_values.h"

#define DLARGE  10000000.0
#define DTEST   100000.0

float measure(struct feature_pts *  , struct feature_pts * );
void alloc_point(struct fet_pt_tr * , struct feature_pts * , int , int );

float trdist(struct tot_tr *tr1 , struct tot_tr *tr2 , float disth, long int is1 , long int is2, int *numm, int *it1, int *it2, int isep, int ty_match, int numst, int upos1, int upos2, int ipos1, int ipos2)
{

    int i;

    int it1s=0, it1e=0, it2s=0, it2e=0;
    int num=0, numt=0, nn=0;
    
    float dist=0.0, ddd=0.0;

    struct feature_pts fpts1, fpts2;
    struct fet_pt_tr *fp1=NULL, *fp2=NULL;


/* find position along track for start and end times */

    if(tr1->time && tr2->time){    

       for(i=0; i < tr1->num; i++){

           fp1 = tr1->trpt + i;

           if(fp1->time == is1) it1s = i;
           if(fp1->time == is2) it1e = i;

       } 



       for(i=0; i < tr2->num; i++){

           fp2 = tr2->trpt + i;

           if(fp2->time == is1) it2s = i;
           if(fp2->time == is2) it2e = i;

       }  


    }

    else {

       for(i=0; i < tr1->num; i++){

           fp1 = tr1->trpt + i;

           if(fp1->fr_id == is1) it1s = i;
           if(fp1->fr_id == is2) it1e = i;

       } 



       for(i=0; i < tr2->num; i++){

           fp2 = tr2->trpt + i;

           if(fp2->fr_id == is1) it2s = i;
           if(fp2->fr_id == is2) it2e = i;

       } 

    }

    numt = it1e - it1s + 1;
    num = it2e - it2s + 1;
    *it1 = it1s;
    *it2 = it2s;


    if(num != numt){
       printf("***ERROR***, track section lengths dont match for TRACK_ID's %d %d\n\n", tr1->trid, tr2->trid);
       exit(1);
    }
    
    if(ty_match) num = numst;
    
    if(num > numt) num = numt;

/* calculate track seperation */

    if(isep){

       dist = DLARGE;

       for(i=0; i < num; i++){

          fp1 = tr1->trpt + it1s + i;
          fp2 = tr2->trpt + it2s + i;
	  
	  alloc_point(fp1, &fpts1, upos1, ipos1);
	  alloc_point(fp2, &fpts2, upos2, ipos2);	  

          if((fpts1.x).xy > ADD_CHECK || (fpts2.x).xy > ADD_CHECK) continue;
	  
          ddd = measure(&fpts1, &fpts2);
          if(ty_match){
            if(ddd / FP_PI > disth) {num = 0; break;}
          }	  
          if(ddd < dist) dist = ddd;
       }

    }

    else {

      for(i=0; i < num; i++){

          fp1 = tr1->trpt + it1s + i;
          fp2 = tr2->trpt + it2s + i;
	  
	  alloc_point(fp1, &fpts1, upos1, ipos1);
	  alloc_point(fp2, &fpts2, upos2, ipos2);
  
	  if((fpts1.x).xy > ADD_CHECK || (fpts2.x).xy > ADD_CHECK) continue;

          ++nn;
	  
	  ddd = measure(&fpts1, &fpts2);
	  
          if(ty_match){
            if(ddd / FP_PI > disth) {num = 0; nn = 0; break;}
          }	  
	  
          dist += ddd;
       } 

       if(nn > 0) dist /= nn;
       else dist = DLARGE;

    }

    *numm = numt;
    if(dist < DTEST) dist /= FP_PI;

    return dist;

}

void alloc_point(struct fet_pt_tr *fp, struct feature_pts *fpts, int upos, int ipos)
{

     if(!upos){
        (fpts->x).xy = fp->xf;
        (fpts->y).xy = fp->yf;
     }
     else if(ipos >= 0){
        (fpts->x).xy = fp->add_fld[ipos];
        (fpts->y).xy = fp->add_fld[ipos + 1]; 
     }
     else{
        printf("****ERROR**** no positional information available for this choice of variable.\n\n");
	exit(1);
     }
     
     return;

}
