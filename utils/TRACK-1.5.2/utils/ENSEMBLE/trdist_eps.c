#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "splice.h"
#include "st_fo.h"
#include "m_values.h"
#include "sqt.h"

#define DLARGE  10000000.0

float measure(struct feature_pts *  , struct feature_pts * );
double ortho_dist(struct fet_pt_tr * , struct fet_pt_tr * , int , int * , int , int , VEC *);

float trdist_eps(struct tot_tr *tr1, struct tot_tr *tr2, float disth, long int is1, long int is2, int *numm, int *it1, int *it2, int isep, int ty_match, int ntrpt, int isepty)
{

    int i;
    int ist=0;
    int ifnd=0;
    int ndf=0;

    int it1s=0, it1e=0, it2s=0, it2e=0;
    int num=0, numt=0;

    long int t1, t2;

    float dist=0.0, ddd=0.0, newd=0.0;
    double dp;

    struct feature_pts fpts1, fpts2;
    struct fet_pt_tr *fp1=NULL, *fp2=NULL;
    
    VEC v1, v2;
    VEC visec;

    *numm = 0;

/* find position along track for start and end times */


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

    numt = it1e - it1s + 1;
    num = it2e - it2s + 1;
    *it1 = it1s;
    *it2 = it2s;

    if(num != numt){
       printf("***ERROR***, track section lengths dont match for TRACK_ID's %d %d\n\n", tr1->trid, tr2->trid);
       exit(1);
    }

    t1 = tr1->trpt->time;
    t2 = tr2->trpt->time;

    if(ty_match == 1){
/*       if(t1 >= t2) num = ntrpt;
       else num = ntrpt - it1s;
       if(num < ntrpt) num = 0; */
       num = ntrpt;
    }
    else if (ty_match == 2){
/*       if(t1 >= t2) num = 1;
       else num = 0; */
       num = 1;
    }

    if(num > numt) num = numt;

/* calculate track seperation */

    if(isep){

       dist = DLARGE;

       if(!isepty){
          for(i=0; i < num; i++){

             fp1 = tr1->trpt + it1s + i;
             fp2 = tr2->trpt + it2s + i;

             (fpts1.x).xy = fp1->xf;
             (fpts1.y).xy = fp1->yf;

             (fpts2.x).xy = fp2->xf;
             (fpts2.y).xy = fp2->yf;

             ddd = measure(&fpts1, &fpts2);

             if(ty_match == 1){
               if(ddd / FP_PI > disth) {num = 0; break;}
             }

             if(ddd < dist) dist = ddd;
         }
      
      }
      
      else {
      
         ifnd = 1;
	 ndf = 0;
         for(i=0; i < num; i++){
	     fp1 = tr1->trpt + it1s + i;
	     if(ist < tr2->num){
	        ddd = ortho_dist(fp1, tr2->trpt, tr2->num, &ist, ifnd, 1, &visec);
		if(ifnd && ist > 0) {
                   fp1 = tr2->trpt + it2s + i;
	           v1.x = fp1->pp[0];
	           v1.y = fp1->pp[1];
	           v1.z = fp1->pp[2];
                   fp2 = tr2->trpt + ist - 1;
	           v2.x = fp2->pp[0];
                   v2.y = fp2->pp[1];
	           v2.z = fp2->pp[2];
                   dp = dotp(&v1, &v2);
	           newd = acos(dp);
	           if(newd / FP_PI > 10.0) {num=0; break;}
		}
	        if(ist > 0){
	           ifnd = 0;
	           ++ndf;
	           if(ty_match == 1){
                     if(ddd / FP_PI > disth) {num = 0; break;}
                   }

                   if(ddd < dist) dist = ddd;
	     
	        }
	     
	     }
	 
	 }

      }

    }

    else {

      if(!isepty){
         for(i=0; i < num; i++){

             fp1 = tr1->trpt + it1s + i;
             fp2 = tr2->trpt + it2s + i;

             (fpts1.x).xy = fp1->xf;
             (fpts1.y).xy = fp1->yf;

             (fpts2.x).xy = fp2->xf;
             (fpts2.y).xy = fp2->yf;

             ddd = measure(&fpts1, &fpts2);

             if(ty_match == 1){
                if(ddd / FP_PI > disth) {num = 0; break;}
             }

             dist += ddd;
         } 
	 
	 if(num > 0) dist /= num;
       
      }
      else {
      
         ifnd = 1;
	 ndf = 0;
         for(i=0; i < num; i++){
	     fp1 = tr1->trpt + it1s + i;
	     v1.x = fp1->pp[0];
	     v1.y = fp1->pp[1];
	     v1.z = fp1->pp[2];
	     if(ist < tr2->num){
	        ddd = ortho_dist(fp1, tr2->trpt, tr2->num, &ist, ifnd, 1, &visec);
		if(ifnd && ist > 0) {
                   fp1 = tr2->trpt + it2s + i;
	           v1.x = fp1->pp[0];
	           v1.y = fp1->pp[1];
	           v1.z = fp1->pp[2];
                   fp2 = tr2->trpt + ist - 1;
	           v2.x = fp2->pp[0];
                   v2.y = fp2->pp[1];
	           v2.z = fp2->pp[2];
                   dp = dotp(&v1, &v2);
	           newd = acos(dp);
	           if(newd / FP_PI > 10.0) {num=0; break;}
		}

	        if(ist > 0){
		   ifnd = 0;
	           ++ndf;
	           if(ty_match == 1){
                     if(ddd / FP_PI > disth) {num = 0; break;}
                   }

                   dist += ddd;
	     
	        }
	     
	     }
	 
	 }
	 
         if(ndf > 0) dist /= ndf;

      }

    }

    if(num > 0){
      *numm = numt;
      dist /= FP_PI;
    }
    else {
      *numm = 0;
      dist = DLARGE;
    }

    return dist;

}
