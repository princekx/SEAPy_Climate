#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "splice.h"
#include "m_values.h"
#include "file_handle.h"
#include "file_cat_out.h"

/* function to sample and compute track smoothness and mean displacement
   distances.                                                            */

double euclid_dev(struct feature_pts * , struct feature_pts * , struct feature_pts * );
double geod_dev(struct feature_pts * , struct feature_pts * , struct feature_pts * );
float measure(struct feature_pts * , struct feature_pts * );
void meantrd(FILE * , struct tot_tr * , int , int );

extern int tom;
extern int tf;
extern float w1, w2;
extern int iext;
extern char *fext;

extern int nfld, nf;
extern int *nfwpos;

void track_sample_dev(struct tot_tr* trs, int trnum)
{

   int i, j;
   int tsmp=1;
   int ipos=0;
   int wrttr=0;
   int npt=0;
   int ntmp=0, noff=0, nstoff=0;
   int nft=0, nfldt=0, *nfwpost=NULL;
   int id1=0, id2=0, id3=0;
   int idf=0, itadd=0;
   int ifreset=0, itrid=0;
   int itsmax=0;

   FILE *fdev=NULL, *fdist=NULL;
   FILE *fout=NULL;

   float dev=0.0, dist=0.0;
   float dis1, dis2;

   struct tot_tr *altr=NULL, *antr=NULL;
   struct tot_tr *newtr=NULL;
   struct fet_pt_tr *atr=NULL, *anr=NULL, *at=NULL;
   struct feature_pts fp0, fp1, fp2;

   char trout[MAXCHR];

   tf = 4;

   printf("**************************************************\n"
          "*           TRACK TIMESTEP SAMPLING              *\n"
          "**************************************************\n\n");
	  
   printf("****WARNING****, missing time steps currently not supported, do you want to continue, 'y' or 'n'\n\n");
   scanf("\n");
   if(getchar() != 'y') exit(1);

   printf("Leave any additional fields as they are or replace with displacements and deviations,\r\n"
          "input '0' to leave alone or '1' for replacement.                                     \n\n");
   scanf("%d", &itadd);

   if(itadd){

     nft = nf;
     nfldt = nfld;
     if(nft){
       nfwpost = (int *)calloc(nf, sizeof(int));
       mem_er((nfwpost == NULL) ? 0 : 1, nf * sizeof(int));
       memcpy(nfwpost, nfwpos, nf * sizeof(int));
     }

     nf = 2;
     nfld = 2;

     if(nfwpos) free(nfwpos);
     nfwpos = (int *)calloc(nf, sizeof(int));
     mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));
     nfwpos[0] = nfwpos[1] = 0;

   }

   printf("What sampling of the tracks is required starting from time step 1?\r\n"
          "Note a sampling value of 1 is every time step.                    \n\n");
   scanf("%d", &tsmp);
   if(tsmp < 1){
      printf("****ERROR****, incorrect sampling value, must be >=1.\n\n");
      exit(1);
   }

   printf("Input the track timestep offset.\n\n");
   scanf("%d", &noff);
   if(noff < 0){
      printf("****ERROR****, offset must be positive.\n\n");
      exit(1);
   }
   
   printf("What is the track truncation time step? Use large value for no truncation\n\n");
   scanf("%d", &itsmax);

   if(!(trs->time)){
      printf("Reset time step id, 'y' or 'n'\n\n");
      scanf("\n");
      if(getchar() == 'y') ifreset = 1;

      printf("What time step reset offset is required?\n\n");
      scanf("%d", &nstoff);
   }

   printf("What are the cost function weights, w1 and w2?\n\n");
   scanf("%f %f", &w1, &w2);

   if(fabs(w1 + w2 - 1.0) > 1.0e-4){

      printf("****WARNING****, weights do not sum to unity.\n\n");

   }

   w1 *= 0.5;

   printf("Do you want to output sampled tracks? 'y' or 'n'.\n\n");
   scanf("\n");
   if(getchar() == 'y') wrttr = 1;

   if(wrttr){
      newtr = (struct tot_tr *)calloc(trnum, sizeof(struct tot_tr));
      mem_er((newtr == NULL) ? 0 : 1, trnum * sizeof(struct tot_tr));
      for(i=0; i < trnum; i++){
          altr = trs + i;
          antr = newtr + i;
	  memcpy(antr, altr, sizeof(struct tot_tr));
	  antr->num = 0;
          ntmp = (altr->num / tsmp) + tsmp;
          antr->trpt = (struct fet_pt_tr * )calloc(ntmp, sizeof(struct fet_pt_tr));
          mem_er((antr->trpt == NULL) ? 0 : 1, ntmp * sizeof(struct fet_pt_tr));
      }

      strncpy(trout, FPTTRS, MAXCHR);
      if(iext) strcpy(strstr(trout, EXTENSION), fext);
      strncat(trout, "_filt", MAXCHR);
   }

   fdev = open_file("dev.dat", "w");
   fdist = open_file("dist.dat", "w");

   for(i=0; i < trnum; i++){

       altr = trs + i;
       atr = altr->trpt;
       if(wrttr) antr = newtr + i;
       npt = 0; 
       
       if(!(altr->num)) continue;
       if(atr->fr_id > itsmax) continue;

/* find first point */

       ipos = altr->num + 1;

       for(j=0; j < altr->num; j++) {
           at = atr + j;
           idf = at->fr_id - 1 - noff;
           if(!(idf % tsmp) && idf >= 0) {ipos = j; itrid = (idf / tsmp) + 1; break;}
       }

       if(ipos >= altr->num || at->fr_id > itsmax) continue;

       fp0.gwky = 0;
       (fp0.x).xy = at->xf;
       (fp0.y).xy = at->yf;
       if(wrttr) {
          memcpy(antr->trpt + npt, at, sizeof(struct fet_pt_tr));
          id1 = npt; ++npt; antr->num = npt;
          anr = antr->trpt + id1;
          if(ifreset) anr->fr_id = itrid - nstoff;
          anr->add_fld = (float *)calloc(nfld, sizeof(float));
          mem_er((anr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
          if(itadd) *(anr->add_fld) = *(anr->add_fld + 1) = ADD_UNDEF;
          else memcpy(anr->add_fld, at->add_fld, nfld * sizeof(float));
          ++itrid;
       }

       ipos += tsmp;
       at += tsmp;
       
       if(ipos >= altr->num || at->fr_id > itsmax) continue;
       fp1.gwky = 0;
       (fp1.x).xy = at->xf;
       (fp1.y).xy = at->yf;
       if(wrttr) {
          memcpy(antr->trpt + npt, at, sizeof(struct fet_pt_tr)); 
          id2 = npt; ++npt; antr->num = npt;
          anr = antr->trpt + id2;
          if(ifreset) anr->fr_id = itrid - nstoff;
          anr->add_fld = (float *)calloc(nfld, sizeof(float));
          mem_er((anr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
          if(itadd)*(anr->add_fld) = *(anr->add_fld + 1) = ADD_UNDEF;
          else memcpy(anr->add_fld, at->add_fld, nfld * sizeof(float));
          ++itrid;
       }

       ipos += tsmp;
       at += tsmp;
       
       if(ipos >= altr->num || at->fr_id > itsmax) continue;
       fp2.gwky = 0;
       (fp2.x).xy = at->xf;
       (fp2.y).xy = at->yf;
       if(wrttr) {
          memcpy(antr->trpt + npt, at, sizeof(struct fet_pt_tr)); 
          id3 = npt; ++npt; antr->num = npt;
          anr = antr->trpt + id3;
          if(ifreset) anr->fr_id = itrid - nstoff;
          anr->add_fld = (float *)calloc(nfld, sizeof(float));
          mem_er((anr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
          if(itadd) *(anr->add_fld) = *(anr->add_fld + 1) = ADD_UNDEF;
          else memcpy(anr->add_fld, at->add_fld, nfld * sizeof(float));
          ++itrid;
       }

       ipos += tsmp;
       at += tsmp;

       switch(tom){
           case 'e':
             dev = euclid_dev(&fp0, &fp1, &fp2);
             break;
           case 'g':
             dev = geod_dev(&fp0, &fp1, &fp2);
             break;
           default:
             printf("***error***, incorrect key used for type of measure in %s\n", __FILE__);
             exit(1);
       }

       dis1 = measure(&fp0, &fp1) / FP_PI;
       dis2 = measure(&fp1, &fp2) / FP_PI;

       dist = 0.5 * (dis1 + dis2);

       fprintf(fdev, "%f %f\n", dist, dev);
       fprintf(fdist, "%f %f %f\n", 0.5 * ((fp0.x).xy + (fp1.x).xy), 0.5 * ((fp0.y).xy + (fp1.y).xy), dis1);
       fprintf(fdist, "%f %f %f\n", 0.5 * ((fp1.x).xy + (fp2.x).xy), 0.5 * ((fp1.y).xy + (fp2.y).xy), dis2);


       if(wrttr && itadd){
          anr = antr->trpt + id1;
          *(anr->add_fld) = dis1;
          *(anr->add_fld + 1) = ADD_UNDEF;
          anr = antr->trpt + id2;
          *(anr->add_fld) = dis2;
          *(anr->add_fld + 1) = dev;
       }

       while(ipos < altr->num){
             fp0 = fp1;
             fp1 = fp2;
             ++id2;
             dis1 = dis2;
             fp2.gwky = 0;
             (fp2.x).xy = at->xf;
             (fp2.y).xy = at->yf;
             if(wrttr) {
                memcpy(antr->trpt + npt, at, sizeof(struct fet_pt_tr)); 
                anr = antr->trpt + npt;
                ++npt; antr->num = npt;
                if(ifreset) anr->fr_id = itrid - nstoff;
                anr->add_fld = (float *)calloc(nfld, sizeof(float));
                mem_er((anr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                if(itadd) *(anr->add_fld) = *(anr->add_fld + 1) = ADD_UNDEF;
                else memcpy(anr->add_fld, at->add_fld, nfld * sizeof(float));
                ++itrid;
             }

             ipos += tsmp;
             at += tsmp;
	     
	     if(ipos >= altr->num || at->fr_id > itsmax) break;

             switch(tom){
                 case 'e':
                    dev = euclid_dev(&fp0, &fp1, &fp2);
                    break;
                 case 'g':
                    dev = geod_dev(&fp0, &fp1, &fp2);
                    break;
                 default:
                    printf("***error***, incorrect key used for type of measure in %s\n", __FILE__);
                    exit(1);
             }

             dis2 = measure(&fp1, &fp2) / FP_PI;
             dist = 0.5 * (dis1 + dis2);

             fprintf(fdev, "%f %f\n", dist, dev);
             fprintf(fdist, "%f %f %f\n", 0.5 * ((fp1.x).xy + (fp2.x).xy), 0.5 * ((fp1.y).xy + (fp2.y).xy), dis2);

             if(wrttr && itadd){
               anr = antr->trpt + id2;
               *(anr->add_fld) = dis2;
               *(anr->add_fld + 1) = dev;
             }

       }

       if(wrttr && itadd){
          anr = antr->trpt + antr->num - 1;
          *(anr->add_fld) = ADD_UNDEF;
          *(anr->add_fld + 1) = ADD_UNDEF;
       } 

   }
 


   close_file(fdev, "dev.dat");
   close_file(fdist, "dev.dat");

   if(wrttr){
      fout = open_file(trout, "w");
      meantrd(fout, newtr, trnum, 's');  
      close_file(fout, trout);
   }

   if(wrttr){
      for(i=0; i < trnum; i++) {
          antr = newtr + i;
          if(nf){
             for(j=0; j < antr->num; j++){
                 atr = antr->trpt + j;
                 free(atr->add_fld);
             }

          }
          free(antr->trpt);
      }
      free(newtr);
   }

   if(itadd){
      nf = nft;
      nfld = nfldt;
      if(nft){
        memcpy(nfwpos, nfwpost, nf * sizeof(int));
        free(nfwpost);
      }
   }

   return;
}
