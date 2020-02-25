#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "complex.h"
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"

#define  MAXCHR   100

/* function to convert time information to frame ID's */

long int new_time(long int , int );

extern int nf;

void time2fr_id(struct tot_tr *all_tr, int trnum)
{

   int i, j;
   int nstime=0;
   int ifr=0;
   int num=0;
   int nff=0;
   int *ioff=NULL;
   int itrunc=0;
   int ndel1=0, ndel2=0, nnew=0;

   long int *stime=NULL, *etime=NULL;
   long int ts=0;
   long int ta=0, tb=0, tst=0, ten=0;
   long int t1=0, t2=0;

   char timef[MAXCHR];

   FILE *ftime=NULL;

   struct tot_tr *altr=NULL;
   struct fet_pt_tr *atr=NULL, *att=NULL;


   printf("Input the time step in the correct units, eg. Hrs.\n\n");
   scanf("%ld", &ts);


   printf("Read starting times from file, 'f', or read from the keyboard, 'k'\n\n");
   scanf("\n");
   if(getchar() == 'f'){

      printf("What is the start time file to be read?\n\n");
      scanf("%s", timef);

      ftime = open_file(timef, "r");

      fscanf(ftime, "%d\n", &nstime);

      stime = (long int *)calloc(nstime, sizeof(long int));
      mem_er((stime == NULL) ? 0 : 1, nstime * sizeof(long int));

      etime = (long int *)calloc(nstime, sizeof(long int));
      mem_er((etime == NULL) ? 0 : 1, nstime * sizeof(long int));

      ioff = (int *)calloc(nstime, sizeof(int));
      mem_er((ioff == NULL) ? 0 : 1, nstime * sizeof(int));

      for(i=0; i < nstime; i++) fscanf(ftime, "%ld %ld %d\n", stime + i, etime + i, ioff + i);

      close_file(ftime, timef);

   }

   else {
       
      printf("How many starting times do you want to specify?\n\n");
      scanf("%d", &nstime);

      stime = (long int *)calloc(nstime, sizeof(long int));
      mem_er((stime == NULL) ? 0 : 1, nstime * sizeof(long int));

      etime = (long int *)calloc(nstime, sizeof(long int));
      mem_er((etime == NULL) ? 0 : 1, nstime * sizeof(long int));

      ioff = (int *)calloc(nstime, sizeof(int));
      mem_er((ioff == NULL) ? 0 : 1, nstime * sizeof(int));

      for(i=0; i < nstime; i++){
          printf("What is the next start and end time and offset: ");
          scanf("%ld %ld %d", stime + i, etime + i, ioff + i);
          printf("\n");

      }


   }

   for(i=0; i < trnum; i++){

       altr = all_tr + i;

       ta = 0;
       t1 = altr->trpt->time;
       t2 = (altr->trpt + altr->num - 1)->time;

       for(j=0; j < nstime; j++){

           if((t1 >= *(stime + j) && t1 <= *(etime + j)) && (t2 >= *(stime + j) && t2 <= *(etime + j))) {

              ta = *(stime + j);
              nff = *(ioff + j);

              break;

           }

       }

       if(!(int)ta){
          if(!itrunc){
/*             printf("Truncate track to match times, 'y' or 'n' or exit\n\n");
             scanf("\n");
             if(getchar() == 'y') itrunc = 1;
             else {
                printf("****ERROR****, track time for tracks ID %d not in specified time periods.\n\n", altr->trid);
                exit(1);
             } */
	     itrunc = 1;
          }
          if(itrunc) {
/* find period which overlaps with track */
             ta = 0;
             for(j=0; j < nstime; j++){
                tst = *(stime + j);
                ten = *(etime + j);
                if((tst >= t1 && tst <= t2) || (ten >= t1 && ten <= t2)){
                   ta = *(stime + j);
                   tb = *(etime + j);
                   nff = *(ioff + j);                   
                }
             }

/* truncate track */

             if(!(int)ta) {
                if(nf){
                   for(j=0; j < altr->num; j++) {
                      atr = altr->trpt + j;
                      free(atr->add_fld); 
                   }
                }
                free(altr->trpt);
                altr->num = 0;
                continue;
             }
             else {

                nnew = altr->num;

                ndel1 = ndel2 = 0;
                for(j=0; j < altr->num; j++){
                   atr = altr->trpt + j;
                   if(atr->time < ta) {
                      if(nf) free(atr->add_fld);
                      ++ndel1;
                      --nnew;
                   }
                   if(atr->time > tb){
                      if(nf) free(atr->add_fld);
                      ++ndel2;
                      --nnew;
                   }
                }
                if(ndel1){
                   att = altr->trpt;
                   for(j=ndel1; j < altr->num; j++){
                       atr = altr->trpt + j;
                       memcpy(att, atr, sizeof(struct fet_pt_tr));
                       ++att;
                   }
                }
                if(ndel1 || ndel2) {
                   if(!nnew) {free(altr->trpt);}
                   else{
                      altr->trpt = (struct fet_pt_tr *)realloc_n(altr->trpt, nnew * sizeof(struct fet_pt_tr));
                      mem_er((altr->trpt == NULL) ? 0 : 1, nnew*sizeof(struct fet_pt_tr));
                   }
                   altr->num = nnew;
                   altr->time = ta;
                }
             }
          }
       }

       ifr = 1 + nff;

       while(ta < altr->time){
            ta = new_time(ta, ts);
            ++ifr;
       }

       num = 0;
       for(j=0; j < altr->num; j++) {

          atr = altr->trpt + j;
          atr->fr_id = ifr + num;
          num += 1 + atr->nfm;

       }

       altr->time = 0;

   }


   free(stime);
   free(etime);
   free(ioff);

   return;

}
