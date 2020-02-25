#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"

#define MAXCHR  150
#define TOLWT   1.0e-6

void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);
float ran3(int * );

int noheader=0;
extern int nfld, nff;
extern int *nfwpos;

extern float sum_wt;
extern float sum_per;

int main(void)

{
   int i, j;
   int nr=0, ns=0;
   int idum;
   int trnum, ntr=0;
   int gpr, ipr;
   int nsamp=0;
   int *isamp=NULL;
   int nu=0, umo=0, umo5=0;
   int awt = 0;
   int iin=0;
   int iwr=0;
   int ntrs1=0, ntrs2=0;
   int trnn=0;
   int *ich=NULL;
   int offs=0;
   int nln=0;

   long int pl1, pl2;

   float alat=0.0, alng=0.0;
   float ran;
   float sc1, sc2;

   float nwt=0.0, swt=0.0;

   FILE *fin=NULL;
   FILE *fout1=NULL, *fout2=NULL;
   FILE *mfil=NULL;

   char trname[MAXCHR];
   char filout1[MAXCHR], filout2[MAXCHR];
   char odir[MAXCHR];
   char fnum[10];
   char mfnm[]="missing.dat";
   
   struct tot_tr *all_tr=NULL, *atr1=NULL, *atr2=NULL, *allout=NULL;
   struct tot_tr *allout2=NULL;
   struct fet_pt_tr *fpt1=NULL, *fpt2=NULL;


   printf("What is the raw track file to be used?\n\n");
   scanf("%s", trname);

   strcpy(filout1, strstr(trname, "tr_trs"));

   fin = fopen(trname, "r");
   if(!fin){
      printf("File %s cannot be opened for read\n\n", trname);
      exit(1);
   }


   all_tr = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

   fclose(fin);

   printf("The number of tracks in this ensemble = %d\n\n", trnum);

   printf("What is the ouput directory for samples?\n\n");
   scanf("%s", odir);

   printf("Do you want to sample with '1' or without '0' replacement.\n\n");
   scanf("%d", &iwr);

   if(!iwr){
      printf("How many tracks are required for the twin samples?\n\n");
      scanf("%d %d", &ntrs1, &ntrs2);
      printf("What are the scaling periods for the two samples?\n");
      scanf("%f %f", &sc1, &sc2);

   }
   
   trnn = (iwr) ? trnum : ntrs1;



   printf("Input a seed value, must be negative\n\n");
   scanf("%d", &idum);

   printf("How many samples do you want?\n\n");
   scanf("%d", &nsamp);

   strcpy(filout1, odir);
   strcat(filout1, "/seed.dat");
   fout1 = fopen(filout1, "w");
   if(!fout1){
      printf("File %s cannot be opened for read\n\n", filout1);
      exit(1);
   }
   fprintf(fout1, "%d\n", idum);
   fclose(fout1);

   printf("What offset is required for numbering?\n");
   scanf("%d", &offs);

   if(nsamp <= 0){
      printf("***ERROR***, number of samples must be positive and non-zero\n\n");
     exit(1);

   }

   allout = (struct tot_tr * )calloc(trnn, sizeof(struct tot_tr));
   mem_er((allout == NULL) ? 0 : 1, trnn * sizeof(struct tot_tr));

   if(!iwr){

      allout2 = (struct tot_tr * )calloc(trnum-trnn, sizeof(struct tot_tr));
      mem_er((allout2 == NULL) ? 0 : 1, (trnum - trnn) * sizeof(struct tot_tr));

   }

   isamp = (int *)calloc(trnum, sizeof(int));
   mem_er((isamp == NULL) ? 0 : 1, trnum * sizeof(int));

   ich = (int *)calloc(trnum, sizeof(int));
   mem_er((ich == NULL) ? 0 : 1, trnum * sizeof(int));

   if(iwr){
     mfil = fopen(mfnm, "w");
     if(!fin){
        printf("File %s cannot be opened for read\n\n", mfnm);
        exit(1);
     }
     fprintf(mfil, "%d %d\n", nsamp, trnum);
   }

   while(ns < nsamp){

      nwt = 0.0;
      swt = 0.0;

      strcpy(filout1, odir);
      strcat(filout1, strstr(trname, "/tr_trs"));
      sprintf(fnum, ".%06d", ns+offs);
      strcat(filout1, fnum);

      printf("Creating sample %d\n", ns+offs);

      if(!iwr){
         strcpy(filout2, filout1);
         strcat(filout1, "_A");
         strcat(filout2, "_B");
      }


      fout1 = fopen(filout1, "w");
      if(!fout1){
         printf("Can't open file %s for write\n\n", filout1);
         exit(1);
      }

      if(!iwr){

         fout2 = fopen(filout2, "w");
         if(!fout2){
            printf("Can't open file %s for write\n\n", filout2);
            exit(1);
         }

      }

      for(i=0; i < trnum; i++) *(ich + i) = 0;


      nr = 0;
      nu = umo = umo5 = 0;
      while(nr < trnn){

         ran = ran3(&idum) * (float) (trnum - 1) + 0.5;

         atr1 = all_tr + (int)ran;

         if(!iwr && *(ich + (int)ran)) continue;

         *(ich + (int)ran) = 1;

         ++(*(isamp + (int)ran));

         atr2 = allout +nr;

         atr2->num = atr1->num;
         atr2->awt = atr1->awt;


         atr2->trpt = (struct fet_pt_tr * )calloc(atr2->num, sizeof(struct fet_pt_tr));
         mem_er((atr2->trpt == NULL) ? 0 : 1, atr2->num * sizeof(struct fet_pt_tr));

         for(i=0; i < atr2->num; i++) {

             fpt1 = atr1->trpt + i;
             fpt2 = atr2->trpt + i;
             memcpy(fpt2, fpt1, sizeof(struct fet_pt_tr));
             if(nff){
                fpt2->add_fld = (float *)calloc(nfld, sizeof(float));
		mem_er((fpt2->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                memcpy(fpt2->add_fld, fpt1->add_fld, nfld * sizeof(float));
             }
 
             if(fpt2->wght > TOLWT) {
                nwt += 1.0;
                swt += fpt2->wght;
             }

         }

         ++nr;

      }


      if(iwr) {pl1 = ftell(mfil); fprintf(mfil, "%5d %5d\n", ns+offs, 0); nln = 0;}
      for(i=0; i < trnum; i++){
         if(*(isamp+i) == 0) {
             ++nu;
             if(iwr){
                fprintf(mfil, "%d ", i);
                ++nln;
                if(nln == 10) {fprintf(mfil, "\n"); nln = 0;}
             }
          }
          if(*(isamp+i) > 1) ++umo;
          if(*(isamp+i) > 5) ++umo5;
      }
      if(iwr){
         if(nln) fprintf(mfil, "\n");
         pl2 = ftell(mfil);
         fseek(mfil, pl1, FSTART);
         fprintf(mfil, "%5d %5d\n", ns+offs, nu);
         fseek(mfil, pl2, FSTART);  
      }

      printf("Number of tracks not used = %d\n", nu);
      printf("Number of tracks used more than once = %d\n", umo);
      printf("Number of tracks used more than 5 times = %d\n", umo5);

      if(!iwr) {sum_wt = sc1 * swt / nwt; sum_per = sc1;}

      meantrd(fout1, allout, trnn, 's', gpr, ipr, alat, alng);

      for(i=0; i < trnn; i++){
          atr2= allout + i;
          if(nff){
             for(j=0; j < atr2->num; j++) free((atr2->trpt + j)->add_fld);
          }
          free(atr2->trpt);
          atr2->trpt = NULL;
          atr2->num = 0;
      }

      fclose(fout1);
      fout1 = NULL;

      if(!iwr){

         nr = 0;
         nwt = 0.0;
         swt = 0.0;

         for(i=0; i < trnum; i++){
             if(! *(ich + i)){
                atr1 = all_tr + i;
                atr2 = allout2 +nr;

                atr2->num = atr1->num;
                atr2->awt = atr1->awt;


                atr2->trpt = (struct fet_pt_tr * )calloc(atr2->num, sizeof(struct fet_pt_tr));
                mem_er((atr2->trpt == NULL) ? 0 : 1, atr2->num * sizeof(struct fet_pt_tr));

                for(j=0; j < atr2->num; j++) {

                   fpt1 = atr1->trpt + j;
                   fpt2 = atr2->trpt + j;
                   memcpy(fpt2, fpt1, sizeof(struct fet_pt_tr));
                   if(nff){
                      fpt2->add_fld = (float *)calloc(nfld, sizeof(float));
		      mem_er((fpt2->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                      memcpy(fpt2->add_fld, fpt1->add_fld, nfld * sizeof(float));
                   }

                   if(fpt2->wght > TOLWT){
                     nwt += 1.0;
                     swt += fpt2->wght;
                   }

                }

                ++nr;

                if(nr == ntrs2) break;
                
             }
         }

         sum_wt = sc2 * swt / nwt;
         sum_per = sc2;

         meantrd(fout2, allout2, nr, 's', gpr, ipr, alat, alng);

         for(i=0; i < trnum - trnn; i++){
             atr2 = allout2 + i;
             if(nff){
                for(j=0; j < atr2->num; j++) free((atr2->trpt + j)->add_fld);
             }
             free(atr2->trpt);
             atr2->trpt = NULL;
             atr2->num = 0;
         }

         fclose(fout2);
         fout2 = NULL;
      }

      ++ns;

      for(i=0; i < trnum; i++) *(isamp + i) = 0;
/*      memset(isamp, NULL, trnum*sizeof(int)); */

   }

   if(mfil) fclose(mfil);

   for(i=0; i < trnum; i++){
      atr1 = all_tr + i;
      if(nff){
        for(j=0; j < atr1->num; j++) free((atr1->trpt + j)->add_fld);
      }
      free(atr1->trpt);
   }
   free(all_tr);
   free(allout);
   free(allout2);
   free(isamp);
   free(ich);

   return 0;

}
