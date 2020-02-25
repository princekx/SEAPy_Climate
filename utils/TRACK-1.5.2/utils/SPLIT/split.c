#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "splice.h"

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );

int noheader=0;

int main(int argc, char **argv)
{

   int i, j;
   int gpr, ipr, trnum;
   int ntr=0;
   int nsplit=0, nn=0;


   float alat, alng;

   FILE *fin=NULL;
   FILE *fout=NULL;

   char outn[200], num[6];

   struct tot_tr *all_tr=NULL, *atr=NULL;
   struct tot_tr *nall_tr=NULL;

   if(argc != 2){
      printf("Usage: split [track file]\n\n");
      exit(1);
   }

   fin = fopen(argv[1], "r");
   if(!fin){
      printf("***ERROR***, file %s does not exist for 'r'\n\n", argv[1]);
      exit(1);

   }

   
   all_tr = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

   fclose(fin);

   for(i=0; i<trnum; i++){

       atr = all_tr + i;

       if(atr->trid == 1) ++nsplit;

   }

   printf("%d\n", nsplit);

   for(i=0; i<nsplit; i++){

       strcpy(outn, argv[1]);
       sprintf(num, "_%02d", i+1);
       strcat(outn, num);


       fout = fopen(outn, "w");
       if(!fout){
          printf("***ERROR***, unable to openw file %s for 'w'\n\n", outn);
          exit(1);
       }

       ntr = 0;
       nn = 0;
       nall_tr = NULL;

       for(j=0; j < trnum; j++){

           atr = all_tr + j;

           if(atr->trid == 1) ++nn;

           if(i+1 == nn){
              if(!nall_tr) nall_tr = atr;
              ++ntr;
           }

       }


       meantrd(fout, nall_tr, ntr, 's', gpr, ipr, alat, alng);

       fclose(fout);
       fout = NULL;

   }


   return 0;
}
