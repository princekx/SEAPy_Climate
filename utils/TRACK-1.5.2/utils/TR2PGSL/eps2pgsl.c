#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include "splice.h"
#include "mem_er.h"

#define NEPS  21
#define LENSTR 100


/* program to convert groups of track files into comma seperated columns for input into postgres */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);

int aniso;
int nff, nfld;
int *nfwpos;

int main(int argc, char **argv)
{
   int i, j, k, l;
   int nmch=0;
   int nfil=0;

   int tr_count=0;
   int gpr, ipr;
   int gprt, iprt;

   int ittyp=0, iclstr=0, itrid=0;

   int yr, mn, dy, hr;

   float alat, alng;
   float alatt, alngt;

   long int datei, datef;

   char lsstr[LENSTR];
   char trfil[LENSTR];
   char fcsv[]="alltr.csv";
   char ffld[5];
   char date[20];

   FILE *fcsvw=NULL;
   FILE *ftr=NULL;

   struct dirent *entry;

   DIR *dir=NULL;

   struct tot_tr *alltr=NULL, *trr=NULL;
   struct fet_pt_tr *atr;

   if(argc != 3){
      printf("Usage: tr2pgsl [Date] [Field]\n\n");
      exit(1);
   } 

   sscanf(argv[1], "%ld", &datei);
   sscanf(argv[2], "%s", ffld);

   sprintf(lsstr, "./%ld_%s", datei, ffld);

   dir = opendir(lsstr);
   if(dir == NULL){
      printf("****ERROR****, cannot open directory %s\n\n", lsstr);
      exit(1);
   }
   entry = readdir( dir );
   while(entry != NULL){
     if( strncmp( entry->d_name, "trmatch", 7 ) == 0 ) ++nfil;
     entry = readdir( dir );
   }
   closedir(dir);

   printf("%d\n", nfil);

   if(chdir(lsstr)){
     printf("****ERROR****, can't change directory to %s\n", lsstr);
     exit(1);
   }

   fcsvw = fopen(fcsv, "w");

/* write matched clusters */

   for(i=0; i < nfil; i++){
      iclstr = i + 1;
      if(i+1 < 10) sprintf(trfil, "trmatch_cntl_tr000%d", i+1);
      else if(i+1 < 100) sprintf(trfil, "trmatch_cntl_tr00%d", i+1);
      else if(i+1 < 1000) sprintf(trfil, "trmatch_cntl_tr0%d", i+1);
      else {
         printf("****ERROR****, not configured for this many clusters.\n\n");
         exit(1);
      }

      ftr = fopen(trfil, "r");
      if(!ftr){
         printf("****ERROR****, can't open file %s for 'r'\n\n", trfil);
         exit(1);
      }
      alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
      fclose(ftr);

      if(!i){gprt = gpr; iprt = ipr; alatt = alat; alngt = alng;}
      else{
         if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
            printf("****ERROR****, track files are incompatable.\n\n");
            exit(1); 
         }
      }

      gprt = gpr;
      iprt = ipr;
      alatt = alat;
      alngt = alng;

      for(j=0; j < tr_count; j++){
         ++itrid;
         trr = alltr + j;
         if(! (trr->trid)) ittyp = 0;
         else if(trr->trid == 1) ittyp = -1;
         else ittyp = trr->trid - 1;

         for(k=0; k < trr->num; k++){
            atr = (trr->trpt) + k; 
/* reformat date */
            yr = atr->time / 1000000;
            mn = (atr->time - yr * 1000000) / 10000;
            dy = (atr->time - yr * 1000000 - mn * 10000) / 100;
            hr = atr->time - yr * 1000000 - mn * 10000 - dy * 100;
            sprintf(date, "%4d-%02d-%02d %02d:00:00", yr, mn, dy, hr);
            fprintf(fcsvw, "%2d, %4d, %4d, \"%18s\", %f, %f, %e, ", ittyp, iclstr, itrid, date, atr->xf, atr->yf, atr->zf);
            for(l=0; l < 12; l++) fprintf(fcsvw, "%e, ", atr->add_fld[l]);
            fprintf(fcsvw, "%e\n", atr->add_fld[12]);
         }

      }

/* free memory */

      for(j=0; j < tr_count; j++) {
         trr = alltr + j;
         for(k=0; k < trr->num; k++) {
             atr = (trr->trpt) + k;
             free(atr->add_fld);
         }
         free(trr->trpt);
      }

      free(alltr);

   }

/* write mean tracks */

   ittyp = -2;

   for(i=0; i < nfil; i++){
      iclstr = i + 1;
      if(i+1 < 10) sprintf(trfil, "trmatch_cntl_tr000%d_mean", i+1);
      else if(i+1 < 100) sprintf(trfil, "trmatch_cntl_tr00%d_mean", i+1);
      else if(i+1 < 1000) sprintf(trfil, "trmatch_cntl_tr0%d_mean", i+1);
      else {
         printf("****ERROR****, not configured for this many clusters.\n\n");
         exit(1);
      }

      ftr = fopen(trfil, "r");
      if(!ftr){
         printf("****WARNING****, can't open file %s for 'r'\n\n", trfil);
         continue;
      }
      alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
      fclose(ftr);

      if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
         printf("****ERROR****, track files are incompatable.\n\n");
         exit(1);
      }

      ++itrid;

      for(k=0; k < alltr->num; k++){
          atr = (alltr->trpt) + k;
/* reformat date */
          yr = atr->time / 1000000;
          mn = (atr->time - yr * 1000000) / 10000;
          dy = (atr->time - yr * 1000000 - mn * 10000) / 100;
          hr = atr->time - yr * 1000000 - mn * 10000 - dy * 100;
          sprintf(date, "%4d-%02d-%02d %02d:00:00", yr, mn, dy, hr);
          fprintf(fcsvw, "%2d, %4d, %4d, \"%18s\", %f, %f, %e, ", ittyp, iclstr, itrid, date, atr->xf, atr->yf, atr->zf);
          for(l=0; l < 12; l++) fprintf(fcsvw, "%e, ", atr->add_fld[l]);
          fprintf(fcsvw, "%e\n", atr->add_fld[12]);
      }


      for(j=0; j < tr_count; j++) {
         trr = alltr + j;
         for(k=0; k < trr->num; k++) {
             atr = (trr->trpt) + k;
             free(atr->add_fld);
         }
         free(trr->trpt);
      }

      free(alltr);

   }


/* write unmatched tracks */

   iclstr = 0;

   for(i=0; i < NEPS; i++){
      if(!i) ittyp = -1;
      else ittyp = i;

      if(i+1 < 10) sprintf(trfil, "trnomatch_ens000%d", i+1);
      else if(i+1 < 100) sprintf(trfil, "trnomatch_ens00%d", i+1); 

      ftr = fopen(trfil, "r");
      if(!ftr){
         printf("****ERROR****, can't open file %s for 'r'\n\n", trfil);
         exit(1);
      }
      alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
      fclose(ftr);

      if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
         printf("****ERROR****, track files are incompatable.\n\n");
         exit(1);
      }

      gprt = gpr;
      iprt = ipr;
      alatt = alat;
      alngt = alng;

      for(j=0; j < tr_count; j++){
         ++itrid;
         trr = alltr + j;

         for(k=0; k < trr->num; k++){
            atr = (trr->trpt) + k;
/* reformat date */
            yr = atr->time / 1000000;
            mn = (atr->time - yr * 1000000) / 10000;
            dy = (atr->time - yr * 1000000 - mn * 10000) / 100;
            hr = atr->time - yr * 1000000 - mn * 10000 - dy * 100;
            sprintf(date, "%4d-%02d-%02d %02d:00:00", yr, mn, dy, hr);
            fprintf(fcsvw, "%2d, %4d, %4d, \"%18s\", %f, %f, %e, ", ittyp, iclstr, itrid, date, atr->xf, atr->yf, atr->zf);
            for(l=0; l < 12; l++) fprintf(fcsvw, "%e, ", atr->add_fld[l]);
            fprintf(fcsvw, "%e\n", atr->add_fld[12]);
         }
      }

      for(j=0; j < tr_count; j++) {
         trr = alltr + j;
         for(k=0; k < trr->num; k++) {
             atr = (trr->trpt) + k;
             free(atr->add_fld);
         }
         free(trr->trpt);
      }

   }

   fclose(fcsvw);

   return 0;
}
