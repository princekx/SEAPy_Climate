#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include "splice.h"

#define LENSTR 100

/* program to convert tracks for plotting in xmgrace */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);

int main(int argc, char **argv)
{
   int i, j, k;
   int hemi=0;
   int nfil=0;
   int ifr=0;
   int iclstr=0;
   int gpr, ipr;
   int gprt, iprt;
   int tr_count=0;
   int yy, mm, dd, hh;

   long int datei;

   float alat, alng;
   float alatt, alngt;

   char ffld[5];
   char lsstr[LENSTR], tmpdir[LENSTR];
   char trfil[LENSTR], trfil2[LENSTR];
   char outfil[LENSTR];

   FILE *fout=NULL;
   FILE *ftr=NULL;

   struct dirent *entry;

   DIR *dir=NULL;

   struct tot_tr *alltr=NULL, *mtr=NULL, *trr=NULL;
   struct fet_pt_tr *atr; 

   if(argc != 4){
      printf("Usage: graceconv [Date] [Field: vor or mslp] [Hemisphere: NH=0; SH=1]\n\n");
      exit(1);
   }

   sscanf(argv[1], "%ld", &datei);
   sscanf(argv[2], "%s", ffld);
   sscanf(argv[3], "%d", &hemi);

   sprintf(lsstr, "%ld_%s", datei, ffld);
   sprintf(tmpdir, "temp_%s", lsstr);

   dir = opendir(lsstr);
   if(dir == NULL){
      printf("****ERROR****, cannot open directory %s\n\n", lsstr);
      exit(1);
   }

   entry = readdir( dir );
   while(entry != NULL){
     if( strncmp( entry->d_name, "trmatch", 7 ) == 0 && ! (strstr(entry->d_name, "mean"))) ++nfil;
     entry = readdir( dir );
   }
   closedir(dir);

   printf("%d \n", nfil);

   if(chdir(lsstr)){
     printf("****ERROR****, can't change directory to %s\n", lsstr);
     exit(1);
   }

   mkdir(tmpdir, S_IRWXU);

   for(i=0; i < nfil; i++){
      iclstr = i + 1;
      if(i+1 < 10) {
         sprintf(trfil, "trmatch_cntl_tr000%d_mean", i+1);
         sprintf(trfil2, "trmatch_cntl_tr000%d", i+1);
      }
      else if(i+1 < 100) {
         sprintf(trfil, "trmatch_cntl_tr00%d_mean", i+1);
         sprintf(trfil2, "trmatch_cntl_tr00%d", i+1);
      }
      else if(i+1 < 1000) {
         sprintf(trfil, "trmatch_cntl_tr0%d_mean", i+1);
         sprintf(trfil2, "trmatch_cntl_tr0%d", i+1);
      }
      else {
         printf("****ERROR****, not configured for this many clusters.\n\n");
         exit(1);
      }

      ftr = fopen(trfil2, "r");
      if(!ftr){
         printf("****ERROR****, file %s does not exist.\n\n", trfil2);
         exit(1);
      }

      alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
      fclose(ftr);

      if(!ifr){gprt = gpr; iprt = ipr; alatt = alat; alngt = alng; ifr=1;}
      else{
         if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
            printf("****ERROR****, track files are incompatable.\n\n");
            exit(1);
         }
      }
      
      for(j=0; j < tr_count; j++){
         trr = alltr + j;
         if(!(trr->trid)){
            sprintf(outfil, "%s/cntl_%d.dat", tmpdir, i+1);
         }
         else if (trr->trid == 1){
            sprintf(outfil, "%s/det_%d.dat", tmpdir, i+1);
         }
         else {
            sprintf(outfil, "%s/eps_%d_%d.dat", tmpdir, i+1, trr->trid - 1);
         }
         fout = fopen(outfil, "w");
         if(!fout){
            printf("****ERROR****, file %s cannot be opened.\n\n", outfil);
            exit(1);
         }

         for(k=0; k < trr->num; k++){
            atr = trr->trpt + k;            
            yy = atr->time / 1000000;
            mm = (atr->time - yy * 1000000) / 10000;
            dd = (atr->time - yy * 1000000 - mm * 10000) / 100;
            hh = atr->time - yy * 1000000 - mm * 10000 - dd * 100;
            fprintf(fout, "%d-%d-%d-%d %e %e %e %e %e %e\n", yy, mm, dd, hh, atr->zf, atr->add_fld[2], atr->add_fld[5], atr->add_fld[8], atr->add_fld[11], atr->add_fld[12]); 

         } 

         fclose(fout); 

      }

      ftr = fopen(trfil, "r");
      if(!ftr){
         printf("****ERROR****, file %s does not exist.\n\n", trfil);
      }
      else{
         alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
         fclose(ftr);

         if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
            printf("****ERROR****, track files are incompatable.\n\n");
            exit(1);
         }
      
         sprintf(outfil, "%s/mean_%d.dat", tmpdir, i+1);

         fout = fopen(outfil, "w");
         if(!fout){
           printf("****ERROR****, file %s cannot be opened.\n\n", outfil);
           exit(1);
         }

         for(k=0; k < alltr->num; k++){
            atr = alltr->trpt + k;
            yy = atr->time / 1000000;
            mm = (atr->time - yy * 1000000) / 10000;
            dd = (atr->time - yy * 1000000 - mm * 10000) / 100;
            hh = atr->time - yy * 1000000 - mm * 10000 - dd * 100;
            fprintf(fout, "%d-%d-%d-%d %e %e %e %e %e %e\n", yy, mm, dd, hh, atr->zf, atr->add_fld[2], atr->add_fld[5], atr->add_fld[8], atr->add_fld[11], atr->add_fld[12]);

         }

         fclose(fout);
      }

   }


   return 0;
}
