#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "splice.h"
#include "mem_er.h"

#define LENSTR 100

/* program to convert track files into comma seperated columns */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);

int aniso;
int nff, nfld;
int *nfwpos;

int main(int argc, char **argv)
{

   int i, j, k, l;
   
   int tr_count=0;
   int gpr, ipr;

   int itrid=0;
  
   int yr, mn, dy, hr;

   float alat, alng;
   
/*   long int datei, datef; */

   char trfil[LENSTR];
   char fcsv[]="alltr.csv";
/*   char ffld[5]; */
   char date[20];

   FILE *fcsvw=NULL;
   FILE *ftr=NULL;

   struct tot_tr *alltr=NULL, *trr=NULL;
   struct fet_pt_tr *atr;

   if(argc != 2){
      printf("Usage: tr2csv [trfile]\n\n");
      exit(1);
   } 

   sscanf(argv[1], "%s", trfil); 

   ftr = fopen(trfil, "r");
   if(!ftr){
      printf("****ERROR****, can't open file %s for 'r'\n\n", trfil);
      exit(1);
   }
   alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
   fclose(ftr);

   fcsvw = fopen(fcsv, "w");

   for(j=0; j < tr_count; j++){
       ++itrid;
       trr = alltr + j;

       for(k=0; k < trr->num; k++){
           atr = (trr->trpt) + k;
/* reformat date */
           if(atr->time) {
              yr = atr->time / 1000000;
              mn = (atr->time - yr * 1000000) / 10000;
              dy = (atr->time - yr * 1000000 - mn * 10000) / 100;
              hr = atr->time - yr * 1000000 - mn * 10000 - dy * 100;
              sprintf(date, "%4d-%02d-%02d %02d:00:00", yr, mn, dy, hr);

              fprintf(fcsvw, "%-8d, \"%18s\", %f, %f, %e, ", itrid, date, atr->xf, atr->yf, atr->zf);

           }
           else {
              fprintf(fcsvw, "%-8d, %-8d, %f, %f, %e, ", itrid, atr->fr_id, atr->xf, atr->yf, atr->zf);
           }
	   if(aniso == 'y'){
	      fprintf(fcsvw, "%f, %f, %f, %e, ",  atr->sh_an, atr->or_vec[0], atr->or_vec[1], atr->area);
	   }
           for(l=0; l < nfld; l++) fprintf(fcsvw, "%e, ", atr->add_fld[l]);
           fprintf(fcsvw, "\n");
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


   fclose(fcsvw);

   return 0;
}




