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

extern int iper_num;

extern float sum_wt, sum_per;

int main(int argc, char **argv)
{

   int i, j, k, l;
   
   int tr_count=0;
   int gpr, ipr;
   int imd=0;
   int iememb=0;

   int itrid=0;
 

   float alat, alng;
   
/*   long int datei, datef; */

   char trfil[LENSTR];
   char fcsv[]="tr_wisc";



   FILE *fcsvw=NULL;
   FILE *ftr=NULL;

   struct tot_tr *alltr=NULL, *trr=NULL;
   struct fet_pt_tr *atr;

   if(argc != 3){
      printf("Usage: tr2csv [trfile] [ensemble memb]\n\n");
      exit(1);
   } 

   sscanf(argv[1], "%s", trfil); 
   sscanf(argv[2], "%d", &iememb);

   ftr = fopen(trfil, "r");
   if(!ftr){
      printf("****ERROR****, can't open file %s for 'r'\n\n", trfil);
      exit(1);
   }
   alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
   fclose(ftr);

   fcsvw = fopen(fcsv, "w");
   
   if(aniso == 'y') fprintf(fcsvw, "%d\n", 1);
   else fprintf(fcsvw, "%d\n", 0);
   
   if(iper_num){
      fprintf(fcsvw, "PER_INFO %1d %12.5f\n", 1, sum_per);
   }
   
   fprintf(fcsvw, "%d %d\n", gpr, ipr);
   
   fprintf(fcsvw, "TRACK_NUM  %8d ADD_FLD  %3d %3d &", tr_count, nff, nfld);
   if(nfwpos){
      for(i=0; i < nff; i++)fprintf(fcsvw, "%1d", *(nfwpos + i));
   }
   fprintf(fcsvw, "\n");

   for(j=0; j < tr_count; j++){
       ++itrid;
       trr = alltr + j;
       
       imd = (trr->num / 2);
       
       if(trr->time) fprintf(fcsvw, "TRACK_ID  %d START_TIME %10ld MEDIAN_TIME %10ld ENSMBLE_MEMB %d NO-NAME\n", trr->trid, trr->time, ((trr->trpt) + imd)->time, iememb);
       else fprintf(fcsvw, "TRACK_ID  %d\n", trr->trid);
       
       fprintf(fcsvw, "POINT_NUM  %d\n", trr->num);

       for(k=0; k < trr->num; k++){
           atr = (trr->trpt) + k;
/* reformat date */
           if(atr->time) {

              fprintf(fcsvw, "%ld %f %f %e ", atr->time, atr->xf, atr->yf, atr->zf);

           }
           else {
              fprintf(fcsvw, "%d %f %f %e ",  atr->fr_id, atr->xf, atr->yf, atr->zf);
           }
	   if(aniso == 'y'){
	      fprintf(fcsvw, "%f %f %f %e ",  atr->sh_an, atr->or_vec[0], atr->or_vec[1], atr->area);
	   }
	   
	   if(nfld){
	      fprintf(fcsvw, "& ");
              for(l=0; l < nfld; l++) fprintf(fcsvw, "%e & ", atr->add_fld[l]);
	   }
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




