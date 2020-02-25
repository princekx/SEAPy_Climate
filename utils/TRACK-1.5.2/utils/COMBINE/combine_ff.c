#include <stdio.h>
#include <stdlib.h>
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"

/* combine sets of tr_trs files */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);

int meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );

int noheader=1;

extern int aniso;
extern int nff, nfld;
extern int *nfwpos;

extern float sum_wt, sum_per;
extern int iper_num;


int main(void)

{


    int i, j, trc=0;
    int nf=0;
    int itsc=0, itim=1;
    int ipert=0;
    int awt=0;
    
    int nffo=0, nfldo=0;

    long int pl=0, plw=0;

    int gpr, ipr;
    float alat, alng;
    float sumper=0.0, sumwt=0.0;

    int trnum=0;
    int tr_count;

    char filnamin[500];
    char filnamout[]="combined_tr_trs";

    FILE *fin=NULL;
    FILE *fout=NULL;


    struct tot_tr *alltr=NULL;

    fout = fopen(filnamout, "w");
    if(!fout){
      printf("****ERROR****, can't open file %s\n", filnamout);
      exit(1);
    }


    printf("How many files do you want to combine?\n");
    scanf("%d", &nf);


    printf("What time scaling is required?   \r\n"
           " '1'    for seasonal.            \r\n"
           " '2'    for monthly.             \r\n"
           " '3'    for daily.               \n\n");
    scanf("%d", &itsc);

    if(itsc == 2){
       printf("How many months per season?\n");
       scanf("%d", &itim);

    }
    else if(itsc == 3){
       printf("How many days per season?\n");
       scanf("%d", &itim);
    }
    

    while(trc < nf){

       printf("What is the next tr_trs file?\n");
       if(scanf("%s", filnamin) < 0){
          printf("****WARNING****, insufficient files for chosen number.\n\n");
          break;
       }
       printf("%s\n", filnamin);

       fin=fopen(filnamin, "r");
       if(!fin){
          printf("****ERROR****, can't open file %s\n", filnamin);
          exit(1);
       }

       alltr = read_tracks(fin, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
       
       if(trc){
         if(nff != nffo || nfld != nfldo){
	    printf("****ERROR****, number of additional fields in track files differs, exiting.\n\n");
	    exit(1);
	 }
       
       }
       
       nffo = nff;
       nfldo = nfld;

       if(!trc){

          if(aniso == 'y') fprintf(fout, "%d\n", 1);
          else fprintf(fout, "%d\n", 0);

          plw = ftell(fout);
	  
	  if(alltr->awt) {
             awt = alltr->awt;
             fprintf(fout, "%-8s %1d %12.5f\n", "WT_INFO", awt, sumwt);
	     ipert = 1;
	  }  

          else {
	     fprintf(fout, "%-8s %1d %12.5f\n", "PER_INFO", 1, sumper);
	  }

          fprintf(fout, "%d %d\n", gpr, ipr);

          if(gpr) fprintf(fout, "%f %f\n", alat, alng);

          pl = ftell(fout);

          fprintf(fout, "TRACK_NUM  %8d ADD_FLD  %3d %3d &", tr_count, nff, nfld);
          if(nfwpos){
             for(i=0; i < nff; i++)fprintf(fout, "%1d", *(nfwpos + i));
          }
          fprintf(fout, "\n");


       }

       trnum += meantrd(fout, alltr , tr_count, 's', gpr, ipr, alat, alng);

       for(i=0; i < tr_count; i++) {
           if(nff){
              for(j=0; j < (alltr+i)->num; j++) free(((alltr+i)->trpt + j)->add_fld);
           }
           free((alltr+i)->trpt);
       }

       free(alltr);

       ++trc;

       fclose(fin);
       fin=NULL;

       sumwt += sum_wt;

       if(iper_num){
	  sumper += sum_per;
       }
       else {
          sumper += (float)itim;
       }

    } 

    fseek(fout, pl, FSTART);
    fprintf(fout, "TRACK_NUM  %8d", trnum);
    fseek(fout, plw, FSTART);
       
    if(!ipert) {
       fprintf(fout, "%-8s %1d %12.5f\n", "PER_INFO", 1, sumper);
    }
    else {
       fprintf(fout, "%-8s %1d %12.5f\n", "WT_INFO", awt, sumwt);
    }

    fclose(fout);



    return 0;

}
