#include <stdio.h>
#include <stdlib.h>
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"

#define  PERIOD  360.0

/* combine sets of tr_trs files */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);

int meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );

int noheader=0;

extern int aniso;
extern int aniso;
extern int nff, nfld;
extern int *nfwpos;

int main(void)

{
    int i, j;
    int gpr, ipr;
    int trnum=0;

    float alat, alng;
    float tlng=0.0;
    float lngtr=0.0;

    FILE *fin=NULL;
    FILE *fout=NULL;

    char filnamin[200];
    char filnamout[]="tr_trans.dat";

    struct tot_tr *alltr=NULL, *trr=NULL;
    struct fet_pt_tr *fp=NULL;

    printf("What is the track file to read?\n\n");
    scanf("%s", filnamin);

    fin=fopen(filnamin, "r");
    if(!fin){
       printf("****ERROR****, can't open file %s\n", filnamin);
       exit(1);
    }

    alltr = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);

    printf("To translate tracks to a common starting longitude, what longitude is required?\r\n"
           "Note value should be in the range 0-360.                          \n\n");
    
    scanf("%f", &tlng);

    if(tlng < 0.0 || tlng > 360.0){
      printf("****ERROR****, incorrect translation value.\n\n");
      exit(1);
    }

    for(i=0; i < trnum; i++){
       trr = alltr + i;
       lngtr = trr->trpt->xf - tlng;
       for(j=0; j < trr->num; j++){
           fp = trr->trpt + j;
           fp->xf -= lngtr;
           if(fp->xf < 0.0) fp->xf += PERIOD;
           if(fp->xf > PERIOD) fp->xf -= PERIOD;
       }
    }

    fout=fopen(filnamout, "w");
    if(!fin){
       printf("****ERROR****, can't open file %s\n", filnamout);
       exit(1);
    }

    meantrd(fout, alltr , trnum, 's', gpr, ipr, alat, alng);

    fclose(fin);    

    return 0;
}
