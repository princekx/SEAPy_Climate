#include <stdio.h>
#include <stdlib.h>
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"

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
    int ilon=0, ilat=0;
    int iflon=0, iflat=0;

    float alat, alng;

    FILE *fin=NULL;
    FILE *fout=NULL;

    char filnamin[100];
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

    printf("What is the added longitude and latitude ids?            \r\n"
           "Note, this should be the last added fields at the moment.\n\n");
    scanf("%d %d", &ilon, &ilat);

    iflon = 0;
    for(i=0; i < ilon - 1; i++){
         if(*(nfwpos + i)) iflon += 3;
         else iflon += 1;
    }

    iflat = 0;
    for(i=0; i < ilat - 1; i++){
         if(*(nfwpos + i)) iflat += 3;
         else iflat += 1;
    }

    for(i=0; i < trnum; i++){
        trr = alltr + i;
        for(j=0; j < trr->num; j++){
            fp = trr->trpt + j;
            if(*(fp->add_fld + iflon) < 0.0) *(fp->add_fld + iflon) = 360.0 + *(fp->add_fld + iflon);
            fp->xf = *(fp->add_fld + iflon);
            fp->yf = *(fp->add_fld + iflat);
            if(fp->xf < 0.0) fp->xf = 360.0 + fp->xf;
        }
    }

    nff =nff-2; nfld=nfld-2;

    fout=fopen(filnamout, "w");
    if(!fin){
       printf("****ERROR****, can't open file %s\n", filnamout);
       exit(1);
    }

     meantrd(fout, alltr , trnum, 's', gpr, ipr, alat, alng);

    fclose(fin);    

    return 0;
}
