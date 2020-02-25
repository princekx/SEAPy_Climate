#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem_er.h"
#include "file_handle.h"

#define  MAXCHR  1000

int main(void )
{
    int i, j;
    int nsamp=0, trnum=0;
    int isamp=0, nnsmp=0;
    int ns=0, nl=0;
    int *imiss=NULL;
    int *nb=NULL;
    int ntot=0;

    long int *cfp=NULL, *cfr=NULL;
    long int fst;

    char mfnm[]="missing.dat";
    char cmfnm[]="conv_missing.dat";
    char line[MAXCHR], *ll=NULL;

    FILE *fin=NULL, *fout=NULL;
    
    fin = fopen(mfnm, "r");
    if(!fin){
       printf("****ERROR****, unable to open file %s for read.\n\n", mfnm);
       exit(1);
    }

    fout = fopen(cmfnm, "w");
    if(!fin){
       printf("****ERROR****, unable to open file %s for write.\n\n", cmfnm);
       exit(1);
    }
    
    fgets(line, MAXCHR, fin);
    sscanf(line, "%d %d", &nsamp, &trnum);

    cfp = (long int *)calloc(trnum, sizeof(long int));
    mem_er((cfp == NULL) ? 0 : 1, trnum*sizeof(long int));

    cfr = (long int *)calloc(trnum, sizeof(long int));
    mem_er((cfr == NULL) ? 0 : 1, trnum*sizeof(long int));

    fst = ftell(fout);

    fprintf(fout, "%6d %10d\n", trnum, ntot); 

    nb = (int *)calloc(trnum, sizeof(int));
    mem_er((nb == NULL) ? 0 : 1, trnum*sizeof(int));

    for(i=0; i < trnum; i++){
        *(cfr + i) = ftell(fout);
        fprintf(fout, "%6d ", 0);
        *(cfp + i) = ftell(fout);
        for(j=0; j < nsamp; j++) fprintf(fout, "      ");
        fprintf(fout, "\n");
    }

    for(i=0; i < nsamp; i++){
        printf("Sample %d\n", i);
        ns = 0;
        fgets(line, MAXCHR, fin);
        sscanf(line, "%d %d", &isamp, &nnsmp);
        if(!imiss){
           imiss = (int *)calloc(nnsmp, sizeof(int));
           mem_er((imiss == NULL) ? 0 : 1, nnsmp*sizeof(int));
        }
        else {
           imiss = (int *)realloc_n(imiss, nnsmp * sizeof(int));
           mem_er((imiss == NULL) ? 0 : 1, nnsmp*sizeof(int));
        }       
        while(ns < nnsmp){
              nl = 0;
              fgets(line, MAXCHR, fin);
              ll = strtok(line, " ");
              sscanf(ll, "%d", imiss + ns);
              ++ns; ++nl;
              if(ns >= nnsmp) break;
              while(nl < 10){
                ll = strtok(NULL, " ");
                sscanf(ll, "%d", imiss + ns);
                ++ns; ++nl;
                if(ns >= nnsmp) break;
              }
              
        }

        for(j=0; j < nnsmp; j++) {
            fseek(fout, *(cfp + *(imiss + j)), FSTART);
            ++(*(nb + *(imiss + j)));
            fprintf(fout, "%d ", i);
            *(cfp + *(imiss + j)) = ftell(fout);
        }
    }

    for(i=0; i < trnum; i++){
        fseek(fout, *(cfr + i), FSTART);
        fprintf(fout, "%6d", *(nb + i));
        ntot += *(nb + i);
    }

    fseek(fout, fst, FSTART);
    fprintf(fout, "%6d %10d\n", trnum, ntot);

    fclose(fin);
    fclose(fout);

    free(imiss);
    free(cfp);
    free(cfr);
    free(nb);

    return 0;

}
