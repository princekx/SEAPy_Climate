#include <Stdio.h>
#include <math.h>
#include <stdlib.h>
#include "splice.h"
#include "mem_er.h"

#define MAXCHR  200

/* vertically integrate chosen additional fields */

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
void sincos(double , double * , double * );

int main(void)
{

    int i, j, k;
    int trnum=0;
    int gpr=0, ipr=0;
    int nlev;

    int *ilev=NULL, *illev=NULL;

    FILE *ftr=NULL, *fout=NULL;

    char filin[MAXCHR], filout[MAXCHR];

    float *flev=NULL, *flint=NULL;
    float alat, alng;
    float sum=0.0;
    float scl=0.0;

    struct tot_tr *tracks=NULL, *trr=NULL;
    struct fet_pt_tr *fp=NULL;

    printf("What is the track file to read?\n\n");
    scanf("%s", filin);

    ftr = fopen(filin, "r");
    if(!ftr){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
       exit(1);
    }

    tracks = read_tracks(ftr, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(ftr);

    if(!nff){
       printf("****ERROR****, no additional fields available.\n\n");
       exit(1);
    }

    printf("How many levels are to be integrated over?\n\n");
    printf("NOTE additional fields should be stored in contiguous order.\n\n");
    scanf("%d", &nlev);

    ilev = (int *)calloc(nlev, sizeof(int));
    mem_er((ilev == NULL) ? 0 : 1, nlev*sizeof(int));
    illev = (int *)calloc(nlev, sizeof(int));
    mem_er((illev == NULL) ? 0 : 1, nlev*sizeof(int));

    flev = (float *)calloc(nlev, sizeof(float));
    mem_er((flev == NULL) ? 0 : 1, nlev*sizeof(float));

    flint = (float *)calloc(nlev, sizeof(float));
    mem_er((flint == NULL) ? 0 : 1, nlev*sizeof(float));

    for(i=0; i < nlev; i++){
        printf("What is the next level Id. and value?\n\n");
        scanf("%d %f", ilev + i, flev + i);
        *(illev + i) = 0;
        for(j=0; j < *(ilev + i); j++){
            if(*(nfwpos + j)) *(illev + i) += 3;
            else *(illev + i) += 1;
        }
        --(*(illev + i));

    }

    ++nff;
    ++nfld;
    nfwpos = (int *)realloc_n(nfwpos, nff * sizeof(int));
    mem_er((nfwpos == NULL) ? 0 : 1, nff * sizeof(int));
    *(nfwpos + nff - 1) = 0;


    for(i=0; i < nlev-1; i++) *(flint + i) = fabs(*(flev + i) - *(flev + i + 1));
    *(flint + nlev - 1) = *(flint + nlev - 2);

    printf("What scaling is required?\n\n");
    scanf("%f", &scl);
        
    for(i=0; i < trnum; i++){
        trr = tracks + i;

        for(j=0; j < trr->num; j++){
            fp = trr->trpt + j;
            fp->add_fld = (float *)realloc_n(fp->add_fld, nfld * sizeof(float));
            mem_er((fp->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
            sum = 0.0;
            for(k=0; k < nlev; k++){
               sum += fp->add_fld[*(illev + k)] * *(flint + k);
            }
            fp->add_fld[nfld - 1] = sum * scl;
        }

    }

    printf("What is the output file required?\n\n");
    scanf("%s", filout);

    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
       exit(1);
    }

    meantrd(fout, tracks, trnum, 's', gpr, ipr, alat, alng);

    fclose(fout);

    free(ilev);
    free(flev);

    return 0;

}
