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

    int i, j;
    int trnum=0;
    int gpr=0, ipr=0;
    int ctype=0;

    int if1=0, if2=0, iff1=0, iff2=0;

    FILE *ftr=NULL, *fout=NULL;

    char filin[MAXCHR], filout[MAXCHR];

    float alat, alng;
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


    printf("Which two additional fields do you want to combine?\n\n");
    scanf("%d %d", &if1, &if2);

    printf("How do you want the two fields combined, \r\n"
           "Addition,                   input '0'                \r\n"
           "Subtraction (2nd - 1st),    input '1'                \r\n"
           "Multiplication,             input '2'                \r\n"
           "Division (2nd/1st),         input '3'                \n\n");
    scanf("%d", &ctype);

    if(ctype < 0 || ctype > 3){
       printf("****ERROR***, operation Id. %d unknown, exiting.\n\n", ctype);
       exit(1);
    }

    for(i=0; i < if1; i++){
       if(*(nfwpos + i)) iff1 += 3;
            else iff1 += 1;
    }
    --iff1;

    for(i=0; i < if2; i++){
       if(*(nfwpos + i)) iff2 += 3;
            else iff2 += 1;
    }
    --iff2;
    
    ++nff;
    ++nfld;
    nfwpos = (int *)realloc_n(nfwpos, nff * sizeof(int));
    mem_er((nfwpos == NULL) ? 0 : 1, nff * sizeof(int));
    *(nfwpos + nff - 1) = 0;


    printf("What scaling is required?\n\n");
    scanf("%e", &scl);

    printf("%e\n", scl);
        
    for(i=0; i < trnum; i++){
        trr = tracks + i;

        for(j=0; j < trr->num; j++){
            fp = trr->trpt + j;
            fp->add_fld = (float *)realloc_n(fp->add_fld, nfld * sizeof(float));
            mem_er((fp->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
            switch (ctype) {
               case 0:
                  fp->add_fld[nfld - 1] = scl * (fp->add_fld[iff2] + fp->add_fld[iff1]);
                  break;
               case 1:
                  fp->add_fld[nfld - 1] = scl * (fp->add_fld[iff2] - fp->add_fld[iff1]);
                  break;
               case 2:
                  fp->add_fld[nfld - 1] = scl * (fp->add_fld[iff2] * fp->add_fld[iff1]);
                  break;
               case 3:
                  fp->add_fld[nfld - 1] = scl * (fp->add_fld[iff2] / fp->add_fld[iff1]);
                  break;
            }
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

    return 0;

}
