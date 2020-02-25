#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "splice.h"

#define MAXCHR     50

/* program to convolve the weights from two different teleconnections */

void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float *, float ** , int *);

int tom='g';
int noheader=0;

extern int aniso;

extern float sum_per, sum_wt;
extern int iper_num;

int main(int argc, char **argv)
{

    int i, j;
    int trnum1, trnum2;
    int gpr1, ipr1, gpr2, ipr2;
    int nw1, nw2;

    float *wght1, *wght2;
    float swght=0.0;

    char trout[MAXCHR];

    char usg[] = "Usage: convolve [track file 1] [track file 2] [outfile]";

    FILE *fin1=NULL, *fin2=NULL;
    FILE *fout=NULL;

    float alat1=0.0, alng1=0.0, alat2=0.0, alng2=0.0;

    struct tot_tr *tr1=NULL, *tr2=NULL;
    struct tot_tr *tm1=NULL, *tm2=NULL;
    struct fet_pt_tr *at1=NULL, *at2=NULL;

    if(argc != 4){

      printf("%s\n", usg);
      exit(1);

    }

    printf("***WRANING***, it is the users responsibility that the two\r\n"
           "               track files are identical apart from their \r\n"
           "               weight values.                             \n\n");

    sscanf(trout, "%s", argv[3]);

    fin1 = fopen(argv[1], "r");
    if(!fin1){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[1]);
       exit(1);
    }

    fin2 = fopen(argv[2], "r");
    if(!fin1){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[2]);
       exit(1);
    }

    tr1 = read_tracks(fin1, &trnum1, &gpr1, &ipr1, 's', &alat1, &alng1, &wght1, &nw1);
    tr2 = read_tracks(fin2, &trnum2, &gpr2, &ipr2, 's', &alat2, &alng2, &wght2, &nw2);

    fclose(fin1);
    fclose(fin2);

    for(i=0; i < nw1; i++){
        swght += *(wght1 + i) * *(wght2 + i);
    }


    if(trnum1 != trnum2 || gpr1 != gpr2 || ipr1 != ipr2 ||
       (fabs(alat1 - alat2) > 1.0e-4 || fabs(alng1 - alng2) > 1.0e-4) ||
       nw1 != nw2){

       printf("****ERROR****, possibly two different sets of tracks.\n\n");
       exit(1);

    }

    if(!(tr1->awt) || !(tr2->awt)){
       printf("****ERROR****, one or both sets of tracks do not contain weight information.\n\n");
       exit(1);

    } 

    for(i=0; i < trnum1; i++){

        tm1 = tr1 + i;
        tm2 = tr2 + i;

        if(tm1->num != tm2->num){
           printf("****ERROR****, tracks differ in number of points\n\n");
           exit(1);
        }

        for(j=0; j < tm1->num; j++){

            at1 = tm1->trpt + j;
            at2 = tm2->trpt + j;

            at1->wght *= at2->wght;

        }


    }

    fout = fopen(argv[3], "w");
    if(!fin1){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", argv[3]);
       exit(1);
    }

    sum_wt = swght;

    meantrd(fout, tr1, trnum1, 's', gpr1, ipr1, alat1, alng1);

    fclose(fout);

    return 0;
}
