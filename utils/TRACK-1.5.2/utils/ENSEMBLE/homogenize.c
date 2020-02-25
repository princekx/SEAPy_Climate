#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "st_fo.h"
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"

/* Program to replace analysis track */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
int toverlap(struct tot_tr * , struct tot_tr * , long int * , long int * );
void truncate(struct tot_tr * , struct tot_tr * , struct tot_tr * );

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int main(int argc, char **argv)
{
    int ftype='s';

    int trnum1=0, trnum2=0, trnum3=0;
    int gpr1, ipr1, gpr2, ipr2, gpr3, ipr3;

    int nff1=0, nff2=0, nff3=0;
    int nfld1=0, nfld2=0, nfld3=0;
    int *nfwpos1=NULL, *nfwpos2=NULL, *nfwpos3=NULL;

    long int tmin=0, tmax=0;

    float alat1=0.0, alng1=0.0, alat2=0.0, alng2=0.0, alat3=0.0, alng3=0.0;

    FILE *fin1=NULL, *fin2=NULL, *fin3=NULL;
    FILE *fout=NULL;

    char fname[50];
    char *ct=NULL;

    struct tot_tr *tr1=NULL, *tr2=NULL, *tr3=NULL;
    struct tot_tr *tr11=NULL, *tr33=NULL;
    struct tot_tr *trmin=NULL, *trmax=NULL;

    if(argc < 4){
       printf("Usage: replace [ensemble track file1] [mean track file2] [detm. track file3] \r\n");
       exit(1);
    }
    
    fin1 = fopen(argv[1], "r");
    if(!fin1){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[1]);
       exit(1);
    }
    
    fin2 = fopen(argv[2], "r");
    if(!fin2){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[2]);
       exit(1);
    }

    fin3 = fopen(argv[3], "r");
    if(!fin2){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[3]);
       exit(1);
    }

    tr1 = read_tracks(fin1, &trnum1, &gpr1, &ipr1, ftype, &alat1, &alng1, NULL, NULL);
    nff1 = nff;
    nfld1 = nfld;
    nfwpos1 = nfwpos;
    nfwpos = NULL;

    tr2 = read_tracks(fin2, &trnum2, &gpr2, &ipr2, ftype, &alat2, &alng2, NULL, NULL);
    nff2 = nff;
    nfld2 = nfld;
    nfwpos2 = nfwpos;
    nfwpos = NULL;

    tr3 = read_tracks(fin3, &trnum3, &gpr3, &ipr3, ftype, &alat3, &alng3, NULL, NULL);
    nff3 = nff;
    nfld3 = nfld;
    nfwpos3 = nfwpos;
    nfwpos = NULL;
    
    fclose(fin1);
    fclose(fin2);
    fclose(fin3);

    tr11 = tr1 + 1;
    tr33 = tr3 + 1; 

    if(tr11->trid != 1){
       printf("****ERROR****, not control forecast.\n\n");
       exit(1);
    }

    tmax = tr11->trpt->time;
    tmin = (tr11->trpt + tr11->num - 1)->time; 
    trmin = tr11;
    trmax = tr11;

    if(tr2->trpt->time > tmax){tmax = tr2->trpt->time; trmax = tr2;}
    if((tr2->trpt + tr2->num - 1)->time < tmin){tmin = (tr2->trpt + tr2->num - 1)->time; trmin = tr2;}
    if(tr33->trpt->time > tmax){tmax = tr33->trpt->time; trmax = tr33;}
    if((tr33->trpt + tr33->num - 1)->time < tmin){tmin = (tr33->trpt + tr33->num - 1)->time; trmin = tr33;}

    truncate(tr11 , trmin , trmax);
    truncate(tr2 , trmin , trmax);
    truncate(tr33 , trmin , trmax); 

    ct = strstr(argv[1], "trmatch");
    strcpy(fname, ct);
       
    fout = fopen(fname, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", fname);
       exit(1);
    }  
       
    nff = nff1;
    nfld = nfld1;
    nfwpos = nfwpos1;
       
    meantrd(fout, tr1, trnum1, ftype, gpr1, ipr1, alat1, alng1);

    fclose(fout);

    ct = strstr(argv[2], "trmatch");
    strcpy(fname, ct);

    fout = fopen(fname, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", fname);
       exit(1);
    }

    nff = nff2;
    nfld = nfld2;
    nfwpos = nfwpos2;

    meantrd(fout, tr2, trnum2, ftype, gpr2, ipr2, alat2, alng2);

    fclose(fout);

    ct = strstr(argv[3], "trmatch");
    strcpy(fname, ct);
    strcat(fname, "_det");

    fout = fopen(fname, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", fname);
       exit(1);
    }

    nff = nff3;
    nfld = nfld3;
    nfwpos = nfwpos3;

    meantrd(fout, tr3, trnum3, ftype, gpr3, ipr3, alat3, alng3);

    fclose(fout);


    return 0;

}


void truncate(struct tot_tr *trr, struct tot_tr *trmin, struct tot_tr *trmax)
{
    int i=0;

    int it1s=0, it1e=0, it2s=0, it2e=0;

    long int is1=0, is2=0;

    struct fet_pt_tr *fp1=NULL, *fp2=NULL;

    if(toverlap(trr, trmax, &is1, &is2)){

       for(i=0; i < trr->num; i++){

           fp1 = trr->trpt + i;

           if(fp1->time == is1) it1s = i;
           if(fp1->time == is2) it1e = i;

       }

       for(i=0; i < trmax->num; i++){

           fp2 = trmax->trpt + i;

           if(fp2->time == is1) it2s = i;
           if(fp2->time == is2) it2e = i;

       }

       trr->num = it1e - it1s + 1;
       trr->trpt = trr->trpt + it1s;
       trr->time = trr->trpt->time;

    }
    else {
       printf("****ERROR****, no overlap between tracks.\n\n");
       exit(1);
    }

    if(toverlap(trr, trmin, &is1, &is2)){

       for(i=0; i < trr->num; i++){

           fp1 = trr->trpt + i;

           if(fp1->time == is1) it1s = i;
           if(fp1->time == is2) it1e = i;

       }

       for(i=0; i < trmin->num; i++){

           fp2 = trmin->trpt + i;

           if(fp2->time == is1) it2s = i;
           if(fp2->time == is2) it2e = i;

       }

       trr->num = it1e - it1s + 1;
       trr->trpt = trr->trpt + it1s;
       trr->time = trr->trpt->time;

    }
    else {
       printf("****ERROR****, no overlap between tracks.\n\n");
       exit(1);
    }


    return;
}
