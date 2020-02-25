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

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int main(int argc, char **argv)
{
    int i;
    
    int ftype='s';

    int trnum1, trnum2;
    int gpr1, ipr1, gpr2, ipr2;

    int it1s=0, it1e=0, it2s=0, it2e=0;

    int nff1=0, nff2=0;
    int nfld1=0, nfld2=0;
    int *nfwpos1=NULL, *nfwpos2=NULL;

    long int is1=0, is2=0;

    float alat1=0.0, alng1=0.0, alat2=0.0, alng2=0.0;

    FILE *fin1=NULL, *fin2=NULL;
    FILE *fout=NULL;

    struct tot_tr *tr1=NULL, *tr2=NULL;
    struct fet_pt_tr *fp1=NULL, *fp2=NULL;

    if(argc < 4){
       printf("Usage: replace [track file1] [track file2] [outfile] \r\n");
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
    
    fclose(fin1);
    fclose(fin2);

    if(toverlap(tr1, tr2, &is1, &is2)){

       for(i=0; i < tr1->num; i++){

           fp1 = tr1->trpt + i;

           if(fp1->time == is1) it1s = i;
           if(fp1->time == is2) it1e = i;

       }

       for(i=0; i < tr2->num; i++){

           fp2 = tr2->trpt + i;

           if(fp2->time == is1) it2s = i;
           if(fp2->time == is2) it2e = i;

       }

       tr1->num = it1e - it1s + 1;
       tr1->trpt = tr1->trpt + it1s;
       tr1->time = tr1->trpt->time;
       
       tr2->num = it2e - it2s + 1;
       tr2->trpt = tr2->trpt + it2s;  
       
       if(tr1->num != tr2->num){
          printf("****ERROR****, track lengths for copy do not match.\n\n");
	  exit(1);
       }
       
       for(i=0; i < tr1->num; i++){
           fp1 = tr1->trpt + i;
	   fp2 = tr2->trpt + i;
	   
	   fp1->add_fld[7] = fp2->xf;
	   fp1->add_fld[8] = fp2->yf;
	   if(fp2->add_fld[0] > ADD_CHECK) fp1->add_fld[9] = ADD_UNDEF;
	   else if(fp2->add_fld[0] <= 0.0) fp1->add_fld[9] = ADD_UNDEF;
	   else fp1->add_fld[9] = fp2->add_fld[0];
       
       }
       
       fout = fopen(argv[3], "w");
       if(!fout){
          printf("***ERROR***, unable to open file %s for 'w'\n\n", argv[1]);
          exit(1);
       }  
       
       nff = nff1;
       nfld = nfld1;
       nfwpos = nfwpos1;
       
       meantrd(fout, tr1, trnum1, ftype, gpr1, ipr1, alat1, alng1);

    }

    else {
       printf("****ERROR****, no overlap between tracks.\n\n");
       exit(1);
    }
    


    return 0;

}
