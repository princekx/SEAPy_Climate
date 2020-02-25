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
    int i=0, j=0;
    
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
    struct tot_tr *trr1=NULL, *trr2=NULL;
    struct fet_pt_tr *fp1=NULL, *fp2=NULL;

    if(argc < 4){
       printf("Usage: trunc [track file1] [track file2] [outfile] \r\n");
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
    nfwpos1 = (int *)calloc(nff1, sizeof(int));
    mem_er((nfwpos1 == NULL) ? 0 : 1, nff1 * sizeof(int));
    memcpy(nfwpos1, nfwpos, nff1*sizeof(int));

    tr2 = read_tracks(fin2, &trnum2, &gpr2, &ipr2, ftype, &alat2, &alng2, NULL, NULL);
    nff2 = nff;
    nfld2 = nfld;
    nfwpos2 = (int *)calloc(nff2, sizeof(int));
    mem_er((nfwpos2 == NULL) ? 0 : 1, nff1 * sizeof(int));
    memcpy(nfwpos2, nfwpos, nff2*sizeof(int));
    
    fclose(fin1);
    fclose(fin2);
    
    if(trnum1 != trnum2){
       printf("****ERROR****, different numbers of track in the two files.\n\n");
       exit(1);
    }
    
    for(i=0; i < trnum1; i++){
    
       trr1 = tr1 + i;
       trr2 = tr2 + i;

       if(toverlap(trr1, trr2, &is1, &is2)){

          for(j=0; j < trr1->num; j++){

              fp1 = trr1->trpt + j;

              if(fp1->time == is1) it1s = j;
              if(fp1->time == is2) it1e = j;

          }

          for(j=0; j < trr2->num; j++){

              fp2 = trr2->trpt + j;

              if(fp2->time == is1) it2s = j;
              if(fp2->time == is2) it2e = j;

          }

          if(it1e - it1s != it2e - it2s) {
             printf("****WARNING****, truncated track links differ.\n\n");
             exit(1);
          }

          trr1->num = it1e - it1s + 1;
          trr1->trpt = trr1->trpt + it1s;
          trr1->time = trr1->trpt->time;

       }
       else {
          printf("****WARNING****. no overlap from track %d\n\n", i+1);
       } 
   } 
       
   fout = fopen(argv[3], "w");
   if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", argv[1]);
       exit(1);
   }  
       
   nff = nff1;
   nfld = nfld1;
   free(nfwpos);
   nfwpos = (int *)calloc(nff, sizeof(int));
   mem_er((nfwpos == NULL) ? 0 : 1, nff * sizeof(int));
   memcpy(nfwpos, nfwpos1, nff1*sizeof(int));
       
   meantrd(fout, tr1, trnum1, ftype, gpr1, ipr1, alat1, alng1); 

   fclose(fout);

   return 0;

}
