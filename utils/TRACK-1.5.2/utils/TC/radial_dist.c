#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "st_fo.h"
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"

/* Program to compute radial distance between additional fields */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
void sincos(double , double * , double * );

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

    int trnum;
    int gpr=0, ipr=0;
    int infl1=0, infl2=0;
    int iff1=0, iff2=0;
    int ifpos1=0, ifpos2=0;

    float alat=0.0, alng=0.0;
    float xx=0.0, yy=0.0;

    double x1, y1, z1;
    double x2, y2, z2;
    double s1=0.0, s2=0.0, c1=0.0, c2=0.0;

    FILE *fin=NULL;
    FILE *fout=NULL;

    struct tot_tr *trr=NULL, *altr=NULL;
    struct fet_pt_tr *atr=NULL;

    if(argc < 5){
       printf("Usage: trunc [track file] [outfile] [1st add field] [2nd add field]\r\n");
       exit(1);
    }


    sscanf(argv[3], "%d", &infl1);
    sscanf(argv[4], "%d", &infl2);
    
    fin = fopen(argv[1], "r");
    if(!fin){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[1]);
       exit(1);
    }
    
    trr = read_tracks(fin, &trnum, &gpr, &ipr, ftype, &alat, &alng, NULL, NULL);
    
    fclose(fin);

    if(infl1 < 0 || infl1 > nff) {
       printf("****ERROR***, requested default or additional field does not exist for first track file.\n\n");
       exit(1);
    }

    if(infl2 < 0 || infl2 > nff) {
       printf("****ERROR***, requested default or additional field does not exist for second track file.\n\n");
       exit(1);
    }

    if(infl1) {
      iff1 = 0;
      for(i=0; i < infl1; i++){
         if(*(nfwpos + i)) iff1 += 3;
         else iff1 += 1;
      }
      --iff1;
      if(*(nfwpos + infl1 - 1)) ifpos1 = iff1 - 2;
      else {
         printf("No location for information additional field.\n\n");
         exit(1);
      }

    }

    if(infl2) {
      iff2 = 0;
      for(i=0; i < infl2; i++){
         if(*(nfwpos + i)) iff2 += 3;
         else iff2 += 1;
      }
      --iff2;
      if(*(nfwpos + infl2 - 1)) ifpos2 = iff2 - 2;
      else {
         printf("No location for information additional field.\n\n");
         exit(1);
      }

    }

    ++nff;
    ++nfld;
    nfwpos = (int *)realloc_n(nfwpos, nff*sizeof(int));
    mem_er((nfwpos == NULL) ? 0 : 1, nff * sizeof(int));    
    *(nfwpos + nff - 1) = 0;

    for(i=0; i < trnum; i++){
       altr = trr + i;

       for(j=0; j < altr->num; j++){

           atr = altr->trpt + j;
           atr->add_fld = (float *)realloc_n(atr->add_fld, nfld * sizeof(float));
           mem_er((atr->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));


           if(!infl1){
	      xx = atr->xf * FP_PI;
	      yy = FP_PI2 - atr->yf * FP_PI;
              sincos(xx, &s1, &c1);
              sincos(yy, &s2, &c2);
           } 
           else {
	      if(atr->add_fld[ifpos1] > ADD_CHECK){
	         atr->add_fld[nfld-1] = ADD_UNDEF;
		 continue;
	      }
	      xx = atr->add_fld[ifpos1] * FP_PI;
	      yy = FP_PI2 - atr->add_fld[ifpos1+1] * FP_PI;
              sincos(xx, &s1, &c1);
              sincos(yy, &s2, &c2);
           }

           x1 = s2 * c1;
           y1 = s2 * s1;
           z1 = c2;

           if(!infl2){
	      xx = atr->xf * FP_PI;
	      yy = FP_PI2 - atr->yf * FP_PI;
              sincos(xx, &s1, &c1);
              sincos(yy, &s2, &c2);
           }
           else {
	       if(atr->add_fld[ifpos2] > ADD_CHECK){
	         atr->add_fld[nfld-1] = ADD_UNDEF;
		 continue;
	      }
	      xx = atr->add_fld[ifpos2] * FP_PI;
	      yy = FP_PI2 - atr->add_fld[ifpos2+1] * FP_PI;
              sincos(xx, &s1, &c1);
              sincos(yy, &s2, &c2);
           }

           x2 = s2 * c1;
           y2 = s2 * s1;
           z2 = c2;
           
           atr->add_fld[nfld-1] = acos(x1 * x2 + y1 * y2 + z1 * z2) / FP_PI;

       }

    }


    fout = fopen(argv[2], "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", argv[2]);
       exit(1);
    }  
       
    meantrd(fout, trr, trnum, ftype, gpr, ipr, alat, alng);

    return 0;

}
