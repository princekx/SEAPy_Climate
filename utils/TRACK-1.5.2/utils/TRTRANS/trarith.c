#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"

#define  PERIOD  360.0

/* apply aritmatic to tr_trs files */

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
    int atyp=0;
    int iadf=0, iff=0;

    float alat, alng;
    float val=1.0;
    float rng1=0.0, rng2=0.0;
    float str=0.0;

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

    printf("What arithetical operation is required:    \r\n"
           "Input '0' for no operation                 \r\n"
	   "Input '1' for addition                     \r\n"
	   "Input '2' for subtraction                  \r\n"
	   "Input '3' for multiplication               \r\n"
	   "Input '4' for division                     \r\n"
           "Input '5' for absolute value               \n\n");
    scanf("%d", &atyp);
    
    printf("Which of the %d additional intensities should be modified. '0' for default.\n\n", nff);
    scanf("%d", &iadf);
    if(iadf < 0 || iadf > nff){
       printf("****ERROR****, additional field does not exist.\n\n");
       exit(1);
    }
    
    if(iadf) {
      iff = 0;
      for(i=0; i < iadf; i++){
         if(*(nfwpos + i)) iff += 3;
         else iff += 1;
      }
      --iff;
    }
    
    printf("Input a valid range, values outside of this range will be set to missing values.\n\n");
    scanf("%e %e", &rng1, &rng2);
    if(rng1 > rng2){
       printf("****ERROR****, invalid range %e %e\n\n", rng1, rng2);
       exit(1);
    }
    
    if(atyp){
       printf("What constant value do you want to apply to the chosen intensity value?\n\n");
       scanf("%e", &val);
    }
    
    for(i=0; i < trnum; i++){
       trr = alltr + i;
       for(j=0; j < trr->num; j++){
           fp = trr->trpt + j;
           str = (iadf) ? fp->add_fld[iff] : fp->zf;
	   if(str > ADD_CHECK) str = ADD_UNDEF;
	   else if(str < rng1 || str > rng2) str = ADD_UNDEF;
           else {	   
              switch(atyp){
                  case 1:
	     	    str += val;
		     break;
		   case 2:
		     str -= val;
		     break;
		   case 3:
		     str *= val;
		     break;
		   case 4:
		     str /= val;
		     break;
                   case 5:
                     str = fabs(str);
                     break;
               }
	    
	    }
	    
	    if(iadf) fp->add_fld[iff] = str;
	    else fp->zf = str;
		
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
