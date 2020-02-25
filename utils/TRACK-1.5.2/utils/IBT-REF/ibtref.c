#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <splice.h>
#include "mem_er.h"
#include "m_values.h"

/* function to reformat tracks into IBTrACS format */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );


int noheader=0;

extern int nfld, nff;
extern int *nfwpos;

int main(int argc, char *argv[])
{

   int i=0, j=0;
   int trnum=0;
   int gpr=0, ipr=0;
   
   int imslp=0, iwind=0;
   int imslpp=0, iwindp=0;
   
   int nff1=0, nfld1=0;
   
   int *nfwpos1=NULL;
   
   float alat=0.0, alng=0.0;

   char com[]="Usage: ibtref trfilein trfileout mslp-id wind-id";
   
   FILE *fin=NULL, *fout=NULL;
   
   struct tot_tr *tracks1=NULL, *tracks2=NULL, *atr1=NULL, *atr2=NULL;
   struct fet_pt_tr *fp1=NULL, *fp2=NULL;

   if(argc != 5){

     printf("****ERROR***, %s \n", com);
     exit(1);


   }  
   
   sscanf(argv[3], "%d", &imslp);
   sscanf(argv[4], "%d", &iwind);


   fin = fopen(argv[1], "r");
   if(!fin){
      printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[1]);
      exit(1);
   }

   tracks1 = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
   
   nff1 = nff;
   nfld1 = nfld;
   nfwpos1 = nfwpos;
   nfwpos = NULL;
   
   fseeko(fin, 0L, SEEK_SET);

   tracks2 = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
   
   nff = 1;
   nfld = 1;
   nfwpos = (int *)realloc_n(nfwpos, nff*sizeof(int));
   *nfwpos = 0;


   fclose(fin);   
   
   
   if((imslp < 0 || imslp > nff1) || (iwind < 0 || iwind > nff1)){
      printf("****ERROR****,incorrect values for additional fields, %d %d\n\n", imslp, iwind);
      exit(1);
   }
   

   imslpp = 0;
   if(imslp){
      for(i=0; i < imslp; i++){
         if(*(nfwpos1 + i)) imslpp += 3;
         else imslpp += 1;
      }
      --imslpp;
   }
   
   iwindp = 0;
   if(iwind){
      for(i=0; i < iwind; i++){
         if(*(nfwpos1 + i)) iwindp += 3;
         else iwindp += 1;
      }
      --iwindp;
   }
   
   for(i=0; i < trnum; i++){
       atr1 = tracks1 + i;
       atr2 = tracks2 + i;
       
       for(j=0; j < atr1->num; j++){
          fp1 = atr1->trpt + j;
	  fp2 = atr2->trpt + j;
	  
	  fp2->add_fld = (float *)realloc_n(fp2->add_fld, nfld * sizeof(float));
	  if(imslpp) fp2->zf = fp1->add_fld[imslpp];
	  if(iwindp) fp2->add_fld[0] = fp1->add_fld[iwindp];
       
       }
       
      
   
   }
   
   fout = fopen(argv[2], "w");
   if(!fout){
      printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[2]);
      exit(1);
   }

   meantrd(fout, tracks2, trnum, 's', gpr, ipr, alat, alng);

   fclose(fout);

   return 0;
}
