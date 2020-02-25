#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <splice.h>
#include "mem_er.h"
#include "m_values.h"

/* function to partition track by reverse shear */ 

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );


int noheader=0;

extern int nfld, nff;
extern int *nfwpos;

int main(void )
{
   int i=0, j=0;
   int ii=0;
   int trnum=0;
   int gpr=0, ipr=0;
   int ithw1=0, ithw2=0;
   int ith1=0, ith2=0;
   int iths1=0, iths2=0;
   int its1=0, its2=0;
   int irs=0;
   int ismp=1, ioff=0;
   
   int itm=0;
   
   int ftype='s';
   
   float alat=0.0, alng=0.0;
   
   float ang1=0.0, ang2=0.0;
   
   float strm=0.0;
   
   double ang=0.0, nrmt=0.0, nrms=0.0;

   FILE *fin=NULL, *fout=NULL;

   char ftrin[MAXCHR];
   
   struct tot_tr *tracks=NULL, *trackn=NULL, *atr=NULL;
   struct fet_pt_tr *fp1=NULL, *fp2=NULL; 

   printf("What is the track file to read?\n\n");
   scanf("%s", ftrin);

   fin = fopen(ftrin, "r");
   if(!fin){
      printf("***ERROR***, unable to open file %s for 'r'\n\n", ftrin);
      exit(1);
   }

   tracks = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
   fseek(fin, 0L, SEEK_SET);
   trackn = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

   fclose(fin);

   printf("Which two additional fields are required for the U and V of the thermal wind?\n\n");
   scanf("%d %d", &ithw1, &ithw2);
   if((ithw1 <= 0 || ithw1 > nff) || (ithw2 <= 0 || ithw2 > nff)){
      printf("****ERROR****,incorrect values for additional fields, %d %d\n\n", ithw1, ithw2);
      exit(1);
   }
   

   ith1 = 0;
   for(i=0; i < ithw1; i++){
      if(*(nfwpos + i)) ith1 += 3;
      else ith1 += 1;
   }
   --ith1;
   
   ith2 = 0;
   for(i=0; i < ithw2; i++){
      if(*(nfwpos + i)) ith2 += 3;
      else ith2 += 1;
   }
   --ith2;  
   
   printf("Which two additional fields are required for the U and V of the steering wind?\n\n");
   scanf("%d %d", &iths1, &iths2);
   if((iths1 <= 0 || iths1 > nff) || (iths2 <= 0 || iths2 > nff)){
      printf("****ERROR****,incorrect values for additional fields, %d %d\n\n", iths1, iths2);
      exit(1);
   }
   

   its1 = 0;
   for(i=0; i < iths1; i++){
      if(*(nfwpos + i)) its1 += 3;
      else its1 += 1;
   }
   --its1;
   
   its2 = 0;
   for(i=0; i < iths2; i++){
      if(*(nfwpos + i)) its2 += 3;
      else its2 += 1;
   }
   --its2;
   
   printf("What angular range to you want to use to define reverse shear systems, input ang1 and ang2 between 0 and 180.\n\n");
   scanf("%f %f", &ang1, &ang2);
   
   if((ang1 < 0.0 || ang1 > 180.0) || (ang2 < 0.0 || ang2 > 180.0)) {
      printf("****ERROR****, angles incorrectly specified.\n\n");
      exit(1);  
   }
   
   ang1 *= FP_PI;
   ang2 *= FP_PI;

/*   
   printf("What time step sampling and offset is required, a value of 1 for sampling is every time step?\n\n");
   scanf("%d %d", &ismp, &ioff);
   
   if(ismp < 1 || ioff < 0){
      printf("****ERROR****, sampling must be >= 1 and offset must be >= 0.\n\n");
      exit(1);
   }
*/
 
   ++nff;
   ++nfld;
   
   nfwpos = (int *)realloc_n(nfwpos, nff*sizeof(int));
   mem_er((nfwpos == NULL) ? 0 : 1, nff * sizeof(int));
   *(nfwpos + nff - 1) = 0;
   
   for(i=0; i < trnum; i++){
       atr = tracks + i;
       
       for(j=0; j < atr->num; j++){
           fp1 = atr->trpt + j;
           fp2 = (trackn + i)->trpt + j;
	   fp1->add_fld = (float *)realloc_n(fp1->add_fld, nfld * sizeof(float));
	   mem_er((fp1->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
	   fp2->add_fld = (float *)realloc_n(fp2->add_fld, nfld * sizeof(float));
	   mem_er((fp2->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));	   	   
       }

/* perform sampling */

/*       if(ismp > 1){ */
       
/* find first point */
/*
          ipos = altr->num + 1;

          for(j=0; j < altr->num; j++) {
              fp1 = atr->trpt + j;
              idf = fp1->fr_id - 1 - ioff;
              if(!(idf % ismp) && idf >= 0) {ipos = j; itrid = (idf / ismp) + 1; break;}
          }

          if(ipos >= altr->num) continue;  
       }
*/
       irs = 0;
       
/* find time of max. intensity */

       strm = 0.0;
       itm = 0;

       for(j=0; j < atr->num; j++){
           fp1 = atr->trpt + j;
           if(fp1->zf > strm) {strm = fp1->zf; itm = j;}     
       }
       
       for(j=0; j < atr->num; j++){
          fp1 = atr->trpt + j;
	  fp2 = (trackn + i)->trpt + j;
	 
	  nrmt = sqrt(fp1->add_fld[ith1] * fp1->add_fld[ith1] + fp1->add_fld[ith2] * fp1->add_fld[ith2]);
	  nrms = sqrt(fp1->add_fld[its1] * fp1->add_fld[its1] + fp1->add_fld[its2] * fp1->add_fld[its2]);
	  ang = acos((fp1->add_fld[ith1] * fp1->add_fld[its1] + fp1->add_fld[ith2] * fp1->add_fld[its2]) / (nrmt * nrms));
	  
	  fp1->add_fld[nfld-1] = ang / FP_PI;
	  fp2->add_fld[nfld-1] = ang / FP_PI;
	  
	  ii = itm - j;
	  
	  if(ang >= ang1 && ang <= ang2 && ii > 0 && ii <= 4 && (j + 1)%2) irs = 1;

       }
       
       if(irs) (trackn + i)->num = 0;    
       else (tracks + i)->num = 0;
   
   }
   
   fout = fopen("ff_trs_yes", "w");
   if(!fout){
      printf("***ERROR***, unable to open file %s for 'w'\n\n", "ff_trs_yes");
      exit(1);
   }

   meantrd(fout, tracks, trnum, ftype, gpr, ipr, alat, alng);

   fclose(fout);   

   fout = fopen("ff_trs_no", "w");
   if(!fout){
      printf("***ERROR***, unable to open file %s for 'w'\n\n", "ff_trs_no");
      exit(1);
   }

   meantrd(fout, trackn, trnum, ftype, gpr, ipr, alat, alng);

   fclose(fout); 
   
   for(i=0; i < trnum; i++){
      atr = tracks + i;
      if(nff){
         for(j=0; j < atr->num; j++) free((atr->trpt + j)->add_fld);
      }
      free(atr->trpt);
      
      atr = trackn + i;
      if(nff){
         for(j=0; j < atr->num; j++) free((atr->trpt + j)->add_fld);
      }
      free(atr->trpt);      
      
   }
   
   free(tracks);
   free(trackn);

   return 0;
}
