#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "st_fo.h"
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"
#include "sqt.h"

#define  NFF   5


struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
int toverlap(struct tot_tr * , struct tot_tr * , long int * , long int * );
void convert_track(struct tot_tr * , int , int , int );
double ortho_dist(struct fet_pt_tr * , struct fet_pt_tr * , int , int * , int , int , VEC * );


int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

/* routine to compute differences between matched tracks assuming first
   track in the control.                                                 */

int main(int argc, char **argv)
{
   int i, j, k;
   int gpr, ipr;
   int trnum=0;
   int it1s=0, it1e=0, it2s=0, it2e=0;
   int ist=0;
   int num=0, numt=0;
   int ifnd=0;
   
   int ifu=0, iffu=0;
   
   int *nfwposu=NULL;
   int nffu=0, nfldu=0;

   long int is1=0, is2=0;

   float alat=0.0, alng=0.0;

   double dd, dp;
   double norm=0.0;
   double sum=0.0;
   double str1=0.0, str2=0.0;

    char infil[MAXCHR], outfil[MAXCHR];

   FILE *fin=NULL, *fout=NULL;

   VEC visec;
   VEC a, b, c, d;

   struct tot_tr *trr=NULL, *tr1=NULL, *tr2=NULL;
   struct fet_pt_tr *fp1=NULL, *fp2;
   struct fet_pt_tr *ff1=NULL, *ff2;


   if(argc != 3){
      printf("Usage: ensdiff [filename] [Field Id.]\n\n");
      exit(1);
   }   

   sprintf(infil, "%s", argv[1]);
   sscanf(argv[2], "%d", &ifu);

   fin = fopen(infil, "r");
   if(!fin){
      printf("Can't open file %s for read.\n\n", infil);
      exit(1);
   }

   trr = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

   fclose(fin);
   
   nffu = nff;
   nfldu = nfld;
   nfwposu = nfwpos;
   
   if(ifu < 0 || ifu > nffu){
      printf("****ERROR****, added field Id. not present in file.\n\n");
      exit(1);
   }
   
   if(ifu){
     iffu = 0;
     for(i=0; i < ifu; i++){
        if(*(nfwposu + i)) iffu += 3;
        else iffu += 1;
     }
     --iffu;
   }

   convert_track(trr, trnum, 0, 0);

   nff = nfld = NFF;
   nfwpos = (int *)calloc(nff, sizeof(int));
   mem_er((nfwpos == NULL) ? 0 : 1, nff * sizeof(int));
   for(i=0; i < nff; i++) *(nfwpos+i) = 0;

   tr1 = trr;
   for(i=0; i < tr1->num; i++){
       fp1 = tr1->trpt + i;
       if(nffu) fp1->atmp = fp1->add_fld;
       fp1->add_fld = (float *)calloc(nfld, sizeof(float));
       mem_er((fp1->add_fld == NULL) ? 0 : 1, nfld*sizeof(float));
       for(j=0; j < nfld; j++) fp1->add_fld[j] = ADD_UNDEF;
   }

   for(i=1; i < trnum; i++){
       tr2 = trr + i;
       for(j=0; j < tr2->num; j++){
           fp2 = tr2->trpt + j;
	   if(nffu) fp2->atmp = fp2->add_fld;
           fp2->add_fld = (float *)calloc(nfld, sizeof(float));
           mem_er((fp2->add_fld == NULL) ? 0 : 1, nfld*sizeof(float));
           for(k=0; k < nfld; k++) fp2->add_fld[k] = ADD_UNDEF;
       }
       
       toverlap(tr1, tr2, &is1, &is2);

/* find position along track for start and end times */

       for(j=0; j < tr1->num; j++){

           fp1 = tr1->trpt + j;

           if(fp1->time == is1) it1s = j;
           if(fp1->time == is2) it1e = j;

       } 

       for(j=0; j < tr2->num; j++){

           fp2 = tr2->trpt + j;

           if(fp2->time == is1) it2s = j;
           if(fp2->time == is2) it2e = j;

       } 

       numt = it1e - it1s + 1;
       num = it2e - it2s + 1;

       if(num != numt){
          printf("***ERROR***, track section lengths dont match for TRACK_ID's %d %d\n\n", tr1->trid, tr2->trid);
          exit(1);
       }

       str1 = str2 = 0.0;

       for(j=0; j < num; j++) {
          fp1 = tr1->trpt + it1s + j;
          fp2 = tr2->trpt + it2s + j; 

          dd = fp1->pp[0] * fp2->pp[0] + fp1->pp[1] * fp2->pp[1] + fp1->pp[2] * fp2->pp[2];
          if(fabs(dd) > 1.) dd = (dd < 0.) ? -1.0 : 1.0;
          fp2->add_fld[0] = acos(dd) / FP_PI;
	  if(!ifu) fp2->add_fld[1] = fp2->zf - fp1->zf; 
	  else fp2->add_fld[1] = fp2->atmp[iffu] - fp1->atmp[iffu];
       }

/* compute along and across track distances */

       ifnd = 1;
       ist = 1;

       for(j=0; j < num; j++){

           fp1 = tr1->trpt + it1s + j;
           fp2 = tr2->trpt + it2s + j;
 
           dd = ortho_dist(fp2, tr1->trpt, tr1->num, &ist, ifnd, 0, &visec);

           if(ist < tr1->num){

              if(dd >= 0.0) {
                 ifnd = 0;
                 ff1 = tr1->trpt + ist - 1;
                 a.x = ff1->pp[0];
                 a.y = ff1->pp[1];
                 a.z = ff1->pp[2];
                 ff2 = tr1->trpt + ist;
                 b.x = ff2->pp[0];
                 b.y = ff2->pp[1];
                 b.z = ff2->pp[2];

                 crosp(&a, &b, &c);
                 norm = sqrt(dotp(&c, &c));
                 normv(&c, norm);

                 d.x = fp2->pp[0];
                 d.y = fp2->pp[1];
                 d.z = fp2->pp[2];

                 dp = dotp(&c, &d);

/* left hand side positive, right hand side negative  */

                 fp2->add_fld[2] = (dp >= 0.0) ? dd / FP_PI : -dd / FP_PI;

                 sum = 0.0;

                 if(ist > it1s + j){
                    k=it1s + j;
                    ff1 = fp1;
                    ff2 = fp1 + 1;
                    while(k < ist - 1){
                        sum += acos(ff1->pp[0] * ff2->pp[0] + ff1->pp[1] * ff2->pp[1] + ff1->pp[2] * ff2->pp[2]);
                        ++k;
                        ++ff1;
                        ++ff2;
                    }
                    sum += acos(ff1->pp[0] * visec.x + ff1->pp[1] * visec.y + ff1->pp[2] * visec.z);

                 }

                 else {
                    k = it1s + j;
                    ff1 = fp1;
                    ff2 = fp1 - 1;
                    while(k > ist){
                        sum -= acos(ff1->pp[0] * ff2->pp[0] + ff1->pp[1] * ff2->pp[1] + ff1->pp[2] * ff2->pp[2]);
                        --k;
                        --ff1;
                        --ff2;
                    }                    
                    sum -= acos(ff1->pp[0] * visec.x + ff1->pp[1] * visec.y + ff1->pp[2] * visec.z);

                 }

                 fp2->add_fld[3] = sum / FP_PI;

                 if(ist > 1) --ist;

              }

           }

       }

       str1 = str2 = 0.0;

       for(j=0; j < num - 1; j++) {
          fp1 = tr1->trpt + it1s + j;
          fp2 = tr1->trpt + it1s + j + 1; 

          dd = fp1->pp[0] * fp2->pp[0] + fp1->pp[1] * fp2->pp[1] + fp1->pp[2] * fp2->pp[2];
          if(fabs(dd) > 1.) dd = (dd < 0.) ? -1.0 : 1.0;
          str1 += acos(dd) / FP_PI;

          fp1 = tr2->trpt + it2s + j;
          fp2 = tr2->trpt + it2s + j + 1; 

          dd = fp1->pp[0] * fp2->pp[0] + fp1->pp[1] * fp2->pp[1] + fp1->pp[2] * fp2->pp[2];
          if(fabs(dd) > 1.) dd = (dd < 0.) ? -1.0 : 1.0;
          str2 += acos(dd) / FP_PI;
          
          fp2->add_fld[4] = str2 - str1;

       }

   }

   sprintf(outfil, "%s", infil);
   strcat(outfil, "_cross-tr");
   fout = fopen(outfil, "w");
   if(!fout){
      printf("****ERROR****, unable to open file %s for write.\n\n", outfil);
      exit(1);
   }   

   meantrd(fout, trr, trnum, 's', gpr, ipr, alat, alng); 

   fclose(fout);
   
   for(i=0; i < trnum; i++){
      tr1 = trr + i;
      for(j=0; j < tr1->num; j++){
          fp1 = tr1->trpt + j;
	  if(nffu) free(fp1->atmp);
	  free(fp1->add_fld);
      }
      free(tr1->trpt);
   }
   free(nfwposu);
   free(nfwpos);
   free(trr);

   return 0;

}
