#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <splice.h>
#include "mem_er.h"
#include "m_values.h"

#define  NSAMP     0.5   /* threshold on sampling if missing values present */
#define  NSAMPLR   0.25  /* threshold for left/right sampling if missing values present */

/* function to compute the Hart cyclone diagnostics */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
void least_squares(float * , float * , float * , float * , int );


int noheader=0;

extern int nfld, nff;
extern int *nfwpos;

int main()
{

   int i, j, k, m;
   int idim=0;
   int trnum=0, trn=0;
   int ntheta=0, nr=0, nwfld=0;
   int idir=0;
   int ifdir=0;
   int ndsmth=0;
   int ifcnt=0, itpadd=0, iprojd=0;
   int ihemi=0;
   int il1=0, il2=0, iu1=0, iu2=0;
   int nwk1=0, nwk2=0, nwk=0;
   int imiss=0;
   int nsamp=0, nright=0, nleft=0;
   int nhf=0;
   int nomiss = 0;
   
   int nfldo=0, nffo=0;

   int gpr=0, ipr=0;

   long int ptnum=0, tptnum=0;

   float **zz=NULL;
   float *flat=NULL, *flng=NULL;
   float *zmax=NULL, *zmin=NULL;
   float *wk1=NULL, *wk2=NULL;
   float *press=NULL, *coslat=NULL;
   float zleft=0.0, zright=0.0, aleft=0.0, aright=0.0;
   float z1=0.0, z2=0.0;
   float hscl=1.0;
   float theta=0.0, rad=0.0;
   float vtl0=0.0, vtu0=0.0, sym0=0.0;
   float s1=0.0, s2=0.0;

   float alat=0.0, alng=0.0;

   FILE *fin=NULL, *fout=NULL;

   char filnam[MAXCHR], filout[MAXCHR];
   char ftrin[MAXCHR], ftrout[MAXCHR];
   char line[200];

   struct tot_tr *tracks=NULL, *atr=NULL;
   struct fet_pt_tr *fp=NULL; 

   printf("What is the track file to read?\n\n");
   scanf("%s", ftrin);

   fin = fopen(ftrin, "r");
   if(!fin){
      printf("***ERROR***, unable to open file %s for 'r'\n\n", ftrin);
      exit(1);
   }

   tracks = read_tracks(fin, &trn, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

   fclose(fin);
   
   fin = NULL;
   
   nffo = nff;
   nfldo = nfld;
   
   nfld += 3;
   nff += 3;
   
   
   nfwpos = (int *)realloc_n(nfwpos, nff * sizeof(int));
   mem_er((nfwpos == NULL) ? 0 : 1, nff * sizeof(int));
   *(nfwpos + nffo) = 0;
   *(nfwpos + nffo + 1) = 0;
   *(nfwpos + nffo + 2) = 0;

   for(i=0; i < trn; i++) tptnum += (tracks + i)->num;

   printf("what is the region sampled file? This needs to be the height field at several levels. \n\n");
   scanf("%s", filnam);

   fin = fopen(filnam, "r");
   if(!fin){
      printf("****ERROR****, can't open file %s for read.\n\n", filnam);
      exit(1);
   } 

   printf("What hemisphere is this? Input '0' for NH and '1' for SH.\n\n");
   scanf("%d", &ihemi);

   if(ihemi) hscl = -1.0;
   
   printf("Are there missing data value, '0' for no, '1' for yes.\n\n");
   scanf("%d", &imiss);
   if(imiss < 0 || imiss > 1){
      printf("****ERROR****, incorrect specifier for missing value option.\n\n");
      exit(1);
   }

   fgets(line, 100, fin);
   sscanf(line, "%d %ld %d %d %d %d %d %d %d %d %d", &trnum, &ptnum, &ntheta, &nr, &nwfld, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd, &ifdir);

   if(trn != trnum || tptnum != ptnum){
      printf("****ERROR****, track file and region file do not match in tracks or points.\n\n");
      exit(1);
   }

   if(!idir){
     printf("***ERROR***, regions must be sampled with storm direction.\n\n");
     exit(1);
   }

   printf("***INFORMATION***, there are %d available levels, numbered from 1.\n\n", nwfld);
   printf("What is the lower pair of levels, il1, il2?\n\n");
   scanf("%d %d", &il1, &il2);
   printf("What is the upper pair of levels, iu1, iu2?\n\n");
   scanf("%d %d", &iu1, &iu2);
   --il1; --il2; --iu1; --iu2;

   if((il1 < 0 || il1 > nwfld - 1) || (il2 < 0 || il2 > nwfld - 1) || 
      (iu1 < 0 || iu1 > nwfld - 1) || (iu2 < 0 || iu2 > nwfld - 1)){ 
      printf("***ERROR***, chosen levels are incompatable with the available data.\n\n");
      exit(1);
   }

   for(i=il1; i <= il2; i++) ++nwk1;
   for(i=iu1; i <= iu2; i++) ++nwk2; 
   nwk = (nwk1 > nwk2) ? nwk1 : nwk2;

   wk1 = (float *)calloc(nwk, sizeof(nwk));
   mem_er((wk1 == NULL) ? 0 : 1, nwk*sizeof(float));
   wk2 = (float *)calloc(nwk, sizeof(nwk));
   mem_er((wk2 == NULL) ? 0 : 1, nwk*sizeof(float));

   idim = nr * ntheta;
   nhf = idim / 2;
   zz = (float **)calloc(nwfld, sizeof(float *));
   mem_er((zz == NULL) ? 0 : 1, nwfld*sizeof(float *));

   for(i=0; i < nwfld; i++){
      zz[i] = (float *)calloc(idim, sizeof(float));
      mem_er((zz[i] == NULL) ? 0 : 1, idim*sizeof(float));
   }

   flng = (float *)calloc(ntheta, sizeof(float));
   mem_er((flng == NULL) ? 0 : 1, ntheta*sizeof(float));

   flat = (float *)calloc(nr, sizeof(float));
   mem_er((flat == NULL) ? 0 : 1, nr*sizeof(float));

   fread(flng, ntheta * sizeof(float), 1, fin);
   fscanf(fin, "%*c");
   fread(flat, nr * sizeof(float), 1, fin);
   fscanf(fin, "%*c");

   zmin = (float *)calloc(nwfld, sizeof(float));
   mem_er((zmin == NULL) ? 0 : 1, nwfld*sizeof(float));

   zmax = (float *)calloc(nwfld, sizeof(float));
   mem_er((zmax == NULL) ? 0 : 1, nwfld*sizeof(float));

   press = (float *)calloc(nwfld, sizeof(float));
   mem_er((press == NULL) ? 0 : 1, nwfld*sizeof(float));
   
   coslat = (float *)calloc(nr, sizeof(float));
   mem_er((coslat == NULL) ? 0 : 1, nr*sizeof(float));
   
   for(i=0; i < nr; i++) *(coslat + i) = cos(*(flat + i) * FP_PI);

   printf("What are the pressure levels?\n\n");
   for(i=0; i < nwfld; i++){
      printf("Input next value: ");
      scanf("%f", press+i);
      printf("\n");
   }
   
   printf("What is the output filename?\n\n");
   scanf("%s", filout);

   fout = fopen(filout, "w");
   if(!fout){
      printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
      exit(1);
   }

   for(i=0; i < trn; i++){
       atr = tracks + i;
       fprintf(fout, "TRACK_ID  %d POINT_NUM %d\n", atr->trid, atr->num);
       for(j=0; j < atr->num; j++){
          fp = atr->trpt + j;
	  
	  fp->add_fld = (float *)realloc_n(fp->add_fld, nfld * sizeof(float));
          mem_er((fp->add_fld == NULL) ? 0 : 1, nfld * sizeof(float)); 
	  
          zright=0.;
          zleft=0.;
          aright=0.;
          aleft=0.;

          for(k=0; k < nwfld; k++){
             fread(zz[k], idim * sizeof(float), 1, fin);
             fscanf(fin, "%*c");
             *(zmax + k) = 0.0;
             *(zmin + k) = 150000.;

             nsamp = 0;
             for(m=0; m < idim; m++){
                 z1 = *(zz[k] + m);  
		 if(imiss && z1 > ADD_CHECK) continue;
                 if(z1 > *(zmax + k)) *(zmax + k) = z1;
                 if(z1 < *(zmin + k)) *(zmin + k) = z1;
		 ++nsamp;
             }
	     if(((float) nsamp) / idim < NSAMP) {*(zmax + k) = ADD_UNDEF; *(zmin + k) = ADD_UNDEF;}
          }
	  
	  nright = 0;
	  nleft = 0;

          for(k=0; k < ntheta; k++){
             theta = *(flng + k);
             for(m=0; m < nr; m++){
                 rad = *(coslat + m);
                 z1 = *(zz[il1] + m * ntheta + k); 
                 z2 = *(zz[il2] + m * ntheta + k);
		 if(imiss && (z1 > ADD_CHECK || z2 > ADD_CHECK)) continue;
                 if(theta > 0.0 && theta < 180.0){
                    zright += rad * (z2 - z1);
                    aright += rad;
		    ++nright; 
                 } 
                 else if (theta > 180.0 && theta < 360.0){
                    zleft += rad * (z2 - z1);
                    aleft += rad; 
		    ++nleft;
                 }
             }
          }
	  
	  if(((float) nright) / nhf < NSAMPLR || ((float) nleft) / nhf < NSAMPLR) sym0 = ADD_UNDEF;
	  else sym0 = hscl * ((zright / aright) - (zleft / aleft));
	  
	  nomiss = 0;
          for(k=0; k < nwk1; k++){
             *(wk1 + k) = *(press + il1 + k);
	     if(*(zmax + il1 + k) > ADD_CHECK || *(zmin + il1 + k) > ADD_CHECK) {nomiss = 1; break;}
             *(wk2 + k) = *(zmax + il1 + k) - *(zmin + il1 + k);
          }

          if(!nomiss){
            least_squares(wk1, wk2, &s1, &s2, nwk1); 
            vtl0 = s1 * (*(press + il1) - *(press + il2)) / log(*(press + il1) / *(press + il2));
	  }
	  else vtl0 = ADD_UNDEF;

          nomiss = 0;
          for(k=0; k < nwk2; k++){
             *(wk1 + k) = *(press + iu1 + k);
	     if(*(zmax + iu1 + k) > ADD_CHECK || *(zmin + iu1 + k) > ADD_CHECK) {nomiss = 1; break;}
             *(wk2 + k) = *(zmax + iu1 + k) - *(zmin + iu1 + k);
          }
	  
	  if(!nomiss){
             least_squares(wk1, wk2, &s1, &s2, nwk2);
             vtu0 = s1 * (*(press + iu1) - *(press + iu2)) / log(*(press + iu1) / *(press + iu2));
	  }
	  else vtu0 = ADD_UNDEF;

          if(atr->time)fprintf(fout, "%10ld %f %f %e %e %e %e\n", fp->time, fp->xf, fp->yf, fp->zf, vtl0, vtu0, sym0);
	  else fprintf(fout, "%d %f %f %e %e %e %e\n", fp->fr_id, fp->xf, fp->yf, fp->zf, vtl0, vtu0, sym0);
	 
	  *(fp->add_fld + nfldo) = vtl0;
	  *(fp->add_fld + nfldo + 1) = vtu0;
	  *(fp->add_fld + nfldo + 2) = sym0; 
	  
       
       }
       
   }

   fclose(fin);
   fclose(fout);
   
   fout = NULL;
   
   strncpy(ftrout, ftrin, MAXCHR);
   strcat(ftrout, ".hart"); 
   
   fout = fopen(ftrout, "w");
   if(!fout){
      printf("***ERROR***, unable to open file %s for 'w'\n\n", ftrout);
      exit(1);
   }

   meantrd(fout, tracks, trn, 's', gpr, ipr, alat, alng);

   fclose(fout);

   for(i=0; i < nwfld; i++) free(*(zz+i));
   free(zz)
   free(zmin);
   free(zmax);
   free(press);
   free(wk1);
   free(wk2);
   free(coslat);
   
   for(i=0; i < trn; i++){
   
   } 

   return 0;

}

void least_squares(float *wk1, float *wk2, float *s1, float *s2, int nk)
{
   int i=0;

   float am=0.0, bm=0.0;
   float aa=0.0, ab=0.0;
   float a=0.0, b=0.0;

   for(i=0; i < nk; i++){
      am += *(wk1 + i);
      bm += *(wk2 + i);
   }
   am /= nk;
   bm /= nk;

   for(i=0; i < nk; i++){
      a = (*(wk1 + i) - am);
      b = (*(wk2 + i) - bm);
      aa += a * a;
      ab += a * b; 
   } 
   aa /= nk;
   ab /= nk;
   *s1 = ab / aa;
   *s2 = bm - *s1 * am;

    return;
}
