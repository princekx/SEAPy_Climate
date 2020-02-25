#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "mem_er.h"
#include "splice.h"

#define MAXCHR 500
#define  SMALL -1.0e+12
#define  LARGE 1.0e+12

/* program to compute minima and maxima over spherical cap sampled regions */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
void maxmin(float * , int * , int *, int , int , int );

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int noheader=0;

int main(void )
{

   int i, j, k;
   
   int trnum=0;
   int gpr=0, ipr=0;
   
   int nffo=0, nfldo=0;
   
   int irtrn=0, irnr=0, irnth=0, irnf=0;
   int ndsmth=0, ifcnt=0, itpadd=0, iprojd=0;
   int idir=0;
   int irdim=0;
   int npt=0;
   int igsmp=0;
   int ifdir=0;
   
   int inth=0, innr=0;
   int irmax=0;
   
   int imxmn=0;
   int ilev=0, nlev=0, ilevt=0;
   
   long int iptnum;
   
   off_t place1, place2, blklen;
   
   float alat=0.0, alng=0.0;
   float *sradf=NULL;
   float *slng=NULL, *slat=NULL;
   
   float rmax=0.0;
   
   FILE *fin=NULL, *fout=NULL;
   FILE *frad=NULL;
   
   char filin[MAXCHR];
   char filout[MAXCHR];
   char filrad[MAXCHR];
   
   char line[MAXCHR];
   
   struct tot_tr *tracks=NULL, *atr=NULL;
   struct fet_pt_tr *fp=NULL;
   
   printf("What is the track file to read?\n\n");
   scanf("%s", filin);

   fin = fopen(filin, "r");
   if(!fin){
      printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
      exit(1);
   }

   tracks = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

   fclose(fin);
   
   nffo = nff;
   nfldo = nfld;
   
   
   printf("Do you want to find minima or maxima, '0' for minima, '1' for maxima.\n\n");
   scanf("%d", &imxmn);
   if(imxmn < 0 || imxmn > 1){
      printf("****ERROR****, incorrect specifier for min or max.\n\n");
      exit(1);
   }
  
   printf("What is the along track radial field file associated with the track file?\n\n");
   scanf("%s", filrad);
   frad = fopen(filrad, "r");
   if(!frad){
     printf("***ERROR***, unable to open file %s for 'r'\n\n", filrad);
     exit(1);
   }
   fgets(line, MAXCHR, frad);    
   sscanf(line, "%d %ld %d %d %d %d %d %d %d %d %d %d", &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd, &ifdir, &igsmp);
   
   if(igsmp){
      printf("****ERROR****, this program does not currently support rotated lat-lng grid.\n\n");
      exit(1);
   }

   printf("****INFORMATION****, number of fields availbale is %d\n", irnf);
   
   if(irtrn != trnum){
      printf("****ERROR****, track file and radial field data not compatable: different numbers of tracks.\n\n");
      exit(1);
   }
   
   if(irnf > 1){
      printf("Are all fields required or a selected field, 'a' for all or 's' for select.\n\n");
      scanf("\n");
      if(getchar() == 's'){
         printf("Which level is required?\n\n");
	 scanf("%d", &ilev);
	 --ilev;
	 nlev = 1;
	 if(ilev < 0 || ilev >= irnf){
	   printf("****ERROR****, level id out of range.\n\n");
	   exit(1);
	 }
	 ilevt = 1;
      }
      else {nlev = irnf; ilevt = 0;}
   }
   
   else {nlev = 1; ilevt = 0;}
   
   printf("What search radius is required in degrees?\n\n");
   scanf("%f", &rmax);

   nff += nlev;
   nfld += 3 *nlev;
   nfwpos = (int *)realloc_n(nfwpos, nff*sizeof(int));
   mem_er((nfwpos == NULL) ? 0 : 1, nff * sizeof(int));
   for(i=0; i < nlev; i++) *(nfwpos + nffo + i) = 1;
   
   irdim = irnr * irnth;
   
   sradf = (float *)calloc(irdim, sizeof(float));
   mem_er((sradf == NULL) ? 0 : 1, irdim*sizeof(float));
   slng = (float *)calloc(irnth, sizeof(float));
   mem_er((sradf == NULL) ? 0 : 1, irnth*sizeof(float));
   slat = (float *)calloc(irnr, sizeof(float));
   mem_er((sradf == NULL) ? 0 : 1, irnr*sizeof(float));  
   
   fread(slng, irnth*sizeof(float), 1, frad);
   fscanf(frad, "%*c");
   fread(slat, irnr*sizeof(float), 1, frad);
   fscanf(frad, "%*c");
   
   for(i=0; i < irnr; i++){
       if(90.0 - *(slat + i) <= rmax) irmax = i;
   }
   
/* calculate data block size */

   place1 = ftello(frad);
   fread(sradf, irdim*sizeof(float), 1, frad);
   fscanf(frad, "%*c");
   place2 = ftello(frad);
   blklen = place2 - place1;
   fseeko(frad, place1, SEEK_SET);
   
   for(i=0; i < trnum; i++){
      atr= tracks + i;

      for(j=0; j < atr->num; j++){
      
          fp = atr->trpt + j;
	  
	  fp->add_fld = (float *)realloc_n(fp->add_fld, nfld*sizeof(float));
	  
          if(!ilevt){
	     for(k=0; k < nlev; k++){
	        fread(sradf, irdim*sizeof(float), 1, frad);
                fscanf(frad, "%*c");
		maxmin(sradf, &inth, &innr, irnth, irmax, imxmn);
		fp->add_fld[nfldo + 3 * k] = *(slng + inth);
		fp->add_fld[nfldo + 3 * k + 1] = *(slat + innr);
		fp->add_fld[nfldo + 3 * k + 2] = *(sradf + innr * irnth + inth);
	     }
	  }
	  else {
	     fseeko(frad, (npt + j) * irnf * blklen + ilev * blklen + place1, SEEK_SET);
	     fread(sradf, irdim*sizeof(float), 1, frad);
             fscanf(frad, "%*c");
             maxmin(sradf, &inth, &innr, irnth, irmax, imxmn);
	     fp->add_fld[nfldo] = *(slng + inth);
             fp->add_fld[nfldo + 1] = *(slat + innr);
	     fp->add_fld[nfldo + 2] = *(sradf + innr * irnth + inth);
	     
	  }
      
      }
      
      npt += atr->num;
   
   }
   
   
   fclose(frad);
   
   printf("What is the output file name?\n\n");
   scanf("%s", filout);

   fout = fopen(filout, "w");
   if(!fout){
      printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
      exit(1);
   }
   
   meantrd(fout, tracks, trnum, 's', gpr, ipr, alat, alng);
   
   fclose(fout);
   
   free(sradf);
   free(slng);
   free(slat);

   return 0;
}

void maxmin(float *sradf, int *inth, int *innr, int nth, int irmax, int imxmn)
{

   int i, j;
   
   float val=0.0, vmxmn=0.0;
   
   if(!imxmn){
      vmxmn = LARGE;
      for(i=0; i <= irmax; i++){
          for(j=0; j < nth; j++){
	      val = *(sradf + i * nth + j);
	      if(val < vmxmn) {
	         vmxmn = val;
	         *inth = j;
	         *innr = i;
	      }
	  }
      }
      
   }
   else {
      vmxmn = SMALL;
      for(i=0; i <= irmax; i++){
          for(j=0; j < nth; j++){
	      val = *(sradf + i * nth + j);
	      if(val > vmxmn) {
	         vmxmn = val;
	         *inth = j;
	         *innr = i;
	      }
	  }
      }
   }

   return;
    
}
