#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "mem_er.h"
#include "splice.h"

/* program to take data from tcident and compute meadian fields */

float *write_fields_netcdf(char * , double ** , float * , float * , float * , int , int , int , int , int , int , int , long int );
float binmedian(int , float * );

int main(void)
{

   int i=0, j=0, k=0;
   int idim=0;
   int ntr=0;
   int nth=0, nrr=0, nnf=0;
   int iext=0, imxmn=0;
   
   int nn=0;
   
   long int ctime=0;
   
   off_t place1, place2, blklen;
   
   FILE *infil=NULL;
   
   char filin[MAXCHR];
   char filout[MAXCHR];
   char line[MAXCHR];
   
   float *slng=NULL, *slat=NULL;
   float **sradl=NULL;
   
   float *dser=NULL;
   
   float val=0.0;
   
   double **sradt=NULL;
   
   printf("What sampling file is required?\n\n");
   scanf("%s", filin);
   
   infil = fopen(filin, "r");
   if(!infil){
      printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
      exit(1);
   }
   fgets(line, MAXCHR, infil);

   sscanf(line, "%d %d %d %d %d %d", &ntr, &nth, &nrr, &nnf, &iext, &imxmn);

   idim = nth * nrr;
   
   slng = (float *)calloc(nth, sizeof(float));
   mem_er((slng == NULL) ? 0 : 1, nth * sizeof(float));
   slat = (float *)calloc(nrr, sizeof(float));
   mem_er((slat == NULL) ? 0 : 1, nrr * sizeof(float));
   sradt = (double **)calloc(nnf, sizeof(double *));
   mem_er((sradt == NULL) ? 0 : 1, nnf * sizeof(double *));
   
   for(i=0; i < nnf; i++){
       *(sradt + i) = (double *)calloc(idim, sizeof(double));
       mem_er((*(sradt + i) == NULL) ? 0 : 1, idim * sizeof(double));
   }
      
   sradl = (float **)calloc(ntr, sizeof(float *));
   mem_er((sradl == NULL) ? 0 : 1, ntr * sizeof(float *));
   
   for(i=0; i < ntr; i++){
       *(sradl + i) = (float *)calloc(idim, sizeof(float));
       mem_er((*(sradl + i) == NULL) ? 0 : 1, idim * sizeof(float));
   }
   
   dser = (float *)calloc(ntr+1, sizeof(float));
   mem_er((dser == NULL) ? 0 : 1, ntr * sizeof(float));
   
   fread(slng, nth*sizeof(float), 1, infil);
   fscanf(infil, "%*c");
   fread(slat, nrr*sizeof(float), 1, infil);
   fscanf(infil, "%*c");
   
   place1 = ftello(infil);
   fread(*sradl, idim*sizeof(float), 1, infil);
   fscanf(infil, "%*c");
   place2 = ftello(infil);
   blklen = place2 - place1;
   fseeko(infil, place1, SEEK_SET);     
   
   for(i=0; i < nnf; i++){         /* loop over number of fields/levels */

      for(j=0; j < ntr; j++){      /* loop over number of time samples to read for each level */
          fseeko(infil, (j * nnf + i) * blklen + place1, SEEK_SET);
          fread(*(sradl + j), idim*sizeof(float), 1, infil); 
      }
      
      for(j=0; j < idim; j++){
          nn = 0;
          for(k=0; k < ntr; k++){
	     val = *(*(sradl + k) + j);
	     if(val > ADD_CHECK) continue;
	     *(dser + nn) = val;
	     ++nn;
	  }
	  
	  if(!nn) *(*(sradt + i) + j) = ADD_UNDEF;
	  else {
	     if(! (nn % 2)) {           /* binmedian needs an odd number to work with at the moment */
	        *(dser + nn) = *(dser + nn - 1);
		++nn;
	     }
  
	     *(*(sradt + i) + j) = binmedian(nn, dser);

	  }
      
      }
      
      
   
   } 
   
   printf("What output file is required?\n\n");
   scanf("%s", filout);
   
   write_fields_netcdf(filout, sradt, slng, slat, NULL, nth, nrr, nnf, iext, imxmn, 0, 0, ctime);

   free(dser);
   for(i=0; i < nnf; i++) free(*(sradt + i));
   free(sradt);   
   
   for(i=0; i < ntr; i++) free(*(sradl + i));
   free(sradl);
   free(sradt);
   free(slng);
   free(slat);


   return 0;
}

float *leveldata(int nlev)
{
    int i;
    float *levdat=NULL;

    levdat = (float *)calloc(nlev, sizeof(float));
    mem_er((levdat == NULL) ? 0 : 1, nlev * sizeof(float));
    printf("Use default level values, 'y' or 'n'.\n\n");
    scanf("\n");
    if(getchar() == 'y'){
      for(i=0; i < nlev; i++){
          *(levdat + i) = i+1;
      }
    }
    else {
       printf("Specify the level values in the correct order.\n\n");
       for(i=0; i < nlev; i++){
           printf("Value for level %d: ", i+1);
           scanf("%f", levdat + i);
           printf("\n");
       }
    }

    return levdat;
}
