#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include <sys/types.h>
#include "mem_er.h"
#include "grid.h"
#include "file_handle.h"
#include "files_out.h"
#include "pp.h"
#include "utf.h"
#include "netcdf_info.h"

#define  NCHRB  30
#define BANDLIM 0.1

/* function to perform a spatial spectral filter on an limited area instanaeous field 
   i.e. fields from Limited Area Models                                                */
   
#ifdef  FFTLIB

#ifdef  NOUNDERSCORE

void setgpfa(double * , int * , int * );

#else

void setgpfa_(double * , int * , int * );

#endif

#endif
   
   
void write_header(GRID * , FILE * );
void dfct2d(complex * , int , int , int , int , int , int , int );
int read_field(float * ,float * , float , FILE * , int , int , int , int , int );
int powi(int , int );
int temptest(int * , int , int );

/* arrays for Temperton FFT */

int *njx=NULL, *njy=NULL;

double *dfreal=NULL, *dfimag=NULL;
double *dtrgx=NULL, *dtrgy=NULL;

/* ------------------------ */

extern GRID *gr;
extern char *chrfld;

extern char *fext;

extern int iext;

extern int form;

void limited_area_filter(FILE *fdat, int fr1, int fri, int frl)
{

    int i=0, j=0, k=0;
    int dim=gr->ix * gr->iy;
    int nxp=0, dimt=0;
    int nf=0;
    int ifr=0;
    
    int nexpx=0, nexpy=0;
    int nx=0, ny=0;
    int nmin=0;
    int fres='n';    
    int nband=0;
    int shf='y';
    int ifftyp=0;
    int dmtgx=0, dmtgy=0;
    int nmax=0;
    int icent=0;
    
    int n1x=0, n2x=0, n1y=0, n2y=0;
    int ndx=0, ndy=0;
    int ndx2=0, ndy2=0;
    int irzo=0;
    
    int itrig=1;
    
    int nxx[3], nyy[3];
    
    int *karray=NULL;
    int *kbands=NULL;
    
    off_t place1, place2, chrnum=0;
    
    float gres=0.0;
    float cut=0.0;
    float fmin=0.0, fmax=0.0;
    float fval=0.0;
    
    float *ap=NULL;
    float *bands=NULL;
    
    double alpha=0.0, py=0.0;
    double k0=-1.0;
    
    double *alphamin=NULL;
    double **fmask=NULL;
    
    double savg=0.0;
    
    complex *ac=NULL, *acc=NULL;
    
    FILE **filspec=NULL;
    
    char **filename=NULL, band_id[NCHRB];

    njx = nxx;
    njy = nyy;

    if(form == 4){

       if(((NETCDF_INFO *)fdat)->imiss){
          printf("****WARNING****, there are possible missing values in the data.\r\n"
	          "                do you want to continue, 'y' or 'n'.          \n\n");
	  scanf("\n");
	  if(getchar() != 'y') return;
       }

    }
    
    nx = gr->ix;
    ny = gr->iy;
    
#ifdef  FFTLIB

    printf("The default FFT requires powers of 2 data lengths greater than or equal to the actual data length, \r\n"
           "this can be expensive for large fields.                                                            \r\n"
	   "A prime factor FFT is also available, should this be used instead of the default, 'y' or 'n'.      \n\n");
    scanf("\n");
    if(getchar() == 'y') ifftyp = 1;
    
    
    if(ifftyp){
       if(nx % 2) ++nx;
       if(ny % 2) ++ny;
       
       if(!temptest(njx, nx, 1)){
	  
	  n1x = nx;
          while(!temptest(njx, n1x, 0)) n1x -= 2;
	  n2x = nx;
          while(!temptest(njx, n2x, 0)) n2x += 2;
	  
	  printf("The next lower and higher valid values of nx for Temperton are %d %d\n", n1x, n2x);
          printf("The transform array dimension nx is reset to %d\n\n", n2x);
	  nx = n2x;
	  if(!temptest(njx, nx, 1)) exit(1);     
       }
       
       if(!temptest(njy, ny, 1)) {
	  
	  n1y = ny;
          while(!temptest(njy, n1y, 0)) n1y -= 2;
	  n2y = ny;
          while(!temptest(njy, n2y, 0)) n2y += 2;
	  
	  printf("The next lower and higher valid values of ny for Temperton are %d %d\n", n1y, n2y);	  
	  printf("The transform array dimension ny is reset to %d\n\n", n2y);
	  ny = n2y;
	  if(!temptest(njy, ny, 1)) exit(1);
	  	       
       }   
       
       printf("Do you want to zero pad the data, 'y' or 'n'\n\n");
       scanf("\n");
       if(getchar() == 'y'){
          printf("What values you do you want for nx and ny?\n\n");
	  scanf("%d %d", &nx, &ny);
          if(nx % 2) ++nx;
          if(ny % 2) ++ny;	  
	  if(nx < gr->ix || ny < gr->iy){
	     printf("****ERROR****, nx or ny too small.\n\n");
	     exit(1);
	  }
	  if(!temptest(njx, nx, 1)) exit(1);
	  if(!temptest(njy, ny, 1)) exit(1);
       } 
       
       if(nx > gr->ix || ny > gr->iy){
         printf("Do you want to replace zero padding with real values, '1' for yes or '0' for no.\n\n");
         scanf("%d", &irzo);
	 if(irzo < 0 || irzo > 1){
	    printf("****ERROR****, incorrect input for option.\n\n");
	    exit(1);
	 }
       }
       
       nmax = (nx < ny) ? ny : nx;
       
       dmtgx = powi(2, njx[0]) + powi(3, njx[1]) + powi(5, njx[2]);
       dmtgx *= 2;
       
       dmtgy = powi(2, njy[0]) + powi(3, njy[1]) + powi(5, njy[2]);
       dmtgy *= 2;
       
       dfreal = (double *)calloc(nmax, sizeof(double));
       mem_er((dfreal == NULL) ? 0: 1, nmax * sizeof(double));
       
       dfimag = (double *)calloc(nmax, sizeof(double));
       mem_er((dfimag == NULL) ? 0: 1, nmax * sizeof(double));
       
       dtrgx = (double *)calloc(dmtgx, sizeof(double));
       mem_er((dtrgx == NULL) ? 0: 1, dmtgx * sizeof(double));
       
       dtrgy = (double *)calloc(dmtgy, sizeof(double));
       mem_er((dtrgy == NULL) ? 0: 1, dmtgy * sizeof(double));
       
#ifdef  NOUNDERSCORE

       setgpfa(dtrgx, &nx, njx);
       setgpfa(dtrgy, &ny, njy);

#else

       setgpfa_(dtrgx, &nx, njx);
       setgpfa_(dtrgy, &ny, njy);

#endif
    
    }

#endif
    
    printf("What is the resolution of the data in kilometers? \r\n\n"
           "Note data grid is assumed to be isotropic.          \n\n");
	   
    scanf("%f", &gres);

    ap=(float *)calloc(dim, sizeof(float));
    mem_er((ap == NULL) ? 0: 1, dim * sizeof(float));
    
/* determine padding for power of 2 FFT */
    
    if(!ifftyp) {
       
       printf("What powers of 2 are required for transform in the X and Y directions?\n\n");
       scanf("%d %d", &nexpx, &nexpy);
       
       nxp = (int)(log((float)nx)/log(2.0));
       if(powi(2, nxp) < nx) {
          ++nxp;
          nx = powi(2, nxp);
       }
       if(nexpx < nxp) nexpx = nxp;
       if(nxp < nexpx) nx = powi(2, nexpx); 
       
       nxp = (int)(log((float)ny)/log(2.0));
       if(powi(2, nxp) < ny) {
          ++nxp;
          ny = powi(2, nxp);
       }
       if(nexpy < nxp) nexpy = nxp;
       if(nxp < nexpy) ny = powi(2, nexpy); 
    
       printf("****INFORMATION****, powers of 2 for transform are %d %d\n", nexpx, nexpy);
    
    }

    nmin = (nx < ny) ? nx : ny;
       
    dimt = nx * ny;     
    
    ndx = (nx - gr->ix) / 2;
    ndy = (ny - gr->iy) / 2; 
    ndx2 = nx - gr->ix - ndx;
    ndy2 = ny - gr->iy - ndy;
    
    if(irzo){
       if(nx - gr->ix > gr->ix) {
          printf("Too much padding to replace in X-direction.\n\n");
	  exit(1);
       }
       if(ny - gr->iy > gr->iy) {
          printf("Too much padding to replace in Y-direction.\n\n");
	  exit(1);
       }       
    }
       
    ac = (complex *)calloc(dimt, sizeof(complex));
    mem_er((ac == NULL) ? 0 : 1, dimt * sizeof(complex)); 
    
    acc = (complex *)calloc(dimt, sizeof(complex));
    mem_er((acc == NULL) ? 0 : 1, dimt * sizeof(complex));      
    
    printf("Center data before applying transforms, 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y') icent=1;
    
    printf("How many spectral bands are required?\n");
    scanf("%d", &nband);

    if(nband <= 0){
       printf("****ERROR****, number of bands must be greater than zero\n\n");
       exit(1);
    }
    
    bands = (float *)calloc(nband+1, sizeof(float));
    mem_er((bands == NULL) ? 0 : 1, (nband+1) * sizeof(float));
    
    kbands = (int *)calloc(nband+1, sizeof(int));
    mem_er((bands == NULL) ? 0 : 1, (nband+1) * sizeof(int));

    printf("Input the band boundaries in kilometers, high to low.\n");
    
    for(i=0; i <= nband; i++){

        printf("Boundary %d = ", i+1);
        scanf("%f", bands+i);
	printf("\n");
	
	if(i){
	   if(*(bands + i) > *(bands + i - 1)){
	      printf("****ERROR****, band boundaries must be in decreasing magnitude. \n\n");
	      exit(1);
	   }
	}
	if(*(bands + i) < BANDLIM) {
	   printf("****WARNING****, band limit to small, resetting to default.\n\n");
	   *(bands + i) = BANDLIM;
	}

/*        if(!i) *(kbands + i) = (int) ceil(2.0 * gres * (double) nmin / *(bands + i));
	else *(kbands + i) = (int) floor(2.0 * gres * (double) nmin / *(bands + i));     */
	
	*(kbands + i) = (int)((2.0 * gres * (double) nmin / *(bands + i)) + 0.5); 
	
	printf("%d %f\n", *(kbands + i), 2.0 * gres * (double) nmin / *(kbands + i));
	
	
	
    } 
    
/* setup filter masks */

    fmask = (double **)calloc(nband, sizeof(double *));
    mem_er((fmask == NULL) ? 0 : 1, nband*sizeof(double *));
    
    for(i=0; i < nband; i++){
        *(fmask + i) = (double *)calloc(dimt, sizeof(double));
        mem_er((*(fmask+i) == NULL) ? 0 : 1, dimt*sizeof(double));
	memset(*(fmask + i), 0, dimt * sizeof(double));
    }

    printf("Do you want to apply further filtering in the form of the Hoskins filter, 'y' or 'n'\n");
    scanf("\n");

    if((shf=getchar()) == 'y'){

       printf("What is the cutoff constant value?\r\n"
              "Typicaly a value of 0.1 is used.  \n\n");

       scanf("%f", &cut);
       
       k0 = (float) (*(kbands + nband)) / pow(-log(cut), 0.25);
       
    }
    
    alphamin = (double *)calloc(nmin, sizeof(double));
    mem_er((alphamin == NULL) ? 0 : 1, nmin*sizeof(double));
    
    for(i=0; i < nmin; i++) *(alphamin + i) = (double) i / (double) nmin;
    
    karray = (int *)calloc(dimt, sizeof(int));
    mem_er((karray == NULL) ? 0 : 1, dimt*sizeof(int));
    
    for(i=0; i < ny; i++){
        py = pow((double) i / (double) ny, 2);
        for(j=0; j < nx; j++){
	    alpha = sqrt(py + pow((double) j / (double) nx, 2));
	    for(k=0; k < nmin; k++){
	        if(alpha >= *(alphamin + k)){ 
		   *(karray + i*nx + j) = k;
		}
	    }
	
	}

    }
    
    if(shf == 'y') {
    
       for(i=0; i < nband; i++){
           if(!i){
             for(j=0; j < dimt; j++){
	         if(*(karray + j) >= *(kbands + i) && *(karray + j) <= *(kbands + i + 1))
	           *(fmask[i] + j) = exp(-pow((double)(*(karray + j)) / k0, 4.0));
                 else *(fmask[i] + j) = 0.0;
	     } 
	   }
	   else{
             for(j=0; j < dimt; j++){
	         if(*(karray + j) > *(kbands + i) && *(karray + j) <= *(kbands + i + 1))
	            *(fmask[i] + j) = exp(-pow((double)(*(karray + j)) / k0, 4.0));
                 else *(fmask[i] + j) = 0.0;
	     }	
	   }
       }
    
    }
    
    else {
       for(i=0; i < nband; i++){
           if(!i){
             for(j=0; j < dimt; j++){
	         if(*(karray + j) >= *(kbands + i) && *(karray + j) <= *(kbands + i + 1))
	           *(fmask[i] + j) = 1.0;
                 else *(fmask[i] + j) = 0.0;
	     } 
	   }
	   else{
             for(j=0; j < dimt; j++){
	         if(*(karray + j) > *(kbands + i) && *(karray + j) <= *(kbands + i + 1))
	            *(fmask[i] + j) = 1.0;
                 else *(fmask[i] + j) = 0.0;
	     }	
	   }
       }   
    }
    
/* assign memory for output file pointers */

    filspec = (FILE **)calloc(nband, sizeof(FILE *));
    mem_er((filspec == NULL) ? 0 : 1, nband * sizeof(FILE *));
    
    filename = (char **)calloc(nband, sizeof(char *));
    mem_er((filename == NULL) ? 0 : 1, nband * sizeof(char *));
    for(i=0; i < nband; i++){
        *(filename + i) = calloc(MAXCHR, sizeof(char));
        mem_er((*(filename + i) == NULL) ? 0 : 1, MAXCHR * sizeof(char *));
    }

    printf("If field values are likely to be large, do you wish to restrict\r\n"
           "them before filtering, 'y' or 'n'?                             \n\n");

    scanf("\n");
    if((fres = getchar()) == 'y'){

       printf("Input the lower and upper bound value in the same units as the current field.\n\n");
       scanf("%f %f", &fmin, &fmax);

       if(fmin > fmax){
          printf("****ERROR****, fmin must be less than fmax, fmin = %e, fmax = %e\n\n", fmin, fmax);
          exit(1);
       }

       printf("What value to do you want to reset large values to?\n");
       scanf("%f", &fval);

    }
    
    
/* start processing fields */


    printf("***INFORMATON***, Computing spectral decomposition of fields....\r\n\n"
           "                                                                \r\n\n"
           " This may take some time ..............                           \n\n");  
	   
/* open files for writing filtered fields */

   for(i=0; i < nband; i++){
       sprintf(band_id, "_band%03d", i);
       strncpy(filename[i], SPECTRAL, MAXCHR);
       if(iext) strcpy(strstr(filename[i], EXTENSION), fext);
       strcat(filename[i], band_id);

       filspec[i] = open_file(filename[i], "w");

/* write header */

       write_header(gr, filspec[i]);

   }  
    
    nf = 0;   
    
/* read and process fields write out in simple binary format for now */
    
    while(nf <= frl){

         if(form != 4){

            if(!nf) {

               place1 = ftello(fdat); 
 
               if(read_field(ap, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;        

               if(fr1 == 1) ++nf;;
 
               place2 = ftello(fdat);
               chrnum = place2 - place1;

               if(fr1 > 1){

                  fseeko(fdat, (fr1-2)*chrnum, ORIGIN);
                  if(read_field(ap, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;

                  nf = fr1;           

               }

               nf += fri;

/* need to reset header strings for new grid */

            }

            else {

               fseeko(fdat, (fri-1)*chrnum, ORIGIN);
               if(read_field(ap, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;
          
               nf += fri;

            }

         }

         else {

           if(!nf) nf = fr1;
           ((NETCDF_INFO *)fdat)->iframe = nf - 1;
           if(read_field(ap, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;
           nf += fri;

         }

/* apply lower and upper bound field value if required */

         if(fres == 'y'){

            for(i=0; i < dim; i++){

                if(*(ap + i) <= fmin || *(ap + i) >= fmax) *(ap + i) = fval;


            }

         }
	 
	 if(icent){
	    savg = 0.0;
	    
	    for(i=0; i < dim; i++) savg += *(ap + i);
	    savg /= dim;
	    for(i=0; i < dim; i++) *(ap + i) -= (float)savg;
	 
	 }

	 
/* zero complex array */

	 memset(ac, 0, dimt * sizeof(complex));
	  
         for(i=0; i < gr->iy; i++){
	 	 
	     for(j=0; j < gr->ix; j++){
	         comp(*(ap + i * gr->ix + j), 0.0, ac + (i + ndy) * nx + j + ndx);	      
	     }
	     
	     if(irzo){
	        for(j=0; j < ndx; j++)
	           comp(*(ap + i * gr->ix + j + 1), 0.0, ac + (i + ndy) * nx - j + ndx - 1);
		   
		for(j=0; j < ndx2; j++)
		    comp(*(ap + i * gr->ix + gr->ix - j - 2), 0.0, ac + (i + ndy) * nx + j + ndx + gr->ix);
	     }
	     
	 }
	 
	 if(irzo){
	    for(i=0; i < nx; i++){
			     
	        for(j=0; j < ndy; j++)
	            comp((ac + (j + ndy + 1) * nx + i)->real, 0.0, ac + (ndy - j - 1) * nx + i);
		    
		for(j=0; j < ndy2; j++)
	            comp((ac + (ny - ndy2 - j - 2) * nx + i)->real, 0.0, ac + (ny - ndy2 + j) * nx + i);
		    
	    }
	    

	 }
	 
	 dfct2d(ac, nx, ny, nexpx, nexpy, 1, itrig, ifftyp);
	  
/* apply filters */

         for(i=0; i < nband; i++) {
	     for(j=0; j < dimt; j++) {
	        cmx(*(fmask[i] + j), *(ac + j), (acc + j));
             }
		 
	     dfct2d(acc, nx, ny, nexpx, nexpy, -1, itrig, ifftyp);
	      
	     itrig = 0;
	      
	     for(j=0; j < gr->iy; j++){
	         for(k=0; k < gr->ix; k++){
       	            *(ap + j * gr->ix + k) = (acc + (j + ndy) * nx + k + ndx)->real; 
		 }	      
	     }

/*printf("\n\n");
for(j=0; j < dim; j++)printf("%f\n", *(ap + j));
exit(1); */
	      
	     fprintf(filspec[i], "FRAME %6d\n", ifr+1);
             fwrite(ap, dim * sizeof(float), 1, filspec[i]);
             fprintf(filspec[i], "\n");	  
	 }
	  
	  
	  
	 ++ifr;

         if(nf > frl) break;
    }
    
    for(i=0; i < nband; i++){

      fseeko(filspec[i], (off_t)0, FSTART);
      fprintf(filspec[i], "%6d %6d %6d\n", gr->ix, gr->iy, ifr);

      close_file(filspec[i], filename[i]);

    }	
    
    if(ifftyp){
       free(dfreal);
       free(dfimag);
       free(dtrgx);
       free(dtrgy);
    }
    
    free(karray);
    free(alphamin);
    for(i=0; i < nband; i++) free(*(fmask + i));
    free(fmask);
    

    if(chrfld) free(chrfld);
    free(bands);
    free(kbands);
    free(ac);
    free(acc);
    free(ap);
    free(filspec);
    for(i=0; i < nband; i++) free(filename[i]);
    free(filename);

    return;
}
