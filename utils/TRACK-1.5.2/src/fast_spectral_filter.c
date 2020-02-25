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
#include "m_values.h"
#include "geo_values.h"

#define TOLGAUSS   2.0e-3
#define  NCHRB  30

/* function to perform spatial spectral filtering on an instantaneous field using the
   fast spectral transform, including conversion from non-gaussian grid to gaussian
   grid prior to applying the transform.                                              */

#ifdef  FFTLIB

#ifdef  NOUNDERSCORE

void setgpfa(double * , int * , int * );
void gpfa(double * , double * , double * , int * , int * , int * , int * , int *, int * );

#else

void setgpfa_(double * , int * , int * );
void gpfa_(double * , double * , double * , int * , int * , int * , int * , int *, int * );

#endif

#endif

int powi(int , int );
double legendre_node(int , int , double * );
int create_spherical_basis(GRID * , int , int );
int read_field(float * ,float * , float , FILE * , int , int , int , int , int );
void write_header(GRID * , FILE * );
double hoskins_filt(float , double );
GRID *gauss_grid(int );
int temptest(int * , int , int );
void del2(double *fn, int n, int ityp);

extern GRID *gr;
extern int form;   
extern complex *yt;
extern int iext;
extern char *fext;

void fast_spectral_filter(FILE *fdat, int fr1, int fri, int frl)
{
    int i, j, k, n;
    int nf=0, nx=0, nxx=0, nxmax=0, nln=0;
    int ncof=0;
    int ifr=0;
    int nfile=0;

    int nj[3], nj_n[3];
    int inc=1, lot=1, idir=1, indir=-1;
    int jump=0;
    int dimtrg=0, dimtrg_n=0;

    int fres='n';
    int gty=0;
    int ng='n';
    int nmax=0;

    int ivrdv=0, ider=0, ideriv=0, iderivt=0; 

    int ngauss=0;
    int dim=gr->ix * gr->iy, dimn=0;
    int nlat2=0, nlng=0;
    int ntrunc=0, ntrunc_n=0, ntr1=0;
    int nttmp1=0, nttmp2=0;

    int nband=0;
    int *bands=NULL;

    off_t place1, place2, chrnum=0;

    char **filename, filnm[MAXCHR], band_id[NCHRB];

    FILE **filspec=NULL, *fg=NULL;

    GRID *newg=NULL;

    float *ap=NULL, *app=NULL, *app2=NULL;
    float fmax, fval;
    float cut;

    float **ifilt=NULL;
    float *hf=NULL;
    
    float rerth2=EARTH_RADIUS_M*EARTH_RADIUS_M;

    double legnod=0.0;
    double *freal=NULL, *fimag=NULL;
    double *trigs=NULL, *trigs_n=NULL;
    double *legwt=NULL;
    double *coslat=NULL;
    double *fn=NULL;
    double pin=0.0, pin1=0.0, pin2=0.0;
    double nsq, cutf;

    complex *cfft1=NULL, *cfft2=NULL, *cfft22=NULL, *cff=NULL, *cff2=NULL;
    complex *cof=NULL, *ctmp=NULL;
    complex plf1, plf2, coft;

    printf("****WARNING****, data must currentely be defined on a global domain          \r\n"
           "                 with no wrap-round included and have no missing data values.\n\n");

    if(form == 4){

       if(((NETCDF_INFO *)fdat)->imiss){
          printf("****WARNING****, there are possible missing values in the data.\n\n");
       }

    }

/* normalization for foreward Legendre-FFT */
    pin1 = ((double)(gr->ix) * sqrt(2.0 * FPI2));
/* normalization for backward Legendre-FFT */
    pin2 = FPI2 / pin1;


    ap=(float *)calloc(dim, sizeof(float));
    mem_er((ap == NULL) ? 0: 1, dim * sizeof(float));

/* check if grid is gaussian */

    if(gr->ixfrm){
       printf("****ERROR****, grid is not uniform in X.\n\n");
       exit(1);
    }

    if(gr->igtyp > 1){
       printf("****ERROR****, not a global grid, can only handle global grids at the moment.\n\n");
       exit(1);
    }

    ngauss = gr->ix / 2;

    if(ngauss == gr->iy && !(ngauss % 2) && gr->iyfrm){

       legwt = (double *)calloc(ngauss, sizeof(double));
       mem_er((legwt == NULL) ? 0 : 1, ngauss*sizeof(double));

/* possible gaussian grid */
       nlat2 = gr->iy / 2;
       for(i=0; i < nlat2; i++){
           legnod = (FP_PI2 - acos(legendre_node(nlat2 - i, gr->iy, legwt + nlat2 - i - 1))) / FP_PI;
           if(fabs(*(gr->ygrid + nlat2 + i) - legnod) > TOLGAUSS){
              printf("****ERROR****, currently only gaussian grids are supported.\n\n");
              exit(1);   
           }
       }

       for(i=0; i < nlat2; i++) *(legwt + ngauss - i - 1) = *(legwt + i);
    }
    else {
       printf("****ERROR****, currently only gaussian grids are supported.\n\n");
       exit(1);
    }

/* check for gaussian grid type */

    nlng = 4;
    nttmp1 = 1;
    while(nlng != gr->ix){
         ++nttmp1;
         nlng = (2 * nttmp1 + 1);
         if(nlng % 2) nlng += 1;
         if(nlng % 4) nlng += 2;
         if(nlng == gr->ix){
            printf("****INFORMATION****, grid looks like a T%d linear gaussian grid.\n\n", nttmp1);
            break;
         }
         else if(nlng > gr->ix) {
            nttmp1 = -1;
            break;
         }
    }

    nlng = 4;
    nttmp2 = 1;
    while(nlng != gr->ix){
          ++nttmp2;
          nlng = (3 * nttmp2 + 1);
          if(nlng % 2) nlng += 1;
          if(nlng % 4) nlng += 2;
          if(nlng == gr->ix){
             printf("****INFORMATION****, grid looks like a T%d quadratic gaussian grid.\n\n", nttmp2);
             break;
          }
          else if(nlng > gr->ix){
            nttmp2 = -1;
            break;
       }
    }

    if(nttmp1 < 0 || nttmp2 < 0){
       printf("****ERROR****, grid incompatable with a gaussian Legendre transform.\n\n");
       exit(1);
    }

    printf("If the input field is relative vorticity or divergence are derived fields required, e.g.                  \r\n"
           "rotational winds or stream function; or divergent winds or velocity potential, 'y' for yes, 'n' for no.\n\n");
    scanf("\n");
    if(getchar() == 'y'){
       ideriv = 1;
       printf("Is the field vorticity or divergence, '0' for vorticity or '1' for divergence.\n\n");
       scanf("%d", &ivrdv);
       if(ivrdv < 0 || ivrdv > 1){
          printf("****ERROR****, incorrect specifier.\n\n");
          exit(1);
       }
       if(ivrdv == 0){ 
          printf("Do you want rotational winds or stream function, '0' for winds, '1' for stream function.\n\n");
	  scanf("%d", &iderivt);
       }
       else if(ivrdv == 1){
          printf("Do you want divergent winds or velocity potential, '0' for winds, '1' for velocity potential.\n\n");
	  scanf("%d", &iderivt);

       }
       if(iderivt < 0 || iderivt > 1){
          printf("****ERROR****, incorrect specifier.\n\n");
	  exit(1);
       }
       if(!iderivt) ider = 1;

    }

    printf("What spectral truncation is required?\n\n");
    scanf("%d", &ntrunc);

    ntr1 = ntrunc + 1;

    if(ntrunc < 0){
       printf("****ERROR****, truncation must be positive.\n\n");
       exit(1);
    }

    if(ntrunc > nttmp1){
       printf("****WARNING****, truncation may not be compatable with the grid.\n\n");
    }

    create_spherical_basis(gr, ntrunc, ider);
    free(gr->sico);

/* only store one half of triangular truncation since transform is real */

    ncof = ((ntrunc + 1) * (ntrunc + 2)) / 2;

    cof = (complex *)calloc(ncof, sizeof(complex));
    mem_er((cof == NULL) ? 0 : 1, ncof * sizeof(complex));

    cfft1 = (complex *)calloc(dim, sizeof(complex));
    mem_er((cfft1 == NULL) ? 0 : 1, dim * sizeof(complex));

/* setup output grid */

   printf("Do you want to output filtered data on a new grid? 'y' or 'n'\n\n");
   scanf("\n");
   if((ng = getchar()) == 'y') {
   
       if(gr->dleng) free(gr->dleng);

       printf("Read grid from file, input '0', or create a Gaussian grid, input '1'\n");
       scanf("%d", &gty);

       if(!gty){

          newg = (GRID *)malloc_initl(sizeof(GRID));
          mem_er((newg == NULL) ? 0 : 1, sizeof(GRID));

          printf("What is the filname for the grid data?\n");
          scanf("%s", filnm);

          fg = open_file(filnm, "r");

          fscanf(fg, "%d %d", &(newg->ix), &(newg->iy));

          newg->xgrid = (float *)calloc(newg->ix, sizeof(float));
          mem_er((newg->xgrid == NULL) ? 0 : 1, newg->ix * sizeof(float));

          newg->ygrid = (float *)calloc(newg->iy, sizeof(float));
          mem_er((newg->ygrid == NULL) ? 0 : 1, newg->iy * sizeof(float));

          for(i=0; i < newg->ix; i++) fscanf(fg, "%f", newg->xgrid + i);
          for(i=0; i < newg->iy; i++) fscanf(fg, "%f", newg->ygrid + i);

          close_file(fg, filnm);
          ntrunc_n = ntrunc;

       }

       else {
          printf("What is the spectral truncation (triangular) for the output gaussian grid?\n");
          scanf("%d", &ntrunc_n);

          if(ntrunc_n < 0){
             printf("****ERROR****, truncation must be positive, exiting.\n\n");
             exit(1);
          }
          if(ntrunc_n > nttmp1){
             printf("****ERROR****, truncation not compatable with the grid.\n\n");
             exit(1);
          }
          newg = gauss_grid(ntrunc_n);

          --(newg->ix);
       }

       create_spherical_basis(newg, ntrunc, ider);
       free(newg->sico);

       dimn = newg->ix * newg->iy;

       app=(float *)calloc(dimn, sizeof(float));
       mem_er((app == NULL) ? 0: 1, dimn * sizeof(float));

       cfft2 = (complex *)calloc(dimn, sizeof(complex));
       mem_er((cfft2 == NULL) ? 0 : 1, dimn * sizeof(complex));

   }

   else {newg = gr; dimn = dim; app = ap; ntrunc_n = ntrunc; cfft2 = cfft1;}
   
 /* calculate cosine of latitude and setup del weights */

   if(ideriv){
   
      fn = (double *)calloc(ntr1, sizeof(double));
      mem_er((fn == NULL) ? 0 : 1, ntr1*sizeof(double));
      del2(fn, ntr1, 1);
   
      if(ider){
         coslat = (double *)calloc(newg->iy, sizeof(double));
         mem_er((coslat == NULL) ? 0 : 1, newg->iy * sizeof(double));
         for(i=0; i < newg->iy; i++){
            *(coslat + i) = cos(*(newg->ygrid + i) * FP_PI);
         }
	 
      }
      
      if(!iderivt){
         cfft22 = (complex *)calloc(dimn, sizeof(complex));
         mem_er((cfft22 == NULL) ? 0 : 1, dimn * sizeof(complex));
      
         app2=(float *)calloc(dimn, sizeof(float));
         mem_er((app2 == NULL) ? 0: 1, dimn * sizeof(float));
      }

   }

/* setup filters */

   printf("How many filter bands are required?\n\n");
   scanf("%d", &nband);

   if(nband <= 0){
      printf("****ERROR****, number of bands must be greater than zero\n\n");
      exit(1);
   }

/* hoskins filter */

   hf = (float *)calloc(ntr1, sizeof(float));
   mem_er((hf == NULL) ? 0 : 1, ntr1 * sizeof(float));

   for(i=0; i<ntr1; i++) *(hf + i) = 1.0;

   printf("Do you want to apply further filtering in the form of the Hoskins filter, 'y' or 'n'\n");
   scanf("\n");

   if(getchar() == 'y'){

      printf("What is the cutoff constant value?\r\n"
             "Typicaly a value of 0.1 is used.  \n\n");

      scanf("%f", &cut);

      nsq = (float) ntrunc * (float) ntr1;
      nsq *= nsq;

      cutf = - log((double) cut) / nsq;

      for(i=0; i<ntr1; i++) *(hf + i) = hoskins_filt((float) i, cutf);

   }

   bands = (int *)calloc(nband+1, sizeof(int));
   mem_er((bands == NULL) ? 0 : 1, (nband+1) * sizeof(int));

   ifilt = (float **)calloc(nband, sizeof(float *));
   mem_er((ifilt == NULL) ? 0 : 1, nband * sizeof(float *));


   printf("Input the band boundaries in terms of total wavenumber.\n");

   for(i=0; i <= nband; i++){

      printf("Boundary %d = ", i+1);
      scanf("%d", bands+i);
      printf("\n");

      if(*(bands+i) < 0) {
         printf("****WARNING****, band boundary cannot be negative, resetting to zero.\n\n");
         *(bands + i) = 0;
      }
      else if(*(bands+i) > ntrunc) {
         printf("****WARNING****, maximum band boundary cannot be greater than truncation, resetting to truncation.\n\n");
         *(bands + i) = ntrunc;
      }

   }

   for(i=0; i < nband; i++){

      *(ifilt + i) = (float *)calloc(ntr1, sizeof(float));
      mem_er((*(ifilt + i) == NULL) ? 0 : 1, ntr1 * sizeof(float));
      memset(*(ifilt + i), 0, ntr1 * sizeof(float));

      for(j=*(bands + i)+((i)?1:0); j<= *(bands + i + 1); j++) *(ifilt[i] + j) = *(hf + j);

   }



/* restrict field values */

   printf("If field values are likely to be large, do you wish to restrict\r\n"
           "them before filtering, 'y' or 'n'?                             \n\n");

   scanf("\n");
   if((fres = getchar()) == 'y'){

      printf("Input the upper bound value in the same units as the current field.\n\n");
      scanf("%f", &fmax);

      printf("What value to do you want to reset large values to?\n");
      scanf("%f", &fval);

   }


/* setup arrays for FFT coefficients */

   nxmax = (gr->ix > newg->ix) ? gr->ix : newg->ix;


#ifdef  FFTLIB

/* test for valid number of longitude points for Temperton */

    if(!temptest(nj, gr->ix, 1)) exit(1);
    if(ng == 'y') {if(!temptest(nj_n, newg->ix, 1)) exit(1);}
    else {nj_n[0] = nj[0]; nj_n[1] = nj[1]; nj_n[2] = nj[2];}

    dimtrg = powi(2, nj[0]) + powi(3, nj[1]) + powi(5, nj[2]);
    dimtrg *= 2;

    freal = (double *)calloc(nxmax, sizeof(double));
    mem_er((freal == NULL) ? 0 : 1, nmax*sizeof(double));

    fimag = (double *)calloc(nxmax, sizeof(double));
    mem_er((fimag == NULL) ? 0 : 1, nmax*sizeof(double));

    trigs = (double *)calloc(dimtrg, sizeof(double));
    mem_er((trigs == NULL) ? 0 : 1, dimtrg*sizeof(double));

    if(ng == 'y'){
       dimtrg_n = powi(2, nj_n[0]) + powi(3, nj_n[1]) + powi(5, nj_n[2]);
       dimtrg_n *= 2;

       trigs_n = (double *)calloc(dimtrg_n, sizeof(double));
       mem_er((trigs == NULL) ? 0 : 1, dimtrg_n*sizeof(double));
    }

#ifdef  NOUNDERSCORE

    setgpfa(trigs, &(gr->ix), nj);
    if(ng == 'y')setgpfa(trigs_n, &(newg->ix), nj_n);
#else

    setgpfa_(trigs, &(gr->ix), nj);
    if(ng == 'y')setgpfa_(trigs_n, &(newg->ix), nj_n);

#endif

    if(! (ng == 'y')){
       dimtrg_n = dimtrg;
       trigs_n = trigs;
    }

#else
    printf("****WARNING****, must build and link fft library to use temperton FFT.\n\n");
    exit(1);
 
#endif

    nx = gr->ix;
    nxx = newg->ix;

/* assign memory for output file pointers */

    if(!ideriv) nfile = nband;
    else {
       if(!iderivt) nfile = 2 * nband;
       else nfile = nband;
    }
    
    filspec = (FILE **)calloc(nfile, sizeof(FILE *));
    mem_er((filspec == NULL) ? 0 : 1, nfile * sizeof(FILE *));
    filename = (char **)calloc(nfile, sizeof(char *));
    mem_er((filename == NULL) ? 0 : 1, nfile * sizeof(char *));
    for(i=0; i < nfile; i++){
        *(filename + i) = calloc(MAXCHR, sizeof(char));
        mem_er((*(filename + i) == NULL) ? 0 : 1, MAXCHR * sizeof(char *));
    }

/* open files for writing filtered fields */

   for(i=0; i < nband; i++){
       sprintf(band_id, "_band%03d", i);
       strncpy(filename[i], SPECTRAL, MAXCHR);
       if(iext) strcpy(strstr(filename[i], EXTENSION), fext);
       strcat(filename[i], band_id);
       
       if(ideriv){
          if(!ivrdv){
	     if(!iderivt){
	        strncpy(filename[i + nband], filename[i], MAXCHR);
	        strcat(filename[i], ".Urot");
		strcat(filename[i + nband], ".Vrot");
	     }
	     else {
	        strcat(filename[i], ".strf");
	     }
	  }
	  else {
	     if(!iderivt){
	        strncpy(filename[i + nband], filename[i], MAXCHR);
	        strcat(filename[i], ".Udiv");
		strcat(filename[i + nband], ".Vdiv");
	     }
	     else {
	        strcat(filename[i], ".velpot");
	     }
	  }
       
       }

       filspec[i] = open_file(filename[i], "w");

/* write header */

       write_header(newg, filspec[i]);
       
       if(ideriv && !iderivt){
         filspec[i + nband] = open_file(filename[i + nband], "w");
	 write_header(newg, filspec[i + nband]);
       }

   }

    nf = 0;

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

/* apply upper bound field value if required */

         if(fres == 'y'){

            for(i=0; i < dim; i++){

                if(*(ap + i) > fmax) *(ap + i) = fval;

            }

         }

/* perform FFT */


#ifdef  FFTLIB

         for(i=0; i < gr->iy; i++){
             cff = cfft1 + i * gr->ix;
             memset(fimag, 0, nxmax*sizeof(double));
             memset(freal, 0, nxmax*sizeof(double));
             for(j=0; j < gr->ix; j++) *(freal + j) = *(ap + i * gr->ix + j);

#ifdef  NOUNDERSCORE

             gpfa(freal, fimag, trigs, &inc, &jump, &(gr->ix), &lot, &idir, nj);

#else

             gpfa_(freal, fimag, trigs, &inc, &jump, &(gr->ix), &lot, &idir, nj);

#endif

             for(j=0; j < gr->ix; j++) comp(*(freal + j), -*(fimag + j), (cff + j));


         }

#endif


/* foreward legendre transform */

         ctmp = cof;
         nln = 0;
         for(i=0; i <= ntrunc; i++){
             pin = pin2 * pow(-1.0, (double) i);
             for(j=i; j <= ntrunc; j++){
                 ctmp->real = ctmp->imag = 0.0;
                 for(k=0; k < gr->iy; k++){       /* Gauss-Legendre */
                     cmx(*(gr->aleng + k * gr->nleng + nln), *(cfft1 + k * nx + i), &plf1);
                     cmx(*(legwt + k), plf1,  &plf2);
                     cadd(plf2, *ctmp, ctmp);
                 }

                 cmx(pin, *ctmp, ctmp);

                 ++nln;
                 ++ctmp;

             }

         }

/* apply del2 operator */
         if(ideriv){
	    ctmp = cof;
            for(j=0; j <= ntrunc; j++){
               for(k=j; k <= ntrunc; k++) {
	           cmx(*(fn + k), *ctmp, ctmp);
                   ++ctmp;
	       }  
            }
	    if(iderivt) cof->real = cof->imag = 0.0;
	 }

/* loop over filters */

         for(n=0; n < nband; n++){

/* inverse legendre transform */

            if(!ideriv || iderivt){

/* normal transform for non-derived fields or del2 operated fields */

               for(i=0; i < newg->iy; i++) {
                   cff = cfft2 + i * nxx;
                   nln = 0;
                   ctmp = cof;
                   memset(cff, 0, nxx * sizeof(complex));
                   for(j=0; j <= ntrunc; j++){
                       pin = pin1 * pow(-1.0, (double) j);
                       for(k=j; k <= ntrunc; k++){
                           cmx(*(ifilt[n] + k), *(ctmp + nln), &coft);
                           cmx(*(newg->aleng + i * newg->nleng + nln), coft, &plf1);
                           cadd(plf1, *(cff + j), cff + j);
                           ++nln;
                       }
                       cmx(pin, *(cff + j), (cff + j));

                   }

               }
	       
	    }
	    
	    else if(ideriv && !iderivt){

/* transform for derived fields, rotational or divergent winds */
	    
               for(i=0; i < newg->iy; i++) {
                   cff = cfft2 + i * nxx;
		   cff2 = cfft22 + i * nxx;
                   nln = 0;
                   ctmp = cof;
                   memset(cff, 0, nxx * sizeof(complex));
		   memset(cff2, 0, nxx * sizeof(complex));
                   for(j=0; j <= ntrunc; j++){
                       pin = pin1 * pow(-1.0, (double) j);
                       for(k=j; k <= ntrunc; k++){
                           cmx(*(ifilt[n] + k), *(ctmp + nln), &coft);
			   if(!ivrdv){
                              cmx(-(*(newg->dleng + i * newg->nleng + nln)), coft, &plf1);
			      cadd(plf1, *(cff + j), cff + j);
			      comp(-(coft.imag), (coft.real), &plf1);
			      comp((plf1.real), (plf1.imag), &coft);      
			      cmx(((double) j) * *(newg->aleng + i * newg->nleng + nln), coft, &plf1);
			      cadd(plf1, *(cff2 + j), cff2 + j);		      
			   }
			   else{
			      cmx(*(newg->dleng + i * newg->nleng + nln), coft, &plf1);
			      cadd(plf1, *(cff2 + j), cff2 + j);
			      comp(-(coft.imag), (coft.real), &plf1);
			      comp((plf1.real), (plf1.imag), &coft);
			      cmx(((double) j) * *(newg->aleng + i * newg->nleng + nln), coft, &plf1);
			      cadd(plf1, *(cff + j), cff + j);
			   }
                           ++nln;
                       }
                       cmx(pin, *(cff + j), (cff + j));
                       cmx(pin, *(cff2 + j), (cff2 + j));
                   }

               }	     


            }
	    
	    memset(ap, 0, dim * sizeof(float));

#ifdef  FFTLIB

            for(i=0; i < newg->iy; i++){
                cff = cfft2 + i * newg->ix;
                memset(fimag, 0, nxmax*sizeof(double));
                memset(freal, 0, nxmax*sizeof(double));
                *freal = cff->real;
                *fimag = cff->imag;
                for(j=1; j <= newg->ix/2; j++) {
                    *(freal + j) =  (cff + j)->real;
                    *(fimag + j) =  -(cff + j)->imag;
                    *(freal + newg->ix - j) = (cff + j)->real;
                    *(fimag + newg->ix - j) = (cff + j)->imag;
                }

#ifdef  NOUNDERSCORE

                gpfa(freal, fimag, trigs_n, &inc, &jump, &(newg->ix), &lot, &indir, nj_n);

#else

                gpfa_(freal, fimag, trigs_n, &inc, &jump, &(newg->ix), &lot, &indir, nj_n);

#endif


/* for inverse transform need to normalise by dividing by number of longitudes */

                for(j=0; j < newg->ix; j++) *(app + i * newg->ix + j) = *(freal + j) / (double)(gr->ix);
		
/* if computing rotational or divergent winds also do the V transform */
		
		if(ideriv && !iderivt){
		    cff = cfft22 + i * newg->ix;
                    memset(fimag, 0, nxmax*sizeof(double));
                    memset(freal, 0, nxmax*sizeof(double));
                    *freal = cff->real;
                    *fimag = cff->imag;
                    for(j=1; j <= newg->ix/2; j++) {
                        *(freal + j) =  (cff + j)->real;
                        *(fimag + j) =  -(cff + j)->imag;
                        *(freal + newg->ix - j) = (cff + j)->real;
                        *(fimag + newg->ix - j) = (cff + j)->imag;
                    }

#ifdef  NOUNDERSCORE

                    gpfa(freal, fimag, trigs_n, &inc, &jump, &(newg->ix), &lot, &indir, nj_n);

#else

                    gpfa_(freal, fimag, trigs_n, &inc, &jump, &(newg->ix), &lot, &indir, nj_n);

#endif


/* for inverse transform need to normalise by dividing by number of longitudes */

                    for(j=0; j < newg->ix; j++) *(app2 + i * newg->ix + j) = *(freal + j) / (double)(gr->ix);
		}

            }

#endif

            if(ideriv){
               if(iderivt){
	          for(i=0; i < dimn; i++) *(app + i) *= rerth2;
	       }
	       else {
	          for(i=0; i < newg->iy; i++) {
		      for(j=0; j < newg->ix; j++){
		         *(app + i * newg->ix + j) *= EARTH_RADIUS_M / *(coslat + i);
		         *(app2 + i * newg->ix + j) *= EARTH_RADIUS_M / *(coslat + i);
		     }
		  }
		  
                  fprintf(filspec[n + nband], "FRAME %6d\n", ifr+1);
                  fwrite(app2, dimn * sizeof(float), 1, filspec[n + nband]);
                  fprintf(filspec[n + nband], "\n");		  
		  
	       }
	    }
	    
            fprintf(filspec[n], "FRAME %6d\n", ifr+1);
            fwrite(app, dimn * sizeof(float), 1, filspec[n]);
            fprintf(filspec[n], "\n");

         }

         ++ifr;

         if(nf > frl) break;

    }

    for(i=0; i < nfile; i++){

       fseeko(filspec[i], (off_t)0, FSTART);
       fprintf(filspec[i], "%6d %6d %6d\n", newg->ix, newg->iy, ifr);

       close_file(filspec[i], filename[i]);
       
       free(filename[i]);

    }
    
    for(i=0; i < nfile; i++) free(*(filename + i));
	
    free(filspec);
    free(filename);
    
    
    if(coslat) free(coslat);
    if(fn) free(fn);
    
    free(bands);
    for(i=0; i < nband; i++) free(*(ifilt+i));
    
    

    free(freal); free(fimag); free(trigs);
    if(ng == 'y'){ free(trigs_n);}
    free(legwt);
    free(cfft1);
    free(cof);

    if(ng == 'y'){
       free(cfft2);
       free(newg->aleng);
       if(newg->dleng) free(newg->dleng);
       free(newg->sico);
       free(newg->xgrid);
       free(newg->ygrid);
       free(newg);

       free(app);

    }
    
    if(cfft22) free(cfft22);
    
    if(app2) free(app2);

    free(gr->aleng);
    if(gr->dleng) free(gr->dleng);
    free(gr->sico);

    free(ap);
    free(hf);

    return;
}

/* checks for Temperton */

int temptest(int *nj, int nn, int ip)
{
    int i;
    int kk=0, ifac=2;
    int itype = 1;

    for(i=1; i <= 3; i++){
       kk = 0;
       while(!(nn % ifac)){
          ++kk;
          nn /= ifac;
       }
       nj[i-1] = kk;
       ifac += i;
    }

    if(nn != 1){
       if(ip) printf("****WARNING****, value of number of grid points or longitudes is not valid for Temperton.\n\n");
       itype = 0;    
    }

    return itype;
}



void del2(double *fn, int n, int ityp)
{

   int i=0;
   
/* itype = 0 for del2 and itype = 1 for del-2 */

   if(!ityp){
      for(i=0; i < n; i++){
         *(fn + i) = -(double)(i * (i + + 1)); 
      }
   }
   else {
      *fn = -1.0; 
      for(i=1; i < n; i++){
         *(fn + i) = -1.0 / (double)(i * (i + 1));
      }
   }
   
   return;

}
