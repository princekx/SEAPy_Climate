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

/* function to perform a spatial spectral filter on an instanaeous field */

GRID *gauss_grid(int );
int create_spherical_basis(GRID * , int , int );
void write_header(GRID * , FILE * );
float sp_harmonic(double * , float * , GRID * , int , int , float * , float * );
double *least_sq_nm(float * , int , int , int , int , int );
double hoskins_filt(float , double );

int read_field(float * ,float * , float , FILE * , int , int , int , int , int );
int powi(int , int );

extern GRID *gr;
extern char *chrfld;

extern char *fext;

extern int iext;

extern int form;


void spectral_filter(FILE *fdat, int fr1, int fri, int frl)

{


    int i, j=0, k;
    int gty;
    int nband, ntr=0, ntr1, nleng=0;
    int dim=gr->ix * gr->iy, dimn;
    int fram=0;
    int immp=0;
    int ntr_new=0;

    int nx, nf;
    int fls=1;
    int ng='n';
    int imk=0;
    int fres='n';
    int bs=0;
    int nmk=0;
    int ider=0;

    off_t place1, place2, chrnum=0;
    long int memuse=0;
 
    int *bands=NULL, *imask=NULL;

    float *filf=NULL, *fsub=NULL;
    float *ap=NULL;
    float rms_error=0.0;
    float cut;
    float fmin, fmax, fval;

    double *cof=NULL;
    double nsq, cutf;

    float **ifilt=NULL, *nfilt=NULL;
    float *hf=NULL;


    GRID *newg=NULL;

    FILE *fg=NULL;
    FILE **filspec=NULL;

    char filename[MAXCHR], band_id[NCHRB];

    printf("****WARNING****, data must currentely be defined on a global domain\r\n"
           "                 and have no missing data values.                  \n\n");

    if(form == 4){

       if(((NETCDF_INFO *)fdat)->imiss){
          printf("****WARNING****, there are possible missing values in the data.\n\n");
       }

    }

    ap=(float *)calloc(dim, sizeof(float));
    mem_er((ap == NULL) ? 0: 1, dim * sizeof(float));

    printf("Do you want a full spectral decomposition, or just band subtraction.\r\n"
           "Band subtraction may be quicker for removing certain spatial scales.\r\n"
           "Input '0' for full decomposition or '1' for band subtraction.       \n\n");

    scanf("%d", &bs);

    printf("What spectral truncation is required (Triangular)?\n");
    scanf("%d", &ntr);

    if(ntr < 0){

       printf("****ERROR****, truncation must be positive, exiting.\n\n");
       exit(1);

    }


    ntr1 = ntr + 1;

    memuse = powi(ntr1,2) * (powi(ntr1, 2) + 1) * sizeof(double) / 2;

    printf("****INFORMATION****, this spectral truncation requires ~%6.1fMB of memory\r\n"
           "                     to perform the least squares spectral computation.  \r\n"
           "                     If your machine has less than or close to this      \r\n"
           "                     amount of memory, consider using the option         \r\n"
           "                     of memory mapped files.                             \n\n",
           (float)memuse / 1.0e+6);

    printf("Do you want to use memory mapped files, needs a working mmap function, input '0' for no and '1' for yes.\n\n");
    scanf("%d", &immp);


    for(i=0; i <= ntr; i++)nleng += (ntr1 - i);

/* create new grid */

    if(!bs){

       printf("Do you want to output filtered fields on a new Gaussian grid, \r\n"
              "or a new grid read from file, 'y' or 'n'. Otherwize existing  \r\n"
              "grid is used\n\n");

       scanf("\n");
       if((ng=getchar()) == 'y'){

          printf("Read grid from file, input '0', or create a Gaussian grid, input '1'\n");
          scanf("%d", &gty);

          if(!gty){

             newg = (GRID *)malloc_initl(sizeof(GRID));
             mem_er((newg == NULL) ? 0 : 1, sizeof(GRID));

             printf("What is the filname for the grid data?\n");
             scanf("%s", filename);

             fg = open_file(filename, "r");

             fscanf(fg, "%d %d", &(newg->ix), &(newg->iy));

             newg->xgrid = (float *)calloc(newg->ix, sizeof(float));
             mem_er((newg->xgrid == NULL) ? 0 : 1, newg->ix * sizeof(float));

             newg->ygrid = (float *)calloc(newg->iy, sizeof(float));
             mem_er((newg->ygrid == NULL) ? 0 : 1, newg->iy * sizeof(float));

             for(i=0; i < newg->ix; i++) fscanf(fg, "%f", newg->xgrid + i);
             for(i=0; i < newg->iy; i++) fscanf(fg, "%f", newg->ygrid + i);


             close_file(fg, filename);


          }

          else {
             printf("Specify the truncation for the new grid\n\n");
             scanf("%d", &ntr_new);
             newg = gauss_grid(ntr_new);
             printf("Do you want to include the longitude wrapround? 'y' or 'n' \n\n");
             scanf("\n");
             if(getchar() == 'n') --(newg->ix);
             
          }

          dimn = newg->ix * newg->iy;

       }

       else {

           printf("****INFORMATION****, using existing grid for filtered output.\n\n");

           newg = gr;

           dimn = dim;

       }

    }

    else {

      ng = 'n';
      newg = gr;
      dimn = dim;

      fsub = (float *)calloc(dimn, sizeof(float));
      mem_er((fsub == NULL) ? 0 : 1, dimn * sizeof(float));


    }


/* assign memory for filtered field */


    filf = (float *)calloc(dimn, sizeof(float));
    mem_er((filf == NULL) ? 0 : 1, dimn * sizeof(float));



/* compute spherical harmonics for current grid. */

    nx = create_spherical_basis(gr, ntr, ider);

    while(nx*gr->iy < (ntr + 1) * (ntr + 1)){ 

       printf("****WARNING****, this choice of the truncation will result in   \r\n"
              "                 underdetermined least squares problem to solve,\r\n"
              "                 due to the number of coefficients being larger \r\n"
              "                 than the number of data points, choose a better\r\n"
              "                 truncation with less coefficients.             \r\n"
              "                 the current number of coefficients is %d and   \r\n"
              "                 the number of data points is %d.\n\n", (ntr + 1) * (ntr + 1), nx*gr->iy);

       printf("What spectral truncation is required (Triangular)?\n");
       scanf("%d", &ntr);

       free(gr->aleng);
       free(gr->sico);
       nx = create_spherical_basis(gr, ntr, ider);
     

    }

/* create spherical basis for new grid */

    if(ng == 'y') create_spherical_basis(newg, ntr, ider);

    printf("How many spectral bands are required (with respect to total wavenumber)?\n");
    scanf("%d", &nband);

    if(nband <= 0){
       printf("****ERROR****, number of bands must be greater than zero\n\n");
       exit(1);
    }


    bands = (int *)calloc(nband+1, sizeof(int));
    mem_er((bands == NULL) ? 0 : 1, (nband+1) * sizeof(int));

    imask = (int *)calloc(nband, sizeof(int));
    mem_er((imask == NULL) ? 0 : 1, nband * sizeof(int));

    printf("Input the band boundaries in terms of total wavenumber.\n");

    for(i=0; i <= nband; i++){

        printf("Boundary %d = ", i+1);
        scanf("%d", bands+i);
        printf("\n");

        if(*(bands+i) < 0) *(bands + i) = 0;
        else if(*(bands+i) > ntr) *(bands + i) = ntr;

    }

    printf("****INFORMATION****, the following bands will be used.\n\n");
    for(i=0; i <= nband; i++) printf("%d ", *(bands + i));
    printf("\n\n");

    if(nband > 1){

       printf("Indicate by typing a 0 any bands you wish to mask.\r\n"
              "For masking only one filtered field is output.    \n\n");

       for(i=0; i < nband; i++){

          printf("Mask for band %d = ", i+1);
          scanf("%d", imask+i);
          printf("\n");
          if(! *(imask+i)) imk = 1;
          else ++nmk;

       }

    }

    else if (nband == 1) {imk = 1; *imask = 1;}

    if(imk && nmk > 1){
       printf("****ERROR****, when masking only one band should be left, exiting.\n\n");
       exit(1);
    }

/* assign coefficient filters */

    nfilt = (float *)calloc(ntr1, sizeof(float));
    mem_er((nfilt == NULL) ? 0 : 1, ntr1 * sizeof(float));

    for(i=0; i<ntr1; i++) *(nfilt + i) = 1.0;

    if(imk){

       ifilt = (float **)calloc(1, sizeof(float *));
       mem_er((ifilt == NULL) ? 0 : 1, sizeof(float *));
       *ifilt = (float *)calloc(ntr1, sizeof(float));
       mem_er((*ifilt == NULL) ? 0 : 1, ntr1 * sizeof(float));

       for(i=0; i < ntr1; i++) *(ifilt[0] + i) = 0.0;

       for(i=0; i < nband; i++){

           if(*(imask + i)){

              for(j=*(bands + i)+((i)?1:0); j<= *(bands + i + 1); j++) *(ifilt[0] + j) = 1.0;
           }


       }

       nband = 1;

       printf("Do you want to apply further filtering in the form of the Hoskins filter, 'y' or 'n'\n");
       scanf("\n");

       if(getchar() == 'y'){

          printf("What is the cutoff constant value?\r\n"
                 "Typicaly a value of 0.1 is used.  \n\n");

          scanf("%f", &cut);

          nsq = (float) ntr * (float) ntr1;
          nsq *= nsq;

          cutf = - log((double) cut) / nsq;

          for(i=0; i < ntr1; i++) *(ifilt[0] + i) *= hoskins_filt((float) i, cutf);

       }

    }
    else {

       ifilt = (float **)calloc(nband, sizeof(float *));
       mem_er((ifilt == NULL) ? 0 : 1, nband * sizeof(float *));

       hf = (float *)calloc(ntr1, sizeof(float));
       mem_er((hf == NULL) ? 0 : 1, ntr1 * sizeof(float));

       for(i=0; i<ntr1; i++) *(hf + i) = 1.0;

       printf("Do you want to apply further filtering in the form of the Hoskins filter, 'y' or 'n'\n");
       scanf("\n");

       if(getchar() == 'y'){

          printf("What is the cutoff constant value?\r\n"
                 "Typicaly a value of 0.1 is used.  \n\n");

          scanf("%f", &cut);

          nsq = (float) ntr * (float) ntr1;
          nsq *= nsq;

          cutf = - log((double) cut) / nsq;

          for(i=0; i<ntr1; i++) *(hf + i) = hoskins_filt((float) i, cutf);

       }


       for(i=0; i < nband; i++){

           ifilt[i] = (float *)calloc(ntr1, sizeof(float));
           mem_er((ifilt[i] == NULL) ? 0 : 1, ntr1 * sizeof(float));

           for(j=0; j < ntr1; j++) *(ifilt[i] + j) = 0.0;

           printf("Use filter for this band %d, 'y' or 'n'\n\n", i+1);
           scanf("\n");
           if(getchar() == 'y'){

              for(j=*(bands + i)+((i)?1:0); j<= *(bands + i + 1); j++) *(ifilt[i] + j) = *(hf + j);

           }

           else {

              for(j=*(bands + i)+((i)?1:0); j<= *(bands + i + 1); j++) *(ifilt[i] + j) = 1.0;

           }

       }

       free(hf);

    }

/* assign memory for output file pointers */

    filspec = (FILE **)calloc(nband, sizeof(FILE *));
    mem_er((filspec == NULL) ? 0 : 1, nband * sizeof(FILE *));

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
           " This may take some time ..............              \n\n");

/* open files for writing filtered fields */

   for(i=0; i < nband; i++){
       sprintf(band_id, "_band%03d", i);
       strncpy(filename, SPECTRAL, MAXCHR);
       if(iext) strcpy(strstr(filename, EXTENSION), fext);
       strcat(filename, band_id);

       filspec[i] = open_file(filename, "w");

/* write header */

       write_header(newg, filspec[i]);

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

/* perform spectral decomposition */

         cof=least_sq_nm(ap, ntr, fls, nx, immp, memuse);
         if(fls == 1) fls = 0;

         if(!cof){
            printf("****ERROR****, failure of least squares, no fit for this field, continuing.\n\n");
            continue;

         }

/* compute RMS error */

         rms_error = sp_harmonic(cof, NULL, gr, ntr, nx, ap, nfilt);
         printf("****INFORMATION****, RMS error for harmonic LS = %e\n", rms_error);

/* write out new field */

         ++fram;

         for(i=0; i<nband; i++){

/* apply filtering */

            sp_harmonic(cof, filf, newg, ntr, newg->ix, NULL, ifilt[i]);

            fprintf(filspec[i], "FRAME %6d\n", fram);

            if(!bs)
              fwrite(filf, dimn * sizeof(float), 1, filspec[i]);
            else{
              for(k=0; k < dimn; k++) *(fsub+k) = *(ap+k) - *(filf + k);

              fwrite(fsub, dimn * sizeof(float), 1, filspec[i]);
            }
            fprintf(filspec[i], "\n");

         }

         if(nf > frl) break;

    }

    printf("\n\n");

    for(i=0; i < nband; i++){

       sprintf(band_id, "_band%03d", i);
       strncpy(filename, SPECTRAL, MAXCHR);
       if(iext) strcpy(strstr(filename, EXTENSION), fext);
       strcat(filename, band_id);

       fseeko(filspec[i], (off_t)0, FSTART);
       fprintf(filspec[i], "%6d %6d %6d\n", newg->ix, newg->iy, fram);

       close_file(filspec[i], filename);

    }


    least_sq_nm(NULL, 0, -1, 0, 0, 0);

    free(bands);
    free(imask);

    for(i=0; i < nband; i++) free(*(ifilt+i));
    free(ifilt);
    free(nfilt);

    if(ng == 'y'){
       free(newg->aleng);
       free(newg->sico);
       free(newg->xgrid);
       free(newg->ygrid);
       free(newg);
    }

    free(gr->aleng);
    free(gr->sico);

    if(chrfld) free(chrfld);


    free(filf);
    free(fsub);
    free(ap);


    free(filspec);

    return;

}
