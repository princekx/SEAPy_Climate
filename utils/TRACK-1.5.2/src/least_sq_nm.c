#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#if defined(AIX)
#define _XOPEN_SOURCE 500
#endif
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "mem_er.h"
#include "grid.h"
#include "file_handle.h"
#include "files_out.h"


#define TOLDIAG 1.0e-10

#define FSPPATH  Add(USER,SPPATH,)

extern GRID *gr;


int cholesky(double ** , int );


double *least_sq_nm(float *ap, int ntr, int im, int nx, int immp, int memuse)

{

    int i, j, m, n, k;
    int ii;

    int ntr1 = ntr + 1;
    int ncof=ntr1 * ntr1;
    int dim = nx * gr->iy;

    int np;
    int fid=-1;

    int ier;
    int iread=0;

    long int ntms=0;

    static int *ind=NULL;

    float *atp=NULL;

    static double **aa=NULL;
    static double *vec=NULL;
    static double *bb=NULL;
    static double *amp=NULL;

/*    double *at=NULL; */

    double sum;

    double *leng=NULL, *cleng=NULL;
    complex *lng=NULL;

    FILE *fchol=NULL;

    char filchol[MAXCHR], filnam[MAXCHR];

    struct stat stt;

/* allocate memory */

    if(im > 0){

      if(ncof <= 0 || ncof > dim){

          printf("***ERROR***, insufficient or too many coefficients to fit.\n\n");
          return NULL;

       }

/* only memory map the least squares array at the moment */

       aa = (double **)calloc(ncof, sizeof(double *));
       mem_er((aa == NULL) ? 0 : 1, ncof * sizeof(double *));

       if(!immp){

         for(i=0; i < ncof; i++){

             *(aa + i) = (double *)calloc(i+1, sizeof(double));
             mem_er((*(aa+i) == NULL) ? 0 : 1, (i+1) * sizeof(double));

             for(j=0; j < i+1; j++) *(aa[i] + j) = 0.0;

         }

       }

       vec = (double *)calloc(ncof, sizeof(double));
       mem_er((vec == NULL) ? 0 : 1, ncof * sizeof(double));

       bb = (double *)calloc(ncof, sizeof(double));
       mem_er((bb == NULL) ? 0 : 1, ncof * sizeof(double));

       ind = (int *)calloc(ncof, sizeof(int));
       mem_er((ind == NULL) ? 0 : 1, ncof * sizeof(int));


       atp = ap;

/* populate the normal equations matrix  */

       printf("Do you want to load an existing cholesky decomposition of the          \r\n"
              "normal equations, 'y' or 'n'. This data must be correct for the        \r\n"
              "current grid and spectral truncation, this is the users responsibility.\r\n"
              "The data must also reside in the directory:\n %s\n", FSPPATH);

       scanf("\n");
       if(getchar() == 'y'){

          printf("What is the filname to read?\n\n"); 
          scanf("%s", filnam); 
          strcpy(filchol, FSPPATH);
          strcat(filchol, filnam);

          if(!immp){

             fchol = open_file(filchol, "r");

             ntms = 0;

             for(i=0; i<ncof; i++)

                 ntms += fread(aa[i], (i+1)*sizeof(double), 1, fchol);

             close_file(fchol, filchol);

             if(ntms != ncof){

                printf("****ERROR****, incorrect number of data read for cholesky decomposition\n\n");
                exit(1);


             }

          }

          else {

             if(stat(filchol, &stt) == -1){

                printf("****ERROR****, returning statistics on file %s\n\n", filchol);
                exit(1);

             }

             if(stt.st_size != memuse){

                printf("****ERROR****, file %s does not contain enough data for the chosen spectral truncation.\n\n", filchol);
                exit(1);

             }

             if((fid = open(filchol, O_RDONLY)) < 0){

                printf("****ERROR****, file %s does not exist or cannot be read.\n\n", filchol);
                exit(1);

             }

             if((amp = (double *)mmap(0, memuse, PROT_READ, MAP_SHARED, fid, 0)) == MAP_FAILED){

                printf("****ERROR****, memory mapped file failed, see mmap documentation for your system.\n\n");
                exit(1);

             }

          }

          iread = 1;

       }

       else if(immp){

          printf("What filname should be used for memory mapping?\n\n"); 
          scanf("%s", filnam); 
          strcpy(filchol, FSPPATH);
          strcat(filchol, filnam);

/* create file for memory mapping */

          fchol = open_file(filchol, "w");

          *aa = (double *)calloc(ncof, sizeof(double));
          mem_er((*aa == NULL) ? 0 : 1, ncof * sizeof(double));    

          for(i=0; i<ncof; i++)
              fwrite(*aa, (i+1)*sizeof(double), 1, fchol);
          close_file(fchol, filchol);

          free(*aa);


          if((fid = open(filchol, O_RDWR)) < 0){

             printf("****ERROR****, file %s does not exist or cannot be opened.\n\n", filchol);
             exit(1);

          }



          if((amp = (double *)mmap(0, memuse, (PROT_READ | PROT_WRITE), MAP_SHARED, fid, 0)) == MAP_FAILED){

             printf("****ERROR****, memory mapped file failed, see mmap documentation for your system.\n\n");
             exit(1);

          }

       }

       if(immp){

          ii = 0;

          for(i=0; i < ncof; i++){

              aa[i] = amp + ii;
              ii += i + 1;

          }

       }

       if(!iread){

          printf("****INFORMATION*****, creating normal equations matrix for least squares spectral fit\n\n");


          for(i=0; i < gr->iy; i++){

              leng = gr->aleng + i * gr->nleng;



              for(j=0; j < nx; j++){

                  lng = gr->sico + j * ntr1;

                  np = 0;

                  k = 0;

                  atp = ap + i * gr->ix + j;

                  for(m=0; m < ntr1; m++){

                      cleng = leng + np;


                      for(n=m; n < ntr1; n++){

                          vec[k++] = *cleng * lng->real;
 
                          if(m > 0) vec[k++] = *cleng * lng->imag; 
  
                          ++cleng;


                      }

                      ++lng;

                      np += (ntr1 - m);


                  }


                  for(m=0; m < ncof; m++){

                      for(n=0; n <=m; n++)

                          *(aa[m] + n) += vec[m] * vec[n];
 

                  }


              }

          }

          printf("****INFORMATION****, completed populating normal equations matrix.\n\n");





/*for(i=0; i<ncof; i++){

for(j=0; j<i+1; j++) printf("%8.4f ", *(aa[i] + j));
printf("\n");



}

printf("\n\n"); */

          ier = cholesky(aa, ncof);


/*for(i=0; i<ncof; i++){

for(j=0; j<i+1; j++) printf("%8.4f ", *(aa[i] + j));
printf("\n");



}

printf("\n\n"); */

         if(!ier){

            printf("****INFORMATION****, normal equations matrix populated and cholesky decomposition completed sucessfully.\n\n");

            if(!immp){

               printf("Do you want to write cholesky decomposition to file for future use, 'y' or 'n'\n\n");

               scanf("\n");
               if(getchar() == 'y'){

                  printf("What is the unique filename?\n\n");
                  scanf("%s", filnam); 
                  strcpy(filchol, FSPPATH);
                  strcat(filchol, filnam);

                  fchol = open_file(filchol, "w");

                  for(i=0; i<ncof; i++)
                     fwrite(aa[i], (i+1)*sizeof(double), 1, fchol);

                  close_file(fchol, filchol);


               }

            }


         }

         else{

            printf("****ERROR****, cholesky decomposition failed\n\n");

         }

       }

    }

    else if(im < 0){       

       if(!immp){
          for(i=0; i < dim; i++) free(*(aa+i));
       }
       else {
         if(munmap((void *)amp, memuse) == -1){
            printf("****ERROR****, unable to unmap memory mapped file.\n\n");
            exit(1);
         }
         close(fid);

       }
       free(aa);
       free(vec);
       free(bb);
       free(ind);

       return NULL;

    }



/* initialize RHS vector */

    for(i=0; i < ncof; i++) *(bb + i) = 0.0;

    atp = ap;

/* compute the RHS */


    for(i=0; i < gr->iy; i++){

        leng = gr->aleng + i * gr->nleng;

        for(j=0; j < nx; j++){

            lng = gr->sico + j * ntr1;

            np = 0;

            k = 0;

            atp = ap + i * gr->ix + j;

            for(m=0; m < ntr1; m++){

                cleng = leng + np;

                for(n=m; n < ntr1; n++){


                    vec[k++] = *cleng * lng->real;
 
                    if(m > 0) vec[k++] = *cleng * lng->imag; 
  
                    ++cleng;


                }

                ++lng;

                np += (ntr1 - m);


            }


            for(m=0; m < ncof; m++)

                bb[m] += *atp * vec[m]; 
 
        }

    }

/* for(i=0; i<ncof; i++)printf("%f\n", bb[i]); */


/* First back substitution to solve Gy=b, G lower triangular. */


    bb[0] /= *aa[0];

    for(i=1; i< ncof; i++){

        sum = 0.0;

        for(j=0; j < i; j++) sum += *(aa[i] + j) * bb[j]; 

        bb[i] -= sum;

        bb[i] /= *(aa[i] + i);

    }

/* Second back substitution to solve G'x=y */

    bb[ncof - 1] /= *(aa[ncof - 1] + ncof - 1);
    for(i=ncof-2; i >= 0; i--){

        sum = 0.0;

        for(j=i+1; j < ncof; j++) sum += *(aa[j] + i) * bb[j];

        bb[i] -= sum;
        bb[i] /= *(aa[i] + i);


    }
 
    return bb;

}






