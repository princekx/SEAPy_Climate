#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include "grid.h"
#include "mem_er.h"
#include "utf.h"
#include "pp.h"
#include "files_out.h"
#include "file_handle.h"
#include "netcdf_info.h"

#define  TOLSCL       1.0e-20

extern float *abuf;
extern PP *pph;

extern char *chrfld;

extern char *ihead[NUMLINES];
extern int nlines;

extern GRID *gr;
extern int form;
extern int eqsw;
extern int utfv;
extern int i_utf4x, i_utf4y;
extern int std_x, std_y;


int rf(FILE * ,FILE * , int , int );
void writef(FILE * , float * , int );
void read_chr_field(FILE * , int , int , char * , int );

/* function to combine weighted statistics 
   Note, reads multiple files with one field to each file        */

void wtm_combine(FILE *fdat)
{
    int i,j;
    int dim=0;
    int istat[3]={0, 0, 0};
    int ifile[5]={0, 0, 0, 0, 0};
    int iuwght=1, iswght=1;

    int nyear=0;

    float uwght=1.0, swght=1.0;

    char *infile[5]={NULL, NULL, NULL};
    char *outfile[5];

    FILE *fin[5]={NULL, NULL, NULL};

    float *avg=NULL, *var=NULL, *varfil=NULL, *wsum=NULL;
    float *avgt=NULL, *vart=NULL, *varfilt=NULL, *wt=NULL;

    float specwt=0.0, specsum=0.0;

    double a1, a2, a3, a4, a5, a6, diff;

    nlines = 0;

    for(i=0; i < NUMLINES; i++){
       ihead[i] = (char *)calloc(MAXCHR, sizeof(char));
       mem_er((ihead[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
    }

    for(i=0; i<5; i++){
      infile[i] = (char *)calloc(MAXCHR, sizeof(char));
      mem_er((infile[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
      outfile[i] = (char *)calloc(MAXCHR, sizeof(char));
      mem_er((outfile[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
    }

/* setup output file paths and names */

    strncpy(outfile[0], Add(USER,PATHO,comp_tavg), MAXCHR);
    strncpy(outfile[1], Add(USER,PATHO,comp_var), MAXCHR);
    strncpy(outfile[2], Add(USER,PATHO,comp_wsum), MAXCHR);
    strncpy(outfile[3], Add(USER,PATHO,comp_varfil), MAXCHR);
    strncpy(outfile[4], Add(USER,PATHO,comp_specwt), MAXCHR);


    printf("How many years are there to combine?\n\n");
    scanf("%d", &nyear);

/* array dimension */

    if(form < 2 || form == 4) dim = std_x * std_y;
    else if (form == 2){
      if(utfv == 3) dim=gr->ix * ((eqsw) ? gr->iy + 1 : gr->iy);
      else if(utfv == 4) dim = i_utf4x * i_utf4y;

    }

    else if(form == 3) dim = pph->lbrow * pph->lbnpt;

    printf("Select statistics to combine, 'a' for all and 's' for selected.\r\n"
           "It is the users responsibilty that all files are available.    \n\n");
    scanf("\n");
    if(getchar() == 'a'){
       istat[0] = istat[1] = istat[2] = 1;
    }
    else{
       printf("Do you want to combine means, '1' for yes, '0' for no.\n\n");
       scanf("%d", istat);
       printf("Do you want to combine variances, '1' for yes, '0' for no.\n\n");
       scanf("%d", istat + 1);
       printf("Do you want to combine filtered variances, '1' for yes, '0' for no.\n\n");
       scanf("%d", istat + 2);
    }

    if(istat[0] || istat[1]){
       avg = (float *)calloc(dim, sizeof(float));
       mem_er((avg == NULL) ? 0 : 1, dim * sizeof(float));
       avgt = (float *)calloc(dim, sizeof(double));
       mem_er((avgt == NULL) ? 0 : 1, dim * sizeof(float));
       wsum = (float *)calloc(dim, sizeof(float));
       mem_er((wsum == NULL) ? 0 : 1, dim * sizeof(float));
       wt = (float *)calloc(dim, sizeof(double));
       mem_er((wt == NULL) ? 0 : 1, dim * sizeof(float));
       ifile[2] = 1;
       printf("Do you want to use a fixed user chosen weight for means and \r\n"
              "variance compositing, 'y' or 'n'                            \n\n");
       scanf("\n");
       if(getchar() == 'y'){
          printf("What is the chosen weight?\n\n");
          scanf("%f", &uwght);
          for(i=0; i < dim; i++){
             *(wsum+i) = *(wt + i) = 1.0;
          }
          iuwght = 0;
          ifile[2] = 0;
       }

       ifile[0] = 1;

    }
    if(istat[1]){
       var = (float *)calloc(dim, sizeof(float));
       mem_er((var == NULL) ? 0 : 1, dim * sizeof(float));
       vart = (float *)calloc(dim, sizeof(double));
       mem_er((vart == NULL) ? 0 : 1, dim * sizeof(float));
       ifile[1] = 1;
    }
    if(istat[2]){
       varfil = (float *)calloc(dim, sizeof(float));
       mem_er((varfil == NULL) ? 0 : 1, dim * sizeof(float));
       varfilt = (float *)calloc(dim, sizeof(float));
       mem_er((varfilt == NULL) ? 0 : 1, dim * sizeof(float));
       ifile[4] = 1;
       printf("Do you want to use a fixed user chosen weight for filtered variance\r\n"
              "compositing, 'y' or 'n'                                            \n\n");
       scanf("\n");
       if(getchar() == 'y'){
          printf("What is the chosen weight?\n\n");
          scanf("%f", &swght);
          specwt = 1.0;
          iswght = 0;
          ifile[4] = 0;

       }
       ifile[3] = 1;

    }
  
    for(i=0; i < nyear; i++){
        if(istat[0] || istat[1]){
           printf("Input next filename for time average.\n\n");
           scanf("%s", infile[0]);
           if(iuwght){
              printf("Input next filename for sum of weights.\n\n");
              scanf("%s", infile[2]);
           }
        }
        if(istat[1]){
           printf("Input next filename for STD.\n\n");
           scanf("%s", infile[1]);
        }
        if(istat[2]){
           printf("Input next filename for the filtered varience.\n\n");
           scanf("%s", infile[3]);
           if(iswght){
              printf("Input next filename for the spectral filter weight.\n\n");
              scanf("%s", infile[4]);
           }
        }

        if(form != 4){
           for(j=0; j < 4; j++){
              if(ifile[j]) fin[j] = open_file(infile[j], "r");
           }
        }
        else {
           for(j=0; j < 4; j++){
              if(ifile[j]){
                fin[j] = (FILE *)nc_clone((NETCDF_INFO *)fdat, infile[j], NC_OPEN_MODE);
                ((NETCDF_INFO *)fin[j])->iframe = 0;
                ((NETCDF_INFO *)fin[j])->levval = 0;
              }
           }
        }

        if(ifile[4]) {
           fin[4] = open_file(infile[4], "r");
           fscanf(fin[4], "%f", &specwt);
        }

        if(i == 0){

           nlines = 0;

           if(avgt){
              abuf = avgt;
              rf(fin[0], NULL, dim, 1);
           }
           if(vart){
              abuf = vart;
              rf(fin[1], NULL, dim, 1);
           }
           if(iuwght){
             abuf = wt;
             rf(fin[2], NULL, dim, 1);
           }
           if(varfilt){
             abuf = varfilt;
             rf(fin[3], NULL, dim, 1);
           }
           specsum = specwt;

        }

        else {

           if(avg){
              abuf = avg;
              rf(fin[0], NULL, dim, 0);
           }
           if(var){
              abuf = var;
              rf(fin[1], NULL, dim, 0);
           }
           if(iuwght){
              abuf = wsum;
              rf(fin[2], NULL, dim, 0);
           }
           if(varfil){
              abuf = varfil;
              rf(fin[3], NULL, dim, 0);
           }

           a6 = specsum + specwt;

           if(vart){
              for(j=0; j<dim; j++){
                 a1 = *(wt + j) * *(avgt + j) + *(wsum + j) * *(avg + j);
                 a2 = *(wt + j) * *(vart + j) * *(vart + j) + *(wsum + j) * *(var + j) * *(var + j);
                 diff = *(avgt + j) - *(avg + j);
                 a3 = *(wt + j) * *(wsum + j) * diff * diff;
                 a4 = *(wt + j) + *(wsum + j);

                 *(avgt + j) = a1 / a4;
                 *(vart + j) = sqrt((a2 / a4) + (a3 / (a4 * a4)));
                 *(wt + j) = a4;

              }
           }

           else if(avgt){
              for(j=0; j<dim; j++){
                 a1 = *(wt + j) * *(avgt + j) + *(wsum + j) * *(avg + j);
                 a4 = *(wt + j) + *(wsum + j);
                 *(avgt + j) = a1 / a4;
                 *(wt + j) = a4;

              }

           }

           if(varfilt){

              for(j=0; j<dim; j++){

                 a5 = (specsum * *(varfilt + j) + specwt * *(varfil + j)) / a6;

                 *(varfilt + j) = a5;
           
              }

           }

           specsum += specwt;

        }

        if(form != 4) {
           for(j=0; j < 4; j++)
              if(fin[j]) close_file(fin[j], infile[j]);
        }
        else {
           for(j=0; j < 4; j++)
               if(fin[j]) netcdf_close((NETCDF_INFO *)fin[j]);
        }

        if(fin[4]) close_file(fin[4], infile[4]);

    }

    if(form != 4){

       for(i=0; i<4; i++) {
           if(ifile[i]){
              fin[i] = open_file(outfile[i], "w");
              for(j=0; j < nlines; j++)
                  fprintf(fin[i], "%s", ihead[j]);
           }
       }

    }
    else{
      for(i=0; i < 4; i++){
         if(ifile[i]){
            fin[i] = (FILE *)nc_define((NETCDF_INFO *)fdat, outfile[i]);
            ((NETCDF_INFO *)fin[i])->nframe = 0;
            ((NETCDF_INFO *)fin[i])->iframe = 0;
         }
      }
    }

    if(ifile[4])fin[4] = open_file(outfile[4], "w");


    if(fin[0]) writef(fin[0], avgt, dim);
    if(fin[1])writef(fin[1], vart, dim);
    if(fin[2])writef(fin[2], wt, dim);
    if(fin[3])writef(fin[3], varfilt, dim);
    if(fin[4])fprintf(fin[4], "%e\n", specsum);

    if(form != 4){
       for(i=0; i < 4; i++) {
           if(fin[i]) close_file(fin[i], outfile[i]);
       }
    }
    else{
       for(i=0; i < 4; i++)
          if(fin[i]) netcdf_close((NETCDF_INFO *)fin[i]);
    }
    if(fin[4])close_file(fin[4], outfile[4]);

    free(avg);
    free(var);
    free(wsum);
    free(avgt);
    free(vart);
    free(wt);
    free(varfil);
    free(varfilt);

    free(chrfld);

    for(i=0; i<NUMLINES; i++) free(ihead[i])
    for(i=0; i < 5; i++){free(infile[i]); free(outfile[i]);}

    return;

}
