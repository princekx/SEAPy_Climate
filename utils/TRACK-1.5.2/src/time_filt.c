#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include <sys/types.h>
#include "mem_er.h"
#include "grid.h"
#include "files_out.h"
#include "file_handle.h"
#include "pp.h"
#include "netcdf_info.h"

#define  NCHRB  30


/* function to compute filtered data in the time domain. */

extern GRID *gr;
extern int frnum;
extern int form;
extern int eqsw;
extern int utfv;
extern int i_utf4x, i_utf4y;
extern int std_x, std_y;

extern float *abuf;
extern char *chrfld;

extern PP *pph;

extern char *fext;

extern int iext;

extern char *ihead[NUMLINES];
extern int nlines;

int rf(FILE * ,FILE * , int , int );
void writef(FILE * , float * , int );
int missing(float, float, int);

void time_filt(FILE *fdat, double **weights, int nfilt, int nor, int fr1, int fri, int frl, int nn)

{

    int i, j, k;
    int m, n;
    int neg, pos;
    int ntot = 0;
    int dim=0;
    int npoint=0;
    int nblock=0;
    int chksum=0;
    int itim=0;
    int nf=0;
    int wf='a';
    int hf=1;
    int imiss='n';
    int icmp=0;

    off_t place1, place2, chrnum=0;
    off_t pl1=0, plw=0, plt=0;

    float **bltim=NULL, *btim=NULL;
    float *bll=NULL;
    float *buft=NULL;
    float *ata=NULL, *lta=NULL;

    float **aband=NULL, **atemp=NULL;

    float mval=MVAL, nval=0.0;

    double sum=0.0;
    double dneg, dpos;

    FILE **ffilt=NULL;
    FILE *timav=NULL;

    char filename[MAXCHR], band_id[NCHRB];

    nlines = 0;

    if(nfilt <= 0){

       printf("****ERROR****, number of bands not >= 1, exiting.\n\n");
       exit(1);

    }

    ffilt = (FILE **)calloc(nfilt, sizeof(FILE *));
    mem_er((ffilt == NULL) ? 0 : 1, nfilt * sizeof(FILE *));

    printf("*******COMPUTING FILTERED FIELD**********\n\n");


    printf("****WARNING****, computing the filtered fields may require\r\n"
           "                 large amounts of memory.                 \n\n");


    printf("Do you want all frames in data file to contribute to the  \r\n"
           "filtering, or just those frames selected or all frames \r\n"
           "at the specified sampling rate from the starting frame,   \r\n"
           " 'a', 's', or 'i'. Make sure to use the same option as for\r\n"
           "the time average calculation.                             \n\n"); 

    scanf("\n");
    wf = getchar();

    if(wf == 'a' ){fr1 = 1; fri = 1; frl = frnum;}

    else if(wf == 's'){
       if(fri < 1) fri = 1;
       if(frl > frnum) frl = frnum;
    }

    else if(wf == 'i'){
       if(fri < 1) fri = 1;
       frl = frnum;
    }

    else {

       printf("***WARNING***, incorrect decriptor for time average calculation\r\n"
             "               defaulting to 'a' (all).                        \n\n");

       wf = 'a';
       fr1 = 1; fri = 1; frl = frnum;

    }

    printf("***INFORMATION***, frame counters for computing time average are:\r\n"
           "                   Start = %d, Interval = %d, Last = %d \n\n", fr1, fri, frl);


/* assign space for header strings */

    nlines = 0;

    for(i=0; i < NUMLINES; i++){
       ihead[i] = (char *)calloc(MAXCHR, sizeof(char));
       mem_er((ihead[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
    }

/* grid dimension */

    if(form < 2 || form == 4) dim = std_x * std_y;
    else if (form == 2){
       if(utfv == 3) dim=gr->ix * ((eqsw) ? gr->iy + 1 : gr->iy);
       else if(utfv == 4) dim = i_utf4x * i_utf4y;

    }

    else if(form == 3) dim = pph->lbrow * pph->lbnpt; 

/* missing value assignment */

    if(form == 4){

       if(((NETCDF_INFO *)fdat)->imiss){
          printf("****INFORMATION****, missing data value is %e\n\n", ((NETCDF_INFO *)fdat)->missing);
          mval = ((NETCDF_INFO *)fdat)->missing;
          imiss = 'y';

       }

    }
    else {

       printf("Are there missing data values, delineated by some extreme value. 'y' or 'n'\n\n");

       scanf("\n");
       imiss=getchar();

    }

    if(imiss == 'y'){

      printf("What is the missing data value?\n\n");
      scanf("%f", &mval);

      printf("How do you want to compare with missing value?\r\n"
             "'0' -- equality.                              \r\n"
             "'1' -- less than.                             \r\n"
             "'2' -- greater than.                          \n\n");
      scanf("%d", &icmp);

      printf("Do you want to replace the missing value by the mean value, 'y' or 'n'\n\n");
      scanf("\n");
      if(getchar() == 'y'){

         ata = (float *)calloc(dim, sizeof(float));
         mem_er((ata == NULL) ? 0 : 1, dim * sizeof(float));

         buft = abuf;
         abuf = ata;

         if(form != 4) timav = open_file(TAVGE, "r"); 
         else {
            timav = (FILE *)nc_clone((NETCDF_INFO *)fdat, TAVGE, NC_OPEN_MODE);
            ((NETCDF_INFO *)timav)->iframe = 0;
            ((NETCDF_INFO *)timav)->nframe = 0;
         }

         rf(timav, NULL, dim, 0);         

         if(form != 4) close_file(timav, TAVGE); 
         else netcdf_close((NETCDF_INFO *)timav);    

         abuf = buft;

      }

      else {

         printf("What value do you want to replace missing value by?\n\n");
         scanf("%f", &nval);

      }

    }

    printf("****INFORMATION****, the dimension of the data is = %d\n\n", dim);

    printf("How many filter evaluations per block, large number requires large memory,\r\n"
           "small number requires more I/O.                                           \n\n");

    scanf("%d", &npoint);

    if(npoint > dim) npoint = dim;

    nblock = (dim % npoint) ? (dim / npoint) + 1 : dim / npoint;

/* allocate memory */

    bltim = (float **)calloc(npoint, sizeof(float *));
    mem_er((bltim == NULL) ? 0 : 1, npoint * sizeof(float *));

    aband = (float **)calloc(nfilt, sizeof(float *));
    mem_er((aband == NULL) ? 0 : 1, nfilt * sizeof(float *));

    atemp = (float **)calloc(nfilt, sizeof(float *));
    mem_er((atemp == NULL) ? 0 : 1, nfilt * sizeof(float *));

    for(i=0; i < nfilt; i++){

       *(aband + i) = (float *)calloc(dim, sizeof(float));
       mem_er((*(aband + i) == NULL) ? 0 : 1, dim * sizeof(float));


    }

    for(i=0; i < npoint; i++){

       *(bltim + i) = (float *)calloc(nn, sizeof(float));
       mem_er((*(bltim + i) == NULL) ? 0 : 1, nn * sizeof(float));

    }


    printf("***INFORMATON***, Computing Time Filtering of Data....\r\n\n"
           "                                                      \r\n\n"
           " Please Wait, this may take some time ..............  \n\n");

    chksum = 0;
    pl1 = ftello(fdat);

/* open files for writing filtered fields */

    for(i=0; i < nfilt; i++){

       sprintf(band_id, "_band%03d", i);
       strncpy(filename, TIME_FILT, MAXCHR);
       if(iext) strcpy(strstr(filename, EXTENSION), fext);
       strcat(filename, band_id);

       if(form != 4) ffilt[i] = open_file(filename, "w+");
       else ffilt[i] = (FILE *)nc_define((NETCDF_INFO *)fdat, filename);

    }

    if(form != 4) plw = ftello(ffilt[0]);

    for(i=0; i < nblock; i++){

      chksum = ((i + 1) * npoint < dim) ? npoint : dim - i * npoint; 

/* position file pointer */

      if(form != 4) fseeko(fdat, pl1, FSTART);

      bll = abuf + i * npoint;
      lta = ata + i * npoint;

      for(j=0; j<nfilt; j++) atemp[j] = *(aband + j) + i * npoint;

      itim = 0;
      nf = 0;

      while(nf <= frl){

           if(form != 4){

              if(!nf) {

                 place1 = ftello(fdat); 

                 if(rf(fdat, NULL, dim, hf)) break;         
 
                 if(fr1 == 1) {

                    if(hf){

                       for(j=0; j<nfilt; j++){
 
                         for(k=0; k<nlines; k++) fprintf(ffilt[j], "%s", ihead[k]);
                         writef(ffilt[j], aband[j], dim);

                       }

                    }

                    ++nf;

                 }
 
                 place2 = ftello(fdat);
                 chrnum = place2 - place1;

                 if(fr1 > 1){

                    fseeko(fdat, (fr1-2)*chrnum, ORIGIN);
                    if(rf(fdat, NULL, dim, hf)) break;;

                    if(hf){

                       for(j=0; j<nfilt; j++){

                         for(k=0; k<nlines; k++) fprintf(ffilt[j], "%s", ihead[k]);
                         writef(ffilt[j], aband[j], dim);

                       }

                    }

                    nf = fr1;           

                 }

                 nf += fri;


              }

              else {

                 fseeko(fdat, (fri-1)*chrnum, ORIGIN);

                 if(rf(fdat, NULL, dim, hf)) break;

                 if(hf){

                   for(j=0; j<nfilt; j++){

                      for(k=0; k<nlines; k++) fprintf(ffilt[j], "%s", ihead[k]);
                      writef(ffilt[j], aband[j], dim);
                      fflush(ffilt[j]);

                    }

                 }
           
                 nf += fri;

              }

           }

           else {


              if(!nf) nf = fr1;

              ((NETCDF_INFO *)fdat)->iframe = nf - 1;

              if(rf(fdat, NULL, dim, hf)) break;

              if(hf){

                for(j=0; j<nfilt; j++){

                   ((NETCDF_INFO *)ffilt[j])->iframe = nf - 1;
                   ((NETCDF_INFO *)ffilt[j])->nframe = itim;

                   writef(ffilt[j], aband[j], dim);

                }

              }
           
              nf += fri;

           }

           if(imiss == 'n'){

              for(j=0; j < chksum; j++){

                  btim = *(bltim + j) + itim;

                  *btim = *(bll + j);

              }

           }

           else{

              for(j=0; j < chksum; j++){

                  btim = *(bltim + j) + itim;

                  *btim = *(bll + j);

                  if(missing(*btim, mval, icmp)) {
                     *btim = (ata) ? *(lta + j) : nval;

                  }

              }               

           }

           ++itim;

           if(nf > frl) break;

      }


      hf = 0;

      if(form != 4) {
        for(j=0; j < nfilt; j++) fseeko(ffilt[j], plw, FSTART);
      }

      buft = abuf;

      for(m=0; m < nn; m++){

         for(k=0; k < nfilt; k++) {
            if(form != 4) plt = ftello(ffilt[k]);
            else ((NETCDF_INFO *)ffilt[k])->iframe = m;

            abuf = aband[k];

            if(rf(ffilt[k], NULL, dim, 1)) break;
         }


         for(k=0; k<chksum; k++){
 
            btim = *(bltim + k);

            for(j=0; j<nfilt; j++){

               sum = **(weights + j) * *(btim+m);

               for(n=1; n <= nor; n++){

                  neg = m - n;
                  if(neg < 0) dneg = *(btim - neg);
                  else dneg = *(btim + neg);
                  pos = m + n;
                  if(pos >= nn) dpos = *(btim + 2 * nn - pos - 2);
                  else dpos = *(btim + pos);

                  sum += *(*(weights + j) + n) * (dneg + dpos); 

               }

               *(atemp[j] + k) = sum;
 

            }


         }

         for(k=0; k < nfilt; k++){

            if(form != 4) fseeko(ffilt[k], plt, FSTART);
            else {
               ((NETCDF_INFO *)ffilt[k])->iframe = m;
               ((NETCDF_INFO *)ffilt[k])->nframe = m;  
            }

            for(j=0; j<nlines; j++) fprintf(ffilt[k], "%s", ihead[j]);          

            writef(ffilt[k], aband[k], dim);

            if(form != 4) fflush(ffilt[k]);
         }

      }


      abuf = buft;

      ntot += chksum;

      printf("Completed %d points so far.\n", ntot);

    }

    for(i=0; i < nfilt; i++){

       sprintf(band_id, "_band%03d", i);
       strncpy(filename, TIME_FILT, MAXCHR);
       if(iext) strcpy(strstr(filename, EXTENSION), fext);
       strcat(filename, band_id);
       if(form != 4) close_file(ffilt[i], filename);
       else netcdf_close((NETCDF_INFO *)ffilt[i]);

    }


    for(i=0; i < npoint; i++) free(*(bltim + i));
    free(bltim);
    for(i=0; i < nfilt; i++) free(*(aband + i));
    free(aband);
    free(atemp);
    free(ffilt);
    free(ata);

    for(i=0; i<NUMLINES; i++) free(ihead[i]);

    if(chrfld) free(chrfld);

    return;

}
