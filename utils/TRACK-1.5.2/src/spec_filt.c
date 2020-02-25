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

#define  LARGE_FREQ   1.0e+6

/* function to compute time filtered data in the spectral domain. */

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
void dfft(complex * , int , int , int );
int powi(int , int );
int missing(float, float, int);


void spec_filt(FILE *fdat, int fr1, int fri, int frl, int nn)

{

    int i, j, k, l;
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
    int nfilt;
    int itrig=1;
    int tsamp=0;
    int nfreq=0, nfreq2=0;
    int nexp=0;
    int **index=NULL;
    int app='n';
    int ident=0, idn=0;
    int icmp=0;


    off_t place1, place2, chrnum=0;
    off_t pl1=0, plw=0, plt=0;

    float **bltim=NULL, *btim=NULL;
    float *bll=NULL;
    float *buft=NULL;
    float *ata=NULL, *lta=NULL;
    float *bands=NULL;
    float ***blband=NULL;
    float *freq=NULL;
    float lowper;

    float **aband=NULL, **atemp=NULL;

    float mval=MVAL, nval=0.0;


    FILE **ffilt=NULL;
    FILE *timav=NULL;

    complex *timser=NULL, *specser=NULL;

    char filename[MAXCHR], band_id[NCHRB];

    nlines = 0;

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


    printf("How many spectral bands are required?\n\n");
    scanf("%d", &nfilt);

    if(nfilt <= 0){

       printf("****ERROR****, number of bands not >= 1, exiting.\n\n");
       exit(1);

    }

    ffilt = (FILE **)calloc(nfilt, sizeof(FILE *));
    mem_er((ffilt == NULL) ? 0 : 1, nfilt * sizeof(FILE *));

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

    bands = (float *)calloc(nfilt+1, sizeof(float));
    mem_er((bands == NULL) ? 0 : 1, (nfilt+1) * sizeof(float));

    printf("What is the number of data samples per time unit (e.g. day) for period?\n\n");
    scanf("%d", &tsamp);

    if(tsamp < 1){

       printf("****ERROR****, samples per unit time must be positive.\n\n");
       exit(1);

    }

    printf("Input the period boundaries for filtering, e.g., 0 10 100,\r\n"
           " for 2 bands between 0-10 time units and 10-100 time units\r\n"
           " ordered as high pass to low pass in frequency.           \r\n"
           "Ouput is in the order high pass to low pass.              \n\n");

    lowper = 2.0 / (float) tsamp;

    printf("****INFORMATION****, the lowest period sampled is %f\n", lowper);

    i = 0;

    while(i < (nfilt + 1)) {
         printf("Boundary %d =", i+1);
         scanf("%f", (bands+i));
         if(*(bands + i) < lowper) *(bands + i) = 0.5;
         else if(*(bands + i) > LARGE_FREQ) *(bands + i) = 0.0;
         else *(bands + i) = 1.0 / ((float) tsamp * *(bands+i));
         printf("Frequency = %f", *(bands + i));
         printf("\n");
         ++i;
    }

    printf("\n");

/* allocate memory */

    bltim = (float **)calloc(npoint, sizeof(float *));
    mem_er((bltim == NULL) ? 0 : 1, npoint * sizeof(float *));


    atemp = (float **)calloc(nfilt, sizeof(float *));
    mem_er((atemp == NULL) ? 0 : 1, nfilt * sizeof(float *));

    aband = (float **)calloc(nfilt, sizeof(float *));
    mem_er((aband == NULL) ? 0 : 1, nfilt * sizeof(float *));

    for(i=0; i < nfilt; i++){

       *(aband + i) = (float *)calloc(dim, sizeof(float));
       mem_er((*(aband + i) == NULL) ? 0 : 1, dim * sizeof(float));


    }


    for(i=0; i < npoint; i++){

       *(bltim + i) = (float *)calloc(nn, sizeof(float));
       mem_er((*(bltim + i) == NULL) ? 0 : 1, nn * sizeof(float));

    }

    blband = (float ***)calloc(nfilt, sizeof(float *));
    mem_er((blband == NULL) ? 0 : 1, nfilt * sizeof(float *));

    for(i=0; i < nfilt; i++){

       *(blband + i) = (float **)calloc(npoint, sizeof(float *));
       mem_er((*(blband + i) == NULL) ? 0 : 1, npoint * sizeof(float *));

       for(j=0; j < npoint; j++){

           *(*(blband + i) + j) = (float *)calloc(nn, sizeof(float));
           mem_er((*(*(blband + i) + j) == NULL) ? 0 : 1, nn * sizeof(float));

       }

    }


    printf("How many frequencies are required for the spectral filtering, usualy some power of 2?\n\n");

    while(nfreq <= 0){
         scanf("%d", &nfreq);
         if(nfreq <= 0)
           printf("****ERROR****, number of frequencies must be positive, try again.\n\n");

    }

    if(nfreq < nn) {

       printf("****WARNING****, number of frequencies less than data length,\r\n"
              "                 resetting to data length.                   \n\n");

       nfreq = nn;

    }


    nexp = (int)(log((float)nfreq)/log(2.0)); 

    if(powi(2, nexp) != nfreq){
       ++nexp;
       nfreq = powi(2, nexp);

       printf("****WARNING****, number of frequencies must be a power of 2,\r\n" 
              "                 resetting to %d\n\n", nfreq);

    }

    nfreq2 = nfreq / 2;

    if(nn < nfreq){

       printf("Do you want to pre and post append data, helps with end effects, 'y' or 'n'.\r\n");
       scanf("\n");
       app = getchar();
       ident = (nfreq - nn) / 2;
       idn = ident;
       if(idn > nn) idn = nn - 1;
    }


    timser = (complex *)calloc(nfreq, sizeof(complex));
    mem_er((timser == NULL) ? 0 : 1, nfreq * sizeof(complex));

    specser = (complex *)calloc(nfreq, sizeof(complex));
    mem_er((specser == NULL) ? 0 : 1, nfreq * sizeof(complex));



    freq = (float *)calloc(nfreq, sizeof(float));
    mem_er((freq == NULL) ? 0 : 1, nfreq * sizeof(float));

    index = (int **)calloc(nfilt, sizeof(int *));
    mem_er((index == NULL) ? 0 : 1, nfilt * sizeof(int *));

    for(i=0; i<nfilt; i++){

       *(index + i) = (int *)calloc(nfreq, sizeof(int));
       mem_er((*(index + i) == NULL) ? 0 : 1, nfreq * sizeof(int));

    }

    for(i=0; i < nfreq; i++)

        *(freq + i) = -0.5 + (float)i / (float)nfreq;

/* apply band masks, ordering is for decimation in time FFT */

    for(i=0; i<nfilt; i++){


       for(j=0; j<nfreq2; j++){

           if(fabs(*(freq + j)) <= *(bands + i) && fabs(*(freq + j)) >= *(bands + i + 1)) *(index[i] + j + nfreq2) = 1;

           else *(index[i] + j + nfreq2) = 0;

          if(fabs(*(freq + j + nfreq2)) <= *(bands + i) && fabs(*(freq + j + nfreq2)) >= *(bands + i + 1)) *(index[i] + j) = 1;
           else *(index[i] + j) = 0;

       }


    }


    printf("***INFORMATON***, Computing Time Filtering of Data....\r\n\n"
           "                                                      \r\n\n"
           " Please Wait, this may take some time ..............  \n\n");

    chksum = 0;
    if(form != 4) pl1 = ftello(fdat);

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

/* position file pointer */

      chksum = ((i + 1) * npoint < dim) ? npoint : dim - i * npoint; 

      if(form != 4) fseeko(fdat, pl1, FSTART);

      bll = abuf + i * npoint;
      if(ata) lta = ata + i * npoint;

      for(j=0; j<nfilt; j++) atemp[j] = *(aband + j) + i * npoint;

      itim = 0;
      nf = 0;

      while(nf <= frl){

           if(form != 4) {

              if(!nf) {

                 place1 = ftello(fdat); 

                 rf(fdat, NULL, dim, hf);         

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
                    rf(fdat, NULL, dim, hf);

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

                  if(missing(*btim, mval, icmp))
                     *btim = (ata) ? *(lta + j) : nval;

              }               

           }

           ++itim;

           if(nf > frl) break;

      }

      if(itim != nn){

         printf("****ERROR****, number of read time steps does not equal\r\n"
                "               previously read number.                 \n\n");
         exit(1);

      }

      hf = 0;


      for(j=0; j<chksum; j++){
 
          btim = *(bltim + j);
          for(k=0; k<nfreq; k++)(timser+k)->real = (timser+k)->imag = 0.0;
          for(k=0; k < nn; k++) comp((*(btim + k)), 0.0, (timser+k+ident));

          if(app == 'y'){
             for(k=0; k < idn; k++){
                 comp((*(btim + k + 1)), 0.0, (timser+ident - k - 1));
                 comp((*(btim + nn - k - 2)), 0.0, (timser + nn + ident + k));
             }

          }


          dfft(timser, nexp, 1, itrig);



          for(k=0; k < nfilt; k++){

              btim = *(*(blband + k) + j);

              for(l=0; l<nfreq; l++) {

                 cmx((float)*(index[k] + l), *(timser + l), (specser + l));

              }

              dfft(specser, nexp, -1, itrig);

              itrig = 0;

              for(l=0; l<nn; l++) *(btim + l) = (specser + l + ident)->real; 


          }


      }

      buft = abuf;

      if(form != 4){
         for(j=0; j < nfilt; j++) fseeko(ffilt[j], plw, FSTART);
      }

      for(j=0; j<nn; j++){

         for(k=0; k<nfilt; k++){

            if(form != 4) plt = ftello(ffilt[k]);
            else ((NETCDF_INFO *)ffilt[k])->iframe = j;

            abuf = aband[k];
            rf(ffilt[k], NULL, dim, 1);

            for(l=0; l<chksum; l++) *(atemp[k] + l) = *(*(blband[k] + l) + j);

            if(form != 4) fseeko(ffilt[k], plt, FSTART);
            else {
               ((NETCDF_INFO *)ffilt[k])->iframe = j;
               ((NETCDF_INFO *)ffilt[k])->nframe = j;
            }
            for(l=0; l<nlines; l++) fprintf(ffilt[k], "%s", ihead[l]);
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

    dfft(NULL, 0, 1, -1);

    for(i=0; i < npoint; i++) free(*(bltim + i));
    free(bltim);
    free(atemp);
    free(ffilt);
    free(timser);
    free(freq);
    free(specser);

    for(i=0; i<nfilt; i++){

        for(j=0; j < npoint; j++) free(*(*(blband + i) + j));
        free(*(blband + i));
        free(*(index + i));
        free(*(aband + i));
    }

    free(blband);
    free(index);
    free(aband);
    free(ata);
    free(bands);

    for(i=0; i<NUMLINES; i++) free(ihead[i]);

    if(chrfld) free(chrfld);

    return;

}
