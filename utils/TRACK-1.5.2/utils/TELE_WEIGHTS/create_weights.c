#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"

#define DNUM  500
#define TOLWT   1.0e-6

/* combine sets of tr_trs files */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int * );

int meantrd_tanh(FILE * , struct tot_tr * , int , int , int * , int * , float * , int , int , int , float , float , float , float );

int aniso;

float tele_sum=0.0;

int noheader=0;

extern int nfld, nff;
extern int *nfwpos;

int main(int argc, char *argv[])

{


    int i, j, trc=0;
    int nf=0;
    int ifabs=0;
    int iwrt=0;
    int dnum=DNUM;
    int ireg=0;
    int iraw=1;
    int itanty=0;
    int ittmp=0;

    int istr[DNUM], iend[DNUM];
    float ival[DNUM];
    float *wval=NULL;

    float vali, hbw, nval;
    float ee=0.0;

    long int pl=0, plw=0, plww=0;

    int gpr, ipr;
    float alat, alng;
    float t_scale=0.0;
    float lat1, lat2, lng1, lng2;

    int trnum=0, trnr=0;
    int tr_count;
    int iwt=1, indx=0;

    char filnamin[MAXCHR];
    char filnamout[]="weights_trs";

    FILE *fin=NULL;
    FILE *fout=NULL;
    FILE *fdat=NULL;

    struct tot_tr *alltr=NULL;


    if(argc < 3){
       printf("USAGE: [dat file] [tanh scaling] [lat1] [lat2] [lng1] [lng2]\r\n"
              "       lat1, lat2, lng1, lng2 are optional.                 \n\n");
       exit(1);

    }

    sscanf(argv[2], "%f", &t_scale);

    if(argc > 3){
       sscanf(argv[3], "%f", &lat1);
       sscanf(argv[4], "%f", &lat2);
       sscanf(argv[5], "%f", &lng1);
       sscanf(argv[6], "%f", &lng2);

       if(lat2 < lat1 || lng2 < lng1){
          printf("****ERROR****, lat2 not greater than lat1, or lng2 not greater than lng1\n\n");
          exit(1);
       }

       printf("Do you want to apply weights as a function of time only, '1',\r\n"
              "or as a function of time and space (region) '2'.             \n\n");
       scanf("%d", &ireg);
    }

    fdat = fopen(argv[1], "r");
    if(!fdat){
      printf("****ERROR****, can't open file %s\n", argv[1]);
      exit(1);
    }

    fscanf(fdat, "%d", &dnum);
    for(i=0; i < dnum; i++)fscanf(fdat, "%d", istr+i);

    for(i=0; i < dnum; i++)fscanf(fdat, "%d", iend+i);


    fclose(fdat);

    printf("****INFORMATION****, weights can be transformed for computing \r\n"
           "                     teleconnection climatologies or left as  \r\n"
           "                     raw values for estimating statistics for \r\n"
           "                     a particular index value.                \n\n");

    printf("Do you want index values transformed or left as raw values, \r\n"
           "input '0' for raw or '1' for transformed.                   \n\n");

    scanf("%d", &iraw);
    if(iraw < 0 || iraw > 1) {
       printf("****ERROR****, incorrect value specified for type of weights.\n\n");
       exit(1);
    }

    fout = fopen(filnamout, "w");
    if(!fout){
      printf("****ERROR****, can't open file %s\n", filnamout);
      exit(1);
    }

    if(iraw){
       printf("Do you want one sided tanh weights (for +/- climatologies) \r\n"
              "or centered gaussian weights (for climatologies centered   \r\n"
              "around a particular value), input '0' or '1'.              \n\n");
       scanf("%d", &itanty);

    }

    if(itanty){
       printf("What is the value for the centered statistics?\n\n");
       scanf("%f", &vali);
       printf("What is the value of the half bandwidth?\n\n");
       scanf("%f", &hbw);

    }

    printf("How many files do you want to combine?\n");
    scanf("%d", &nf);

    printf("%d\n", nf);

    printf("Do you want to use absolute values for the index values, \r\n"
           "used for significance testing, 'y' or 'n',               \n\n");

    scanf("\n");
    if(getchar() == 'y') {

       if(t_scale < 0){
          printf("****WARNING****, scaling factor is negative for use with \r\n"
                 "                 absolute values of indicies.            \n\n");

       }

       ifabs = 1;

    }

    if(iraw){

       printf("There are %d time intervals for weights, do you want to use\r\n"
              "them all or just one interval, 'a' or 's'.                 \n\n", dnum);

       scanf("\n");
       if(getchar() == 's'){

          iwt = 0;

          printf("Which index do you want?\n\n");
          scanf("%d", &indx);

          if(indx < 1 || indx > dnum){
             printf("****ERROR****, chosen index is not in range.\n\n");
             exit(1);
          }

          --indx;
       }

       else iwt = 1;

       printf("Do you want the individual time interval weights written to the header, 'y' or 'n'\n\n");
       scanf("\n");
       if(getchar() == 'y') {
          iwrt = 1;
          wval = (float *)calloc(dnum*nf, sizeof(float));
          mem_er((wval == NULL) ? 0 : 1, dnum*nf*sizeof(float));
       }

    }



    while(trc < nf){

       printf("What is the next tr_trs file?\n");
       if(scanf("%s", filnamin) < 0){
          printf("****WARNING****, insufficient files for chosen number.\n\n");
          break;
       }
       printf("%d %s\n", trc, filnamin);

       fin=fopen(filnamin, "r");
       if(!fin){
          printf("****ERROR****, can't open file %s\n", filnamin);
          exit(1);
       }

       printf("What are the %d teleconnection indicies?\n", dnum);

       if(iraw) {

          if(!itanty){

             for(i=0; i<dnum; i++) {

                 scanf("%f", ival+i);

                 if(ifabs) *(ival + i) = fabs(*(ival + i));
  
                 *(ival + i) *= t_scale;

                 if(iwt){
                    if(*(ival + i) > TOLWT) tele_sum += (*(ival + i) = tanh(*(ival+i)));
                    else *(ival + i) = 0.0;
              
                 }
                 else {
                    if(*(ival + i) > TOLWT && i == indx) tele_sum += (*(ival + i) = tanh(*(ival+i)));
                    else *(ival + i) = 0.0;
                 }

                 if(wval) *(wval + trc * dnum + i) = *(ival + i);

             }

          }

          else {
             for(i=0; i<dnum; i++) {
                 scanf("%f", ival+i);
                 nval = (*(ival + i) - vali) / hbw;
                 if(fabs(nval) < 1.0) {
                    ee = 1.0 - nval * nval;
                    tele_sum += (*(ival + i) = ee * ee);
                 }
                 else *(ival + i) = 0.0;
                 
                 if(wval) *(wval + trc * dnum + i) = *(ival + i);

             }
          }

       }
       else {
          for(i=0; i<dnum; i++) {scanf("%f", ival+i); tele_sum += 1.0; ittmp += *(ival+i);}
       }

       alltr = read_tracks(fin, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);


       if(!trc){

          if(aniso == 'y') fprintf(fout, "%d\n", 1);
          else fprintf(fout, "%d\n", 0);

          plw = ftell(fout);

          fprintf(fout, "WT_INFO %1d %12.5f\n", 1, tele_sum);

          if(iwrt){


             fprintf(fout, "IND_WTS %d %d\n", nf, dnum);
             plww = ftell(fout);
             for(i=0; i < nf; i++) {
                for(j=0; j < dnum; j++){
                   fprintf(fout, "%8.4f ", *ival);
                }
                fprintf(fout, "\n");
             }

          }

          fprintf(fout, "%d %d\n", gpr, ipr);

          if(gpr) fprintf(fout, "%f %f\n", alat, alng);

          pl = ftell(fout);

          fprintf(fout, "TRACK_NUM  %8d ADD_FLD  %3d %3d &", tr_count, nff, nfld);
          if(nfwpos){
             for(i=0; i < nff; i++)fprintf(fout, "%1d", *(nfwpos + i));
          }
          fprintf(fout, "\n");

       }


       trnr += meantrd_tanh(fout, alltr , tr_count, 's', istr, iend, ival, dnum, iraw, ireg, lat1, lat2, lng1, lng2);

       trnum += tr_count;

       for(i=0; i < tr_count; i++) free((alltr+i)->trpt);

       free(alltr);

       ++trc;

       fclose(fin);
       fin=NULL;




    } 

    printf("%d %d\n", trnum, trnr);

    printf("%d\n", ittmp);

    if(ireg) trnum = trnr;


    fseek(fout, pl, FSTART);
    fprintf(fout, "TRACK_NUM  %8d", trnum);

    fseek(fout, plw, FSTART);
    if(iraw) 
      fprintf(fout, "WT_INFO %1d %12.5f\n", 1, tele_sum);
    else
      fprintf(fout, "WT_INFO %1d %12.5f\n", 2, tele_sum);

    if(iwrt){
       fseek(fout, plww, FSTART);
       for(i=0; i < nf; i++){
           for(j=0; j < dnum; j++){
              fprintf(fout, "%8.4f ", *(wval + i * dnum + j));
           }
           fprintf(fout, "\n");
       } 


    } 

    free(wval);

    fclose(fout);



    return 0;

}
