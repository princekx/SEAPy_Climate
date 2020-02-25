#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"

#define TOLMATCH    2.0
#define DISTLARGE  180.0
#define TOLNUM     0.6
#define MAXCHR     50
#define FILENAME    "match_"
#define ID1        -10.0
#define ID2        10.0
#define NID        40
#define IDD1       0.
#define IDD2       2.
#define NIDD       10
#define TPBW       5.0


/* program to compare two track ensembles */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
int toverlap(struct tot_tr * , struct tot_tr * , long int * , long int * );
float trdist(struct tot_tr * , struct tot_tr * , long int , long int , int * , int * , int * , int , int , int , int , int );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;

int main(int argc, char **argv)
{

    int i,j;
    int nump=0, numpm=0;
    int totnum=0;
    int jtrmin=0;
    int nmat=0;

    int *nm1=NULL;

    int trnum1, trnum2;
    int gpr1, ipr1, gpr2, ipr2;
    int nbin=0;
    int iper_num_t;
    int str_type=0, sep_type=0;
    int nmbin=0;
    int nid=NID, nidd=NIDD;
    int in1, in2, it1, it2;
    int ior=0, iorr=0;


    long int is1, is2;

    char filout[MAXCHR];

/* Histograms:- 1st match, 1st no match, 2nd match, 2nd no match */

    float *hist[4]={NULL, NULL, NULL, NULL};
    float *nptm=NULL;
    float *npid=NULL;
    float *npdd=NULL;
    float **ntpdt=NULL;

    float id1=ID1, id2=ID2;
    float idd1=IDD1, idd2=IDD2;

    float alat1=0.0, alng1=0.0, alat2=0.0, alng2=0.0;
    float dist=0.0;
    float distmin;
    float tolnum=TOLNUM, tolmatch=TOLMATCH;
    float trpm1=0.0, trpm2=0.0;

    float frat=0.0, ffrat=0.0;
    float ditty;
    float t1, t2, dt, did, ddd;
    float tpbw=TPBW;

    float *temp=NULL;
    float *temp2=NULL;
    float *temp3=NULL;
    float *temp4=NULL;

    float str1, str2, str;
    float ntrm=0.0;
    float nmatch=0.0;

    float sum_per_t;

    FILE *fin1=NULL, *fin2=NULL;
    FILE *fout=NULL;

    FILE *fcmp=NULL;

    FILE *ftpdt=NULL;

    struct tot_tr *tr1=NULL, *tr2=NULL;
    struct tot_tr *tm1=NULL, *tm2=NULL;
    struct tot_tr *atr1=NULL, *atr2=NULL;

    struct fet_pt_tr *fp1=NULL, *fp2=NULL;

    if(argc < 7){
       printf("Usage: ensemble [track file1] [track file2] [temp1] [temp2] [nbin]               \r\n"
              "                [mean (0) or max (1) intensity] [Mean (0)  or Min (1) seperation]\r\n"
              "                [Dist. Thresh. (optional)] [Time. Thresh. (optional)]            \n\n");
       exit(1);
    }


    fcmp = fopen("cmp.dat", "r");
    if(fcmp){
       printf("****WARNING****, file exists with subsidary bining data for \r\n"
              "                 intensity and seperation distributions     \r\n"
              "                 this will be used instead of defaults.     \n\n");
       fscanf(fcmp, "%f", &tpbw);
       fscanf(fcmp, "%f %f %d", &id1, &id2, &nid);
       fscanf(fcmp, "%f %f %d", &idd1, &idd2, &nidd);

    }
    
    fin1 = fopen(argv[1], "r");
    if(!fin1){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[1]);
       exit(1);
    }

    tr1 = read_tracks(fin1, &trnum1, &gpr1, &ipr1, 's', &alat1, &alng1, NULL, NULL);
    fseek(fin1, 0L, SEEK_SET);
    tm1 = read_tracks(fin1, &trnum1, &gpr1, &ipr1, 's', &alat1, &alng1, NULL, NULL);

    iper_num_t = iper_num;
    sum_per_t = sum_per;

    fclose(fin1);

    fin2 = fopen(argv[2], "r");
    if(!fin2){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[2]);
       exit(1);
    }

    tr2 = read_tracks(fin2, &trnum2, &gpr2, &ipr2, 's', &alat2, &alng2, NULL, NULL);
    fseek(fin2, 0L, SEEK_SET);
    tm2 = read_tracks(fin2, &trnum2, &gpr2, &ipr2, 's', &alat2, &alng2, NULL, NULL);

    fclose(fin2);


    sscanf(argv[3], "%f", &t1);
    sscanf(argv[4], "%f", &t2);
    sscanf(argv[5], "%d", &nbin);
    sscanf(argv[6], "%d", &str_type);
    sscanf(argv[7], "%d", &sep_type);

    if(argc > 8){
       if(argc != 10){
          printf("***ERROR***, must specify both TOLMATCH and TOLNUM\n\n");
          exit(1);
       }
       sscanf(argv[8], "%f", &tolmatch);
       sscanf(argv[9], "%f", &tolnum);
    }

    if(str_type < 0 || str_type > 1){
       printf("****ERROR****, wrong intensity ID\n\n");
       exit(1);

    }

    temp = (float *)calloc(nbin, sizeof(float));
    mem_er((temp == NULL) ? 0 : 1, nbin * sizeof(float));


    nmbin = (int)floor((100.0 * (1.0 - tolnum) / tpbw) + 0.5);

    for(i=0; i < 4; i++){
        hist[i] = (float *)calloc(nbin + 1, sizeof(float));
        mem_er((hist[i] == NULL) ? 0 : 1, (nbin + 1) * sizeof(float));
    }

    temp2 = (float *)calloc(nmbin, sizeof(float));
    mem_er((temp2 == NULL) ? 0 : 1, nmbin * sizeof(float));

    temp3 = (float *)calloc(nid, sizeof(float));
    mem_er((temp3 == NULL) ? 0 : 1, nid * sizeof(float));

    temp4 = (float *)calloc(nidd, sizeof(float));
    mem_er((temp4 == NULL) ? 0 : 1, nidd * sizeof(float));

    nptm = (float *)calloc(nmbin, sizeof(float));
    mem_er((nptm == NULL) ? 0 : 1, nmbin * sizeof(float));

    npid = (float *)calloc(nid, sizeof(float));
    mem_er((npid == NULL) ? 0 : 1, nid * sizeof(float));

    npdd = (float *)calloc(nidd, sizeof(float));
    mem_er((npdd == NULL) ? 0 : 1, nidd * sizeof(float));

    ntpdt = (float **)calloc(nidd, sizeof(float *));
    mem_er((ntpdt == NULL) ? 0 : 1, nidd * sizeof(float *));

    for(i=0; i < nidd; i++){
        *(ntpdt + i) = (float *)calloc(nmbin, sizeof(float));
        mem_er((*(ntpdt + i) == NULL) ? 0 : 1, nmbin * sizeof(float));        
    }
    

/* assign temperature bins */

    dt = (t2 - t1) / nbin;
    did = (id2 - id1) / nid;
    ddd = (idd2 - idd1) / nidd;

    for(i=0; i < nbin; i++) *(temp + i) = t1 + 0.5 * dt + i * dt; 


    nm1 = (int *)calloc(trnum1, sizeof(int));
    mem_er((nm1 == NULL) ? 0 : 1, trnum1 * sizeof(int));

/* assign bins for the percentage of points that match */

    tpbw /= 100.0;

    for(i=0; i < nmbin; i++) *(temp2 + i) = tolnum + tpbw * ((float)i + 0.5);

    for(i=0; i < nid; i++) *(temp3 + i) = id1 + 0.5 * did + i * did;
    
    for(i=0; i < nidd; i++) *(temp4 + i) = idd1 + 0.5 * ddd + i * ddd;

    if((gpr1 != gpr2 || ipr1 != ipr2) ||
        (fabs(alat1 - alat2) > 1.0e-4 || fabs(alng1 - alng2) > 1.0e-4)){
        printf("***ERROR***, projection parameters do not match\n\n");
        exit(1);
    }

    if((tr1->time && !tr2->time) || (!tr1->time && tr2->time)){
       printf("***ERROR***, temporal identities do not match\n\n");
       exit(1);
    }

    if(iper_num != iper_num_t || fabs(sum_per - sum_per_t) > 1.0e-6){
       printf("****ERROR****, number of time periods do not match\n\n");
       exit(1);

    }


    for(i=0; i < trnum1; i++){

        atr1 = tr1 + i;


        if(!str_type){

          str1 = 0.0;

          for(j=0; j < atr1->num; j++) str1 += (atr1->trpt + j)->zf;

          str1 /= atr1->num;

        }

        else if(str_type == 1){

          str1 = atr1->trpt->zf;
          for(j=0; j < atr1->num; j++) {
              str = (atr1->trpt + j)->zf;
              if(str > str1) str1 = str;
          }

        }

        if(str1 > t2){
          printf("***ERROR***, mean strength outside range, %f.\n\n", str1);
          exit(1);
        }

        distmin = DISTLARGE;
        jtrmin = -1;
        numpm = 0;
        frat = 0.0;

        nmat = 0;

        if(atr1->num){

           for(j=0 ; j < trnum2; j++){

              atr2 = tr2 + j;

              if(atr2->num){

/* Check for overlap in time */

                 if(toverlap(atr1, atr2, &is1, &is2)){

                    dist = trdist(atr1, atr2, is1, is2, &nump, &it1, &it2, sep_type, 0, 0, -1, -1);
                    ffrat = 2.0 * nump / (float)(atr1->num + atr2->num);

                    if(dist < tolmatch && ffrat >= tolnum){

                       ++nmat;

                       if(jtrmin < 0 ){

                          jtrmin = j;
                          distmin = (sep_type == 0) ? dist : trdist(atr1, atr2, is1, is2, &nump, &it1, &it2, 0, 0, 0, -1, -1);
                          numpm = nump;
                          in1 = it1;
                          in2 = it2;

                       }

                       else if(jtrmin >= 0 && nump > numpm) {

                            jtrmin = j;
                            distmin = (sep_type == 0) ? dist : trdist(atr1, atr2, is1, is2, &nump, &it1, &it2, 0, 0, 0, -1, -1);;
                            numpm = nump;
                            in1 = it1;
                            in2 = it2;


                       }


                    }


                 }


              }

           }

           if(jtrmin >= 0) frat = 2.0 * numpm / (float)(atr1->num + (tr2 + jtrmin)->num);

           if(jtrmin >= 0 && frat >= tolnum){

              *(nm1 + i) = nmat;
              nmatch += 1.0;

              if(!str_type){

                 str2 = 0.0;

                 for(j=0; j < (tr2 + jtrmin)->num; j++) str2 += ((tr2 + jtrmin)->trpt + j)->zf;

                 str2 /= (tr2 + jtrmin)->num;

              }

              else if(str_type == 1){

                 str2 = (tr2 + jtrmin)->trpt->zf;
                 for(j=0; j < (tr2 + jtrmin)->num; j++) {
                     str = ((tr2 + jtrmin)->trpt + j)->zf;
                     if(str > str2) str2 = str;
                 }

              }

              if(str2 > t2){
                printf("***ERROR***, mean strength outside range %f.\n\n", str2);
                exit(1);
              }

              ++(*(hist[0] + (int)((str1 - t1) / dt)));
              ++(*(hist[2] + (int)((str2 - t1) / dt)));
              ++(*(nptm + (int)((frat - tolnum) / tpbw)));
              atr1->num = 0;
              (tr2 + jtrmin)->num = 0;

              totnum += numpm;

              for(j=0; j < numpm; j++){

                  fp1 = atr1->trpt + in1 + j;
                  fp2 = (tr2 + jtrmin)->trpt + in2 + j;
                  ditty = fp2->zf - fp1->zf;
                  if(ditty < id1 || ditty > id2){
                     if(!ior){
                        printf("****WARNING****, track point intensity difference out of range for histogram.\n\n");

                        ior = 1;
                     }
                     printf("Value = %f\n", ditty);
                  }

                  else  ++(*(npid + (int)((ditty - id1) / did)));

              }

              if(distmin > idd2){

                 if(!iorr){
                    printf("****WARNING****, track point seperation distance out of range for histogram.\n\n");

                    iorr = 1;
                 }
                     printf("Value = %f\n", distmin);
              }

              else  {
                   ++(*(npdd + (int)((distmin - idd1) / ddd)));
                   ++(*(*(ntpdt + (int)((distmin - idd1) / ddd)) + (int)((frat - tolnum) / tpbw)));
              }


              
           }

           else {

              ++(*(hist[1] + (int)((str1 - t1) / dt)));
              

           }

        }

    }

    for(i=0; i < trnum2; i++){

        atr2 = tr2 + i;

        if(atr2->num){

           if(!str_type){

              str2 = 0.0;

              for(j=0; j < atr2->num; j++) str2 += (atr2->trpt + j)->zf;

              str2 /= atr2->num;

            }

            else if(str_type == 1){

                str2 = atr2->trpt->zf;
                for(j=0; j < atr2->num; j++) {
                   str = (atr2->trpt + j)->zf;
                   if(str > str2) str2 = str;
                }

            }


            if(str2 > t2){
              printf("***ERROR***, mean strength outside range %f.\n\n", str2);
              exit(1);
            }
            ++(*(hist[3] + (int)((str2 - t1) / dt)));

        }

    }

    for(i=0; i < trnum1; i++){
        atr1 = tr1 + i;
        atr2 = tm1 + i;

        if(atr1->num > 0) atr2->num = 0;

    }

    for(i=0; i < trnum2; i++){
        atr1 = tr2 + i;
        atr2 = tm2 + i;

        if(atr1->num > 0) atr2->num = 0;

    }


/* write out sub-sets of track ensembles */

    strcpy(filout, FILENAME);
    strcat(filout, "ens1_yes.dat");
    
    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    meantrd(fout, tm1, trnum1, 's', gpr1, ipr1, alat1, alng1);

    fclose(fout);

    strcpy(filout, FILENAME);
    strcat(filout, "ens1_no.dat");
    
    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    meantrd(fout, tr1, trnum1, 's', gpr1, ipr1, alat1, alng1);

    fclose(fout);


    strcpy(filout, FILENAME);
    strcat(filout, "ens2_yes.dat");
    
    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    meantrd(fout, tm2, trnum2, 's', gpr2, ipr2, alat2, alng2);

    fclose(fout);

    strcpy(filout, FILENAME);
    strcat(filout, "ens2_no.dat");
    
    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    meantrd(fout, tr2, trnum2, 's', gpr2, ipr2, alat2, alng2);

    fclose(fout);



    for(i=0; i < trnum1; i++){
        if(*(nm1 + i) > 1) ntrm += 1.0;
    }


    if(iper_num){

       ntrm /= sum_per;


       for(i=0; i < nbin; i++){
           for(j=0; j < 4; j++)
               *(hist[j] + i) /= sum_per;
       }

/*       for(i=0; i < nmbin; i++)
           *(nptm + i) /= sum_per; */

/*       for(i=0; i < nid; i++)
           *(npid + i) /= nmatch; */

   

    }

    for(i=0; i < nbin; i++){
        for(j=0; j < 4; j++)
            *(hist[j] + nbin) += *(hist[j] + i);
    }


    trpm1 = (float)trnum1 / sum_per;
    trpm2 = (float)trnum2 / sum_per;

    printf("************ SUMMARY STATISTICS ********************\n\n");

    printf("First file %s, second file %s\n\n", argv[1], argv[2]);
    printf("Distance threshold is %f deg.\n", tolmatch);
    printf("Relative No. of points threshold is %f\n", tolnum);
    printf("No. of periods is %f \n", sum_per);
    printf("No. of tracks in ensemble 1 is %d\n", trnum1);
    printf("No. of tracks per month in ensemble 1 is %f\n", trpm1); 
    printf("No. of tracks in ensemble 2 is %d\n", trnum2);
    printf("No. of tracks per month in ensemble 1 is %f\n", trpm2);
    printf("No of tracks that match for distance and No. of points thresholds is %f\n", nmatch / sum_per);
    printf("No. of tracks in ensemble 1 that match more than one track \r\n"
           "in ensemble 2 for distance and No. of points thresholds is %f\n\n", ntrm);

    printf("Percentage of tracks that match in ensemble 1 is %f\n", *(hist[0] + nbin) / trpm1);
    printf("Percentage of tracks that don't match in ensemble 1 is %f\n", *(hist[1] + nbin) / trpm1);

    printf("Percentage of tracks that match in ensemble 2 is %f\n", *(hist[2] + nbin) / trpm2);
    printf("Percentage of tracks that don't match in ensemble 2 is %f\n", *(hist[3] + nbin) / trpm2);

    printf("*****************Histograms**************************\n\n");
    printf("          1st Match   1st No Match   2nd Match   2nd No Match\n\n");
    for(i=0; i < nbin; i++){
        printf("%9.6f    %9.4f          %9.4f        %9.4f           %9.4f\n", *(temp + i),  *(hist[0] + i), *(hist[1] + i), *(hist[2] + i), *(hist[3] + i));
    }
    printf("Sum          %9.4f          %9.4f        %9.4f           %9.4f\n",   *(hist[0] + nbin), *(hist[1] + nbin), *(hist[2] + nbin), *(hist[3] + nbin));

    printf("-----------------------------------------------------\n");
    printf("      Number of Points That Match Distribution       \n");
    printf("Fract.  ");
    for(i=0; i<nmbin; i++)printf("%6.3f ", *(temp2 + i));
    printf("\n");
    printf("No.     ");
    for(i=0; i<nmbin; i++)printf("%6.2f ", *(nptm + i) / (nmatch * tpbw));
    printf("\n");

    printf("\n\n");

    printf("-----------------------------------------------------\n");
    printf("      Distribution for intensity differences         \n");
    printf("Intens.  ");
    for(i=0; i<nid; i++)printf("%8.4f ", *(temp3 + i));
    printf("\n");
    printf("No.      ");
    for(i=0; i<nid; i++)printf("%8.5f ", *(npid + i) / (totnum * did));

    printf("\n\n");

    printf("\n\n");

    printf("-----------------------------------------------------\n");
    printf("      Distribution for seperation distances          \n");
    printf("Seperat..  ");
    for(i=0; i<nidd; i++)printf("%8.4f ", *(temp4 + i));
    printf("\n");
    printf("No.      ");
    for(i=0; i<nidd; i++)printf("%8.5f ", *(npdd + i) / (nmatch * ddd));

    printf("\n\n");

/* write 2D statistics */

    ftpdt = fopen("temp-dist.stat", "w");

/*     fprintf(ftpdt, "%d %d\n", nidd, nmbin); */
    for(i=0; i < nidd; i++){

/*        for(j=0; j < nmbin; j++) fprintf(ftpdt, "%8.4f ", *(*(ntpdt + i) + j) / (nmatch * nmbin * ddd * tpbw));
        fprintf(ftpdt, "\n"); */
          for(j=0; j < nmbin; j++) 
             fprintf(ftpdt, "%8.4f %8.4f %8.4f \n", *(temp2 + j), *(temp4 + i), *(*(ntpdt + i) + j) / (nmatch * nmbin * ddd * tpbw));
        

    }

    fclose(ftpdt);

    return 0;

}

