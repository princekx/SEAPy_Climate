#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "splice.h"
#include "mem_er.h"

#define  FP_PI   0.017453292519943295   
#define  FP_PI2  1.570796326794896619

#define  LARGE_DIST  1.0e+12

#define BOXTOL   0.00

#define ADD_H    1.0e+10


struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
long int new_time(long int );

int noheader=0;
extern int nfld;

int main(int argc, char *argv[])

{

    int i,j, k;
    int year;
    int gpr, ipr, trnum;
    int inn=0;
    int time=0;
    int negate=0;
    int ig=0;
    int nbox=0;
    int ihist=0;
    int nhbin=0;
    int htype=2;
    int nval=0;
    int hoty=0;

    int igen=0;

    long int ntime[1000], ttt;

    float alat, alng;
    float alath, alngh;
    float alat1, alat2, alng1, alng2;
    float ssum=LARGE_DIST;
    float tr_thrs=BOXTOL;
    float maxstr;
    float bmin, bmax, bint;
    float *hist=NULL, *binv=NULL;
    float *add=NULL;

    double sum1, sum2, sum;
    double d1, d2;
    

    char newf[MAXCHR];
    char com[]="count [filname] [Lat1.] [Lat2.] [Lng1.] [Lng2.] [Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)] [Negate (1)] [Track Thresh] [Start Time, YYYYMMDDHH]";

    int ntr=0;

    struct tot_tr *all_tr=NULL, *atr=NULL;
    struct fet_pt_tr *pt=NULL;

    FILE *fin=NULL;
    FILE *fout=NULL;
    FILE *ftime=NULL;


    if(argc < 8){

       printf("%s \n", com);
       exit(1);


    }

    sscanf(argv[2], "%f", &alat1);
    sscanf(argv[3], "%f", &alat2);
    sscanf(argv[4], "%f", &alng1);
    sscanf(argv[5], "%f", &alng2);
    sscanf(argv[6], "%d", &igen);
    sscanf(argv[7], "%d", &negate);

    alath = 0.5 * (alat1 + alat2);

    if(alng2 < alng1){

       printf("****WARNING****, stradeling the Greenwich meridian\n\n");
       ig = 1;
       alngh = 0.5 * (alng2 + alng1 - 360.0);
       if(alngh < 0.0) {alngh += 360.0; ig = 2;}


    }

    else alngh = 0.5 * (alng1 + alng2);


    fin = fopen(argv[1], "r");
    if(!fin){

       printf("No such file exists\n");
       exit(1);


    }

    if(argc >= 9) sscanf(argv[8], "%f", &tr_thrs);

    printf("Track Threshold = %f\n", tr_thrs);

    if(igen == 3 || igen == 4){

       if(argc == 10){

          sscanf(argv[9], "%ld", &ttt);
          printf("%ld\n", ttt);
          ntime[0] = ttt;
          for(i=1; i < 1000; i++)ntime[i] = new_time(ntime[i-1]);
       }

       else{

          printf("****ERROR****, to use actual times require a starting time\n\n");
          printf("%s \n", com);
          exit(1);

       }


    }

    all_tr = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);

    if(igen == 4){

       for(i=0; i < trnum; i++){

          atr = all_tr + i;
          atr->time = ntime[atr->trpt->fr_id - 1];

          for(j=0; j< atr->num; j++){

              pt = atr->trpt + j;
              pt->time = ntime[pt->fr_id - 1];
          }

       }

    }

    else if(igen == 3){

       ftime = fopen("time.dat", "a");

       for(i=0; i < trnum; i++){

          atr = all_tr + i;

          inn = 0;

          ssum = LARGE_DIST;

          maxstr = 0.0;


          for(j=0; j< atr->num; j++){

              pt = atr->trpt + j;


              alat = pt->yf;
              alng = pt->xf;

              sum1 = (alat2 - alat) * (alat - alat1);
              d1 = alath - alat;
              if(!ig) {
                 sum2 = (alng2 - alng) * (alng - alng1);
                 d2 = alngh - alng;

              }
              else{
                 sum2 = (360.0 - alng) * (alng - alng1);

                 d2 = (ig == 1) ? alngh + 360.0 - alng : alngh - alng; 

                 if(sum2 < 0.0) {
                    sum2 = (alng2 - alng) * alng;
                    d2 = (ig == 2) ? alngh - alng - 360.0 : alngh - alng;

                 }

              }
              


              if(sum1 >= 0. && sum2 >= 0.0){

                sum = sqrt(d1 * d1 + d2 * d2);

                inn=1;
                if(sum < ssum) {ssum = sum; time = pt->fr_id;}

              }


           }


           if(inn){
              if(negate) atr->num = 0;
              else {
                atr->time = ntime[time-1];

                for(j=0; j< atr->num; j++){

                    pt = atr->trpt + j;

                    if(pt->zf > maxstr) maxstr = pt->zf;

                }

                fprintf(ftime, "%10ld %10.4e\n", ntime[time-1], maxstr);
              }
           }

           else if(!inn) {
              if(!negate) atr->num = 0;
           }

        }

        fclose(ftime);

    }

    else if(igen == 2){

       

       for(i=0; i < trnum; i++){

          atr = all_tr + i;

          inn = 0;

          nbox = 0;

          for(j=0; j< atr->num; j++){

              pt = atr->trpt + j;


              alat = pt->yf;
              alng = pt->xf;

              sum1 = (alat2 - alat) * (alat - alat1);

              if(!ig) {
                 sum2 = (alng2 - alng) * (alng - alng1);


              }
              else{
                 sum2 = (360.0 - alng) * (alng - alng1);

                 if(sum2 < 0.0) {
                    sum2 = (alng2 - alng) * alng;

                 }

              }

              if(sum1 >= 0.0 && sum2 >= 0.0){inn=1; ++nbox;}

           }

           if( ! (((float)nbox) / atr->num  > tr_thrs) ) atr->num = 0;

           if(inn && negate) atr->num = 0;

           else if(!inn && !negate) atr->num = 0;


       }

    }

    else {

       for(i=0; i < trnum; i++){

          atr = all_tr + i;



          if(!igen) pt = atr->trpt;
          else pt = atr->trpt + atr->num - 1;

          alat = pt->yf;
          alng = pt->xf;

          sum1 = (alat2 - alat) * (alat - alat1);

          if(!ig) {
             sum2 = (alng2 - alng) * (alng - alng1);


          }
          else{
             sum2 = (360.0 - alng) * (alng - alng1);

             if(sum2 < 0.0) {
                sum2 = (alng2 - alng) * alng;

             }

          }

          if(sum1 >= 0. && sum2 >= 0.0 && negate) atr->num=0;
          else if((sum1 < 0. || sum2 < 0.) && !negate) atr->num=0;



       }

    }


/*    printf("%d %d\n", year, ntr); */

    strcpy(newf, argv[1]);
    strcat(newf, ".new");

    fout = fopen(newf, "w");

    meantrd(fout, all_tr, trnum, 's', gpr, ipr, alat, alng);

    fclose(fout);

    if(nfld){
       printf("Additional field values are available do you want to produce histograms, 'y' or 'n'.\n\n");
       scanf("\n");
       if(getchar() == 'y') ihist = 1;

       if(ihist){

          printf("For genesis '0' \r\n"
                 "For lysis   '1' \r\n"
                 "For all     '2' \n\n");
          scanf("%d", &htype);
          if(htype < 0 || htype > 2){
             printf("****ERROR****, histogram type %d not known.\n\n", htype);
             exit(1);
          }

       }

    }


    if(ihist){

       for(i=0; i < nfld; i++){

           printf("Input number of bins and the histogram range for additional value %d.\n\n", i+1);
           scanf("%d %f %f", &nhbin, &bmin, &bmax);

           binv = (float *)calloc(nhbin, sizeof(float));
           mem_er((binv == NULL) ? 0 : 1, nhbin * sizeof(float));

           hist = (float *)calloc(nhbin, sizeof(float));
           mem_er((hist == NULL) ? 0 : 1, nhbin * sizeof(float));

           bint = (bmax - bmin) / nhbin;

           for(j=0; j < nhbin; j++) *(binv + j) = bmin + ((float) j + 0.5) * bint;  

           nval = 0;

           for(j=0; j < trnum; j++){

               atr = all_tr + j;

               if(atr->num){

                  if(!htype){
                     add = atr->trpt->add_fld + i;
                     if(*add < ADD_H){
                        if(*add < bmin || *add > bmax){
                           printf("Value %f outside of range.\n\n", *add); 
                        }
                        else {
                           *(hist + (int)((*add - bmin) / bint)) += 1.0; 
                           ++nval;
                        }

                     }

                  }

                  else if(htype == 1){
                     add = (atr->trpt + atr->num - 1)->add_fld + i;
                     if(*add < ADD_H){
                        if(*add < bmin || *add > bmax){
                           printf("Value %f outside of range.\n\n", *add); 
                        }
                        else {
                           *(hist + (int)((*add - bmin) / bint)) += 1.0; 
                           ++nval;
                        }

                     }

                  }

                  else if(htype == 2){

                      for(k=0; k < atr->num; k++){
                          add = (atr->trpt + k)->add_fld + i;

                          pt = atr->trpt + k;

                          alat = pt->yf;
                          alng = pt->xf;

                          sum1 = (alat2 - alat) * (alat - alat1);

                          if(!ig) {
                             sum2 = (alng2 - alng) * (alng - alng1);
                          }
                          else{
                             sum2 = (360.0 - alng) * (alng - alng1);

                             if(sum2 < 0.0) {
                                sum2 = (alng2 - alng) * alng;

                             }

                          } 

                          if(*add < ADD_H && sum1 >= 0. && sum2 >= 0.0){
                             if(*add < bmin || *add > bmax){
                               printf("Value %f outside of range.\n\n", *add); 
                             }
                             else {
                               *(hist + (int)((*add - bmin) / bint)) += 1.0; 
                               ++nval;

                             }

                          }

                      }

                  }

               }

           }

           printf("Use raw histograms,       '0'\r\n"
                  "Use frequency histograms, '1'\r\n"
                  "Use pdf's,                '2'\n\n");
           scanf("%d", &hoty);

           if(hoty == 1) {
              for(j=0; j < nhbin; j++) *(hist + j) /= nval;
           }
           else if(hoty == 2){
              for(j=0; j < nhbin; j++) *(hist + j) /= (nval * bint);
           }

           printf("Distribution for additional field %d\n", i+1);

           for(j=0; j < nhbin; j++){
               printf("%8.2f %10.4f\n", *(binv + j), *(hist + j));
           }

           free(binv);
           free(hist);

       }

    }


    return ntr;

}
