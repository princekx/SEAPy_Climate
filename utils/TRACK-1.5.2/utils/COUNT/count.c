#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "splice.h"

#define  FP_PI   0.017453292519943295   
#define  FP_PI2  1.570796326794896619

#define  LARGE_DIST  1.0e+12
#define  NTIM        100000


struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
long int new_time(long int , int );
void sincos(double , double * , double * );

int noheader=0;

int main(int argc, char *argv[])

{

    int i,j;
    int jj=0, iwt=0;
/*    int year; */
    int gpr, ipr, trnum;
    int inn=0;
    int time=0;
    int negate=0;
    int itim=0;
    int tstep=6;

    int igen=0;

    long int ntime[NTIM], ttt;

    float alat, alng, rad;
    float ssum=LARGE_DIST;
    float wthresh=0.0;

    double x1, y1, z1;
    double x2, y2, z2;
    double s1, c1, s2, c2;
    double con, sum;
    

    char newf[MAXCHR];
    char com[]="count [filname] [Lat.] [Lng.] [Rad.] [Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)] [Negate (1)] [Start Time, YYYYMMDDHH] [tstep]";

    int ntr=0;

    struct tot_tr *all_tr=NULL, *atr=NULL;
    struct fet_pt_tr *pt=NULL;

    FILE *fin=NULL;
    FILE *fout=NULL;
    FILE *ftime=NULL;


    if(argc < 7){

       printf("%s \n", com);
       exit(1);


    }

    sscanf(argv[2], "%f", &alat);
    sscanf(argv[3], "%f", &alng);
    sscanf(argv[4], "%f", &rad);
    sscanf(argv[5], "%d", &igen);
    sscanf(argv[6], "%d", &negate);

    con = 1.0 / cos(rad * FP_PI);

    alng = alng * FP_PI;
    alat = FP_PI2 - alat * FP_PI;

    sincos(alng, &s1, &c1);
    sincos(alat, &s2, &c2);

    x1 = s2 * c1;
    y1 = s2 * s1;
    z1 = c2;

/*    sscanf(argv[1], "%*14c%d", &year); */

    fin = fopen(argv[1], "r");
    if(!fin){

       printf("No such file exists\n");
       exit(1);


    }

    all_tr = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    if(all_tr->awt){
       printf("****INFORMATION****, tracks have weights associated with them.\n\n");
       printf("Do you want to use the weights as a means to exclude tracks, 'y' or 'n'\n\n");
       scanf("\n");
       if(getchar() == 'y'){
          printf("****WARNING****, this only applies to options 0, 1, and 2.\n\n");
          printf("What threshold is required for the weights?\n\n");
          scanf("%f", &wthresh);
          iwt = 1;
       }

    }

    fclose(fin);

    if(igen == 3 || igen == 4){

       if(argc >= 8){

          sscanf(argv[7], "%ld", &ttt);
          printf("%ld\n", ttt);
          ntime[0] = ttt;
          sscanf(argv[8], "%d", &tstep);
printf("%d\n", tstep);
          for(i=1; i < NTIM; i++) ntime[i] = new_time(ntime[i-1], tstep);
       }

       else {

          printf("****ERROR****, to use actual times require a starting time\n\n");
          printf("%s \n", com);
          exit(1);

       }

    }


    if(igen == 4){

       for(i=0; i < trnum; i++){

          atr = all_tr + i;
          atr->time = ntime[atr->trpt->fr_id - 1];
          itim = atr->trpt->fr_id - 1;

          for(j=0; j< atr->num; j++){

              pt = atr->trpt + j;
              pt->time = ntime[itim];
              itim += 1 + pt->nfm;
          }

       }

    }

    else if(igen == 3){

       ftime = fopen("time.dat", "a");

       for(i=0; i < trnum; i++){

          atr = all_tr + i;

          inn = 0;

          ssum = LARGE_DIST;


          for(j=0; j< atr->num; j++){

              pt = atr->trpt + j;


              alat = FP_PI2 - pt->yf * FP_PI;
              alng = pt->xf * FP_PI;


              sincos(alng, &s1, &c1);
              sincos(alat, &s2, &c2);

              x2 = s2 * c1;
              y2 = s2 * s1;
              z2 = c2;
  
              sum = con * (x1 * x2 + y1 * y2 + z1 * z2) - 1.0;


              if(sum > 0.){

                inn=1;
                if(sum < ssum) {ssum = sum; time = pt->fr_id;}
              }


           }

           if(inn){
              if(negate) atr->num = 0;
              else {
                atr->time = ntime[time-1];
                fprintf(ftime, "%ld\n", ntime[time-1]);
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


          for(j=0; j< atr->num; j++){

              pt = atr->trpt + j;


              alat = FP_PI2 - pt->yf * FP_PI;
              alng = pt->xf * FP_PI;


              sincos(alng, &s1, &c1);
              sincos(alat, &s2, &c2);

              x2 = s2 * c1;
              y2 = s2 * s1;
              z2 = c2;
  
              sum = con * (x1 * x2 + y1 * y2 + z1 * z2) - 1.0;

              if(sum >= 0.){inn=1; jj = j; break;}

           }

           if(inn && negate) atr->num = 0;

           else if(!inn && !negate) atr->num = 0;

           if(iwt){
              pt = atr->trpt + jj;
              if(fabs(pt->wght) < wthresh) atr->num = 0;
           }


       }

    }

    else {

       for(i=0; i < trnum; i++){

          atr = all_tr + i;


          if(!igen) pt = atr->trpt;
          else pt = atr->trpt + atr->num - 1;

          alat = FP_PI2 - pt->yf * FP_PI;
          alng = pt->xf * FP_PI;


          sincos(alng, &s1, &c1);
          sincos(alat, &s2, &c2);

          x2 = s2 * c1;
          y2 = s2 * s1;
          z2 = c2;
  
          sum = con * (x1 * x2 + y1 * y2 + z1 * z2) - 1.0;

          if(sum >= 0. && negate) atr->num=0;
          else if(sum < 0. && !negate) atr->num=0;

          if(iwt){
             if(fabs(pt->wght) < wthresh) atr->num = 0;
          }


       }

    }


/*    printf("%d %d\n", year, ntr); */

    strcpy(newf, argv[1]);
    strcat(newf, ".new");

    fout = fopen(newf, "w");

    meantrd(fout, all_tr, trnum, 's', gpr, ipr, alat, alng);

    fclose(fout);


    return ntr;

}
