#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "splice.h"

/* function to iterrogate the location data attached to tracks */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;


int main(int argc, char **argv)
{
    int i, j;
    int trnum=0;
    int gpr=0, ipr=0;
    int ftype='s';

    int ip1=0, ip2=0;
    int iff1=0, iff2=0, iff3=0;
    int idd=0;
    int nth=0;
    int timthresh=0;

    float alat=0.0, alng=0.0;
    float diff=0.0;
    float dmax=0.0;
    float dist=0.0;
    float pthresh=0.0, dthresh=0.0;

    FILE *trin=NULL;
    FILE *fout=NULL;

    char trfil[MAXCHR];
    char trout[MAXCHR];

    struct tot_tr *trr=NULL;
    struct tot_tr *atr=NULL;
    struct fet_pt_tr *fp=NULL;

    if(argc < 9){
       printf("Usage: location [input file] [output file] [1st loc] [2nd loc] [dist loc] [pthress] [distthresh] [timthresh]\n");
       printf("1st loc = Id of first location point.                   \r\n"
              "2nd loc = Id of second location point.                  \r\n"
              "dist loc = Id of location for distance check.           \r\n"
              "pthresh = threshold for pressure difference             \r\n"
              "distthresh = threshold for distance from target point   \r\n"
              "timthresh = number of time steps that conditions prevail\r\n"
              "Difference criteria loc2 - loc1                         \n\n");
       exit(1);
    }

    sscanf(argv[1], "%s", trfil);
    sscanf(argv[2], "%s", trout);
    sscanf(argv[3], "%d", &ip1);
    sscanf(argv[4], "%d", &ip2);
    sscanf(argv[5], "%d", &idd);
    sscanf(argv[6], "%f", &pthresh);
    sscanf(argv[7], "%f", &dthresh);
    sscanf(argv[8], "%d", &timthresh);

    trin = fopen(trfil, "r");
    if(!trin){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", trfil);
       exit(1);
    }

    trr = read_tracks(trin, &trnum, &gpr, &ipr, ftype, &alat, &alng, NULL, NULL);

    if(ip1 < 0 || ip1 > nff) {
       printf("****ERROR***, requested default or additional field does not exist for track file.\n\n");
       exit(1);
    }

    if(ip2 < 0 || ip2 > nff) {
       printf("****ERROR***, requested default or additional field does not exist for track file.\n\n");
       exit(1);
    }

    if(idd < 0 || idd > nff) {
       printf("****ERROR***, requested default or additional field does not exist for track file.\n\n");
       exit(1);
    }

    if(ip1){
      iff1 = 0;
      for(i=0; i < ip1; i++){
         if(*(nfwpos + i)) iff1 += 3;
         else iff1 += 1;
      }
      --iff1;
    }

    if(ip2){
      iff2 = 0;
      for(i=0; i < ip2; i++){
         if(*(nfwpos + i)) iff2 += 3;
         else iff2 += 1;
      }
      --iff2;
    }

    if(idd){
      iff3 = 0;
      for(i=0; i < idd; i++){
         if(*(nfwpos + i)) iff3 += 3;
         else iff3 += 1;
      }
      --iff3;
    }

    printf("%d %d %d\n", iff1, iff2, iff3);

    for(i=0; i < trnum; i++){

        atr = trr + i;

        if(atr->num > 0){

           dmax = 0.0;

           nth = 0;

           for(j=0; j < atr->num; j++){

               fp = atr->trpt + j;

               diff = fp->add_fld[iff2] - fp->add_fld[iff1];
               if(diff > dmax) {dmax = diff; dist = fp->add_fld[iff3];}
               if(diff >= pthresh && dist <= dthresh) ++nth;
           }

           if(dmax < pthresh || dist > dthresh || nth < timthresh) {atr->num = 0; free(atr->trpt);}
           else
              printf("pdif = %f, dist = %f\n", dmax, dist);

        } 


    }

    fout = fopen(trout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", trout);
       exit(1);
    }

    meantrd(fout, trr, trnum, ftype, gpr, ipr, alat, alng);

    fclose(fout);


    return 0;
}
