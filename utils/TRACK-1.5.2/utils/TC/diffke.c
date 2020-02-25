#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"

#define  GTOL  1.0e-4

/* program to compute the KE derivative from the spatially sampled U and V fields */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
FILE *open_radial_file(char * , int * , long int * , int * , int * , int * , int * , int * , int * , int * , int * );
void kecalc(float * , float * , float * , int );
double anomaly(float * , double * , int , int , int );

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int main(void )
{

    int i, j, k, n;

    int dtype=0;
    int trnum=0;
    int irtrn=0;
    int irnth=0, irnr=0;
    int irnf=0;
    int idir=0;
    int ndsmth=0, ifcnt=0, itpadd=0, iprojd=0;
    int irdim=0;
    int nwtr=0;
    int iptc=0;
    int ianom=0, irad=0;

    int gpr=0, ipr=0;

    int istrt=0, iend=0;

    long int iptnum=0;

    float alat=0.0, alng=0.0;
    float ke1, ke2;
    float arad=0.0;

    float *sradf1=NULL, *sradf2=NULL;
    float *ske1=NULL, *ske2=NULL;
    float *sdiff=NULL, *sdum=NULL;
    float *slng1=NULL, *slat1=NULL;
    float *slng2=NULL, *slat2=NULL;

    double *coslat=NULL;

    off_t place1, place2, blklen;

    FILE *fintr=NULL;
    FILE *frad1=NULL, *frad2=NULL;
    FILE *frado=NULL;

    char filtr[MAXCHR];
    char filrad[MAXCHR];
    char filout[MAXCHR];

    struct tot_tr *tracks=NULL, *atr=NULL;

    printf("What type of derivative is required, '0' for foreward difference or '1' for centered difference.?\n\n");
    scanf("%d", &dtype);

    if(dtype < 0 || dtype > 1){
       printf("****ERROR****, invalid option %d\n", dtype);
       exit(1);
    }

    printf("What is the track file to read?\n\n");
    scanf("%s", filtr);

    fintr = fopen(filtr, "r");
    if(!fintr){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filtr);
       exit(1);
    }

    tracks = read_tracks(fintr, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fintr);

    printf("What is the filename for the U radial field?\n\n");
    scanf("%s", filrad);

    frad1 = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);

    if(irtrn != trnum){
       printf("****ERROR****, track file and radial field data not compatable: different numbers of tracks.\n\n");
       exit(1);
    }

    printf("What is the filename for the V radial field?\n\n");
    scanf("%s", filrad);
    frad2 = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);

    printf("%ld\n", iptnum);

/* setup arrays and read file headers */

    irdim = irnr * irnth;

    sradf1 = (float *)calloc(irdim, sizeof(float));
    mem_er((sradf1 == NULL) ? 0 : 1, irdim * sizeof(float));
    slng1 = (float *)calloc(irnth, sizeof(float));
    mem_er((slng1 == NULL) ? 0 : 1, irnth * sizeof(float));
    slat1 = (float *)calloc(irnr, sizeof(float));
    mem_er((slat1 == NULL) ? 0 : 1, irnr * sizeof(float));

    fread(slng1, irnth*sizeof(float), 1, frad1);
    fscanf(frad1, "%*c");
    fread(slat1, irnr*sizeof(float), 1, frad1);
    fscanf(frad1, "%*c");

    sradf2 = (float *)calloc(irdim, sizeof(float));
    mem_er((sradf2 == NULL) ? 0 : 1, irdim * sizeof(float));
    slng2 = (float *)calloc(irnth, sizeof(float));
    mem_er((slng2 == NULL) ? 0 : 1, irnth * sizeof(float));
    slat2 = (float *)calloc(irnr, sizeof(float));
    mem_er((slat2 == NULL) ? 0 : 1, irnr * sizeof(float));

    fread(slng2, irnth*sizeof(float), 1, frad2);
    fscanf(frad2, "%*c");
    fread(slat2, irnr*sizeof(float), 1, frad2);
    fscanf(frad2, "%*c");

/* KE arrays */

    ske1 = (float *)calloc(irdim, sizeof(float));
    mem_er((ske1 == NULL) ? 0 : 1, irdim * sizeof(float));
    ske2 = (float *)calloc(irdim, sizeof(float));
    mem_er((ske2 == NULL) ? 0 : 1, irdim * sizeof(float))

/* output array */

    sdiff = (float *)calloc(irdim, sizeof(float));
    mem_er((sdiff == NULL) ? 0 : 1, irdim * sizeof(float));

/* missing data array */

    sdum = (float *)calloc(irdim, sizeof(float));
    mem_er((sdum == NULL) ? 0 : 1, irdim * sizeof(float));

    for(i=0; i < irdim; i++) *(sdum + i) = ADD_UNDEF;

    for(i=0; i < irnr; i++){
        if(fabs(*(slat1 + i) - *(slat2 + i)) > GTOL){
           printf("****ERROR****, grids do not match in radial direction\n\n");
           exit(1);
        }
    }

    for(i=0; i < irnth; i++){
       if(fabs(*(slng1 + i) - *(slng2 + i)) > GTOL){
          printf("****ERROR****, grids do not match in tangential direction\n\n");
          exit(1);
       }
    }

/* determine averging radius if needed */

    printf("Do you want anomaly by subtracting an area average? '0' for no '1' for yes.\n\n");
    scanf("%d", &ianom);
    if(ianom < 0 || ianom > 1){
       printf("****ERROR****, incorrect option %d specified.\n\n", ianom);
       exit(1);
    }

    if(ianom){

       coslat = (double *)calloc(irnr, sizeof(double));
       mem_er((coslat == NULL) ? 0 : 1, irnr * sizeof(double));

       for(i=0; i < irnr; i++) *(coslat + i) = cos(*(slat1 + i) * FP_PI);

       printf("Specify the radius required for area averaging in degrees.\n\n");
       scanf("%f", &arad);
       arad = 90.0 - arad;
       irad = irnr;
       if(arad <= *(slat1 + irnr - 1)) {arad = *(slat1 + irnr - 1); irad = irnr;}
       else{
          for(i=0; i < irnr; i++){
              if(*(slat1 + i) < arad) {arad = *(slat1 + i - 1); irad = i; break;}
          }
       }
printf("%d %f\n", irad, arad);

    }
/* determine data block size */

    place1 = ftello(frad1);
    fread(sradf1, irdim*sizeof(float), 1, frad1);
    fscanf(frad1, "%*c");
    place2 = ftello(frad1);
    blklen = place2 - place1;
    fseeko(frad1, place1, SEEK_SET);

/* open output file and write header */

    printf("What is the name of the output file?\n\n");
    scanf("%s", filout);

    frado = fopen(filout, "w");
    if(!frado){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

/* write header */

    fprintf(frado, "%6d %10ld %3d %3d %2d %2d %3d %2d %2d %2d\n", irtrn, iptnum, irnth, irnr, irnf, idir, ndsmth, ifcnt, itpadd, iprojd);

    fwrite(slng1, irnth * sizeof(float), 1, frado);
    fprintf(frado, "\n");
    fwrite(slat1, irnr * sizeof(float), 1, frado);
    fprintf(frado, "\n");

    nwtr = 0;
    iptnum = 0;

    for(i=0; i < trnum; i++){

printf("%d\n", i);

       atr = tracks + i;

       iend = atr->num - 1;

       if(dtype) {istrt = 1;}
       else istrt = 0;

       iptc = 0;

       for(j=0; j < irnf; j++){
           fwrite(sdum, irdim * sizeof(float), 1, frado);
           fprintf(frado, "\n");
       }

       ++iptnum;

       for(j=istrt; j < iend; j++){

           for(k=0; k < irnf; k++){

              fseeko(frad1, ((nwtr + iptc) * irnf + k) * blklen + place1, SEEK_SET);
              fseeko(frad2, ((nwtr + iptc) * irnf + k) * blklen + place1, SEEK_SET);
              fread(sradf1, irdim*sizeof(float), 1, frad1);
              fscanf(frad1, "%*c");
              fread(sradf2, irdim*sizeof(float), 1, frad2);
              fscanf(frad2, "%*c");

              if(ianom){
                 anomaly(sradf1, coslat, irnr, irnth, irad);
                 anomaly(sradf2, coslat, irnr, irnth, irad);
              }

              kecalc(sradf1, sradf2, ske1, irdim);

              if(dtype){
                 fseeko(frad1, ((nwtr + iptc + 2) * irnf + k) * blklen + place1, SEEK_SET);
                 fseeko(frad2, ((nwtr + iptc + 2) * irnf + k) * blklen + place1, SEEK_SET);
                 fread(sradf1, irdim*sizeof(float), 1, frad1);
                 fscanf(frad1, "%*c");
                 fread(sradf2, irdim*sizeof(float), 1, frad2);
                 fscanf(frad2, "%*c");

                 if(ianom){
                    anomaly(sradf1, coslat, irnr, irnth, irad);
                    anomaly(sradf2, coslat, irnr, irnth, irad);
                 }

                 kecalc(sradf1, sradf2, ske2, irdim);

                 for(n=0; n < irdim; n++) {
                    ke1 = *(ske1 + n);
                    ke2 = *(ske2 + n);
                    if(ke1 < ADD_CHECK && ke2 < ADD_CHECK)
                       *(sdiff + n) = 0.5 * (ke2 - ke1);
                    else
                       *(sdiff + n) = ADD_UNDEF;
                 }

              }

              else{
                 fseeko(frad1, ((nwtr + iptc + 1) * irnf + k) * blklen + place1, SEEK_SET);
                 fseeko(frad2, ((nwtr + iptc + 1) * irnf + k) * blklen + place1, SEEK_SET);
                 fread(sradf1, irdim*sizeof(float), 1, frad1);
                 fscanf(frad1, "%*c");
                 fread(sradf2, irdim*sizeof(float), 1, frad2);
                 fscanf(frad2, "%*c");

                 if(ianom){
                    anomaly(sradf1, coslat, irnr, irnth, irad);
                    anomaly(sradf2, coslat, irnr, irnth, irad);
                 }

                 kecalc(sradf1, sradf2, ske2, irdim);

                 for(n=0; n < irdim; n++) {
                    ke1 = *(ske1 + n);
                    ke2 = *(ske2 + n);

                    if(ke1 < ADD_CHECK && ke2 < ADD_CHECK)
                       *(sdiff + n) = ke2 - ke1;
                    else
                       *(sdiff + n) = ADD_UNDEF;
		       
                 }
              }

              fwrite(sdiff, irdim * sizeof(float), 1, frado);
              fprintf(frado, "\n");

           }

           ++iptc;
           ++iptnum;
       }

       if(dtype){

          for(j=0; j < irnf; j++){
              fwrite(sdum, irdim * sizeof(float), 1, frado);
              fprintf(frado, "\n");
          }

          ++iptnum;

       }

       nwtr += atr->num;

    }

    printf("%ld\n", iptnum);

    fseeko(frado, 0L, SEEK_SET);
    fprintf(frado, "%6d %10ld %3d %3d %2d %2d %3d %2d %2d %2d\n", irtrn, iptnum, irnth, irnr, irnf, idir, ndsmth, ifcnt, itpadd, iprojd);

    fclose(frado);

    free(sdum);
    free(sdiff);
    free(sradf1);
    free(sradf2);
    free(ske1); free(ske2);
    free(slat1); free(slat2);
    free(slng1); free(slng2);
    free(coslat);

    return 0;
}

void kecalc(float *u, float *v, float *ke, int irdim)
{
    int i;
    float uu, vv;

    for(i=0; i < irdim; i++) {
        uu = *(u + i);
        vv = *(v + i);
        if(uu < ADD_CHECK && vv < ADD_CHECK)
           *(ke + i) = 0.5 * (uu * uu + vv * vv);
        else
           *(ke + i) = ADD_UNDEF;
    }

}

/* open radial data file for read */

FILE *open_radial_file(char *filrad, int *irtrn, long int *iptnum, int *irnth, int *irnr, int *irnf, int *idir, int *ndsmth, int *ifcnt, int *itpadd, int *iprojd)
{

    int irtrn2=0;
    int irnth2=0, irnr2=0, irnf2=0;
    int idir2=0, ndsmth2=0;
    int ifcnt2=0, itpadd2=0, iprojd2=0;

    long int iptnum2=0;

    static int ifrst=0;

    FILE *frad=NULL;
    char line[MAXCHR];

    frad = fopen(filrad, "r");
    if(!frad){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filrad);
       exit(1);
    }
    fgets(line, MAXCHR, frad);

    if(!ifrst){
       sscanf(line, "%d %ld %d %d %d %d %d %d %d %d", irtrn, iptnum, irnth, irnr, irnf, idir, ndsmth, ifcnt, itpadd, iprojd);
       ifrst = 1;
    }
    else{
       sscanf(line, "%d %ld %d %d %d %d %d %d %d %d", &irtrn2, &iptnum2, &irnth2, &irnr2, &irnf2, &idir2, &ndsmth2, &ifcnt2, &itpadd2, &iprojd2);

       if(*irtrn != irtrn2 || *iptnum != iptnum2 || *irnth != irnth2 ||
          *irnr != irnr2   || *irnf != irnf2 || *idir != idir2  || *ndsmth != ndsmth2 ||
          *ifcnt != ifcnt2 || *itpadd != itpadd2 || *iprojd != iprojd2                   ){
           printf("****ERROR****, incompatability between radial field files.\n\n");
           exit(1);
       }
    }

    return frad;
}

/* determine the deviation from the areal mean */

double anomaly(float *sradf, double *coslat, int irnr, int irnth, int iran)
{
     int i, j;
     int irdim=0;

     double sum=0.0, sumcl=0.0; 

     irdim = irnr * irnth;

     for(i=0; i < iran; i++) {
        for(j=0; j < irnth; j++){
           if(*(sradf + i * irnth + j) > ADD_CHECK) continue;
           sum += *(coslat + i) * *(sradf + i * irnth + j);
           sumcl += *(coslat + i);
        }
     }

     if(sumcl > 0.0)
        sum /= sumcl;
     else
        sum = ADD_UNDEF;

     for(i=0; i < irdim; i++) {
        if(*(sradf + i) > ADD_CHECK) continue;
        *(sradf + i) -= sum;
     }

    return sum;
}
