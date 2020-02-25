#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "st_fo.h"
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"

#define DISTLARGE  180.0
#define TIMINT     6
#define HUGENUM    100000

/* Program to compare an EPS ensemble */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
int toverlap(struct tot_tr * , struct tot_tr * , long int * , long int * );
float trdist_eps(struct tot_tr * , struct tot_tr * , float , long int , long int , int * , int * , int * , int , int , int , int );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
long int new_time(long int , int );
float measure(struct feature_pts *  , struct feature_pts * );
void convert_track(struct tot_tr * , int , int , int );

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int main(int argc, char **argv)
{

   int i, j, k, l;
   int tnum_c, tnum_e;
   int gpr1, ipr1, gpr2, ipr2;
   int iper_num_t;
   int ty_match=0;
   int nens=0;
   int isep=0;
   int jtrmin=0;
   int numpm=0, nump=0;
   int it1=0, it2=0;
   int in1=0, in2=0;
   int ntrpt=100000000;
   int nsamp=0;
   int nhisd=0, nhiss=0;
   int nmat=0;
   int ifr=0;
   int sumhr=0;
   int idty=0;
   int numtr=0;
   int ngg=0, ngtr=0;
   int inorm=0;
   int asty=0;
   int ntims=0;
   int ireg=0, nreg=0;
   int nr1=0, nr2=0, nrs=0;
   int ntpta=HUGENUM;
   int nmaxpt=0, nmaxal=0;
   int ipt=0;
   int isepty=0;
   int itcheck=1;
   
   int timint=TIMINT;

   int *tsamp=NULL;
   int *nspd=NULL, *nsps=NULL;
   int ths=0, the=0;
   int insamp=0;
   int ndt=0, nst=0;
   int *trmatn=NULL;
   int nmth=HUGENUM, ifw=0;
   int iskip=0, niskp=0;

   int nfft, nfldt=0;
   int *nfwpost=NULL;

   int **napt=NULL;

   long int is1=0, is2=0, dt;
   long int stdate=0, stdate1=0, ctim=0, newt=0;
   long int tmin=0;

   float alat1=0.0, alng1=0.0, alat2=0.0, alng2=0.0;
   float sum_per_t=0.0;
   float disth=0.0;
   float timth=0.0;
   float distm=0.0, dist=0.0;
   float sdif=0.0;
   float frat=0.0, ffrat=0.0;
   float hisds, hisde, hisss, hisse;
   float dtd=0, dts=0;
   float rlng1, rlng2, rlat1, rlat2;
   float lngtmp=0.0;

   double **ddist=NULL, **dstr=NULL;

   float *dmean=NULL, *smean=NULL, *absmean=NULL;
   float *dstd=NULL, *sstd=NULL;
   float *bind=NULL, *bins=NULL;

   char cntlfil[MAXCHR], ensmfil[MAXCHR];
   char trout[MAXCHR];
   char timchr[MAXCHR];

   FILE *fin=NULL, *fout=NULL;
   FILE *fcmp=NULL;

   struct tot_tr *trcntl=NULL, *trens=NULL;
   struct tot_tr *atcntl=NULL, *atens=NULL;
   struct tot_tr **trmat=NULL, *trtmp=NULL;
   struct tot_tr *mtrack=NULL;
   struct fet_pt_tr *fp1=NULL, *fp2=NULL, *fpp=NULL;
   struct feature_pts fpts1, fpts2;

   fcmp = fopen("region.dat", "r");
   if(fcmp){
      printf("****WARNING****, file exists with region data for sub-setting tracks,   \r\n"
             "                 this information will be used to restrict the analysis \r\n"
             "                 to specified region.                                   \n\n");
      fscanf(fcmp, "%f %f", &rlng1, &rlng2);
      fscanf(fcmp, "%f %f", &rlat1, &rlat2);
      if(rlng1 > rlng2){
         printf("****ERROR****, regions can't currently stradle the Greenwich meridion.\n\n");
         exit(1);
      }
      ireg = 1;
      fclose(fcmp);
   }

/* read in histogram bin information */

   fin = fopen("hist.dat", "r");
   if(!fin){
      printf("****ERROR****, histogram info file, hist.dat, cannot be read.\n\n");
      exit(1);
   }
   ifr = fscanf(fin, "%d %d %d", &nsamp, &ths, &the);
   if(ifr != 3) {printf("****ERROR****, insufficient input in hist.dat\n\n"); exit(1);}
   ifr = fscanf(fin, "%d %f %f", &nhisd, &hisds, &hisde);
   if(ifr != 3) {printf("****ERROR****, insufficient input in hist.dat\n\n"); exit(1);}
   ifr = fscanf(fin, "%d %f %f", &nhiss, &hisss, &hisse);
   if(ifr != 3) {printf("****ERROR****, insufficient input in hist.dat\n\n"); exit(1);}
   fclose(fin);

   printf("Perform a verification against analysis, '0', or a Lorenz style analysis '1'\n\n");
   scanf("%d", &asty);

   if(!asty){
      printf("What is the control (analysis) track file to compare?\n\n");         
      scanf("%s", cntlfil);
   }

   else {
      printf("What is the first ensemble track file and starting data?\n\n");
      scanf("%s %ld", cntlfil, &stdate);
      printf("What is the target time seperation of track ensembles, in number of time steps?\n\n");
      scanf("%d", &ntims);
      printf("Check number of time periods are consistent, '0' for no and '1' for yes.\n\n");
      scanf("%d", &itcheck);
      if(itcheck < 0 || itcheck > 1) itcheck = 1;
   }
   
   printf("What is the timestep of the data?\n\n");
   scanf("%d", &timint);
   if(timint < 0) {
      printf("****ERROR****, timestep must be positive.\n\n");
      exit(1);
   }

/* read control track ensemble */

   printf("Do you want actual geodesic seperation distance or the orthogonal seperation distance, \r\n"
          "Input '0' or '1'                                                                       \n\n");
   scanf("%d", &isepty);

   if(isepty < 0 || isepty > 1){
     printf("****ERROR****, incorrect input for seperation distance measure.\n\n");
     exit(1);
   }

   fin = fopen(cntlfil, "r");
   if(!fin){
      printf("Can't open file %s for read.\n\n", cntlfil);
      exit(1);
   }

   trcntl = read_tracks(fin, &tnum_c, &gpr1, &ipr1, 's', &alat1, &alng1, NULL, NULL);

   fclose(fin);

   nfft = nff;
   nfldt = nfld;
   nfwpost = (int *)calloc(nff, sizeof(int));
   mem_er((nfwpost == NULL) ? 0 : 1, nff * sizeof(int));
   for(i=0; i < nff; i++) *(nfwpost + i) = *(nfwpos + i);

   if(isepty) convert_track(trcntl, tnum_c, 0, 0);

   if(ireg){
      nr1 = 0;
      for(i=0; i < tnum_c; i++){
          atcntl = trcntl + i;
          nreg = 0;
          for(j=0; j < atcntl->num; j++){
             fp1 = atcntl->trpt + j;
             if((rlng2 - fp1->xf) * (fp1->xf - rlng1) >= 0.0 &&
                (rlat2 - fp1->yf) * (fp1->yf - rlat1) >= 0.0    ) ++nreg;
          } 
          if((float)nreg / (float)(atcntl->num) < 0.6) {
             if(!asty) {
                if(nff){
                  for(j=0; j < atcntl->num; j++) {
                     fpp = atcntl->trpt + j;
                     free(fpp->add_fld);
                  }
                }
                atcntl->num = 0;
                free(atcntl->trpt);
             }
             continue;
          }
          ++nr1; 

      }

   }

/* write to file all analysis tracks that ar in region */

   if(!asty){
      fout = fopen("tr_analysis_reg", "w");
      meantrd(fout, trcntl, tnum_c, 's', gpr1, ipr1, alat1, alng1);
      fclose(fout);
   }

   printf("Number of tracks in first file is %d\n", tnum_c);
   printf("Number of tracks in region in first file is %d\n", nr1);

   trmat = (struct tot_tr **)calloc(tnum_c, sizeof(struct tot_tr *));
   mem_er((trmat == NULL) ? 0 : 1, tnum_c * sizeof(struct tot_tr *));
   trmatn = (int *)calloc(tnum_c, sizeof(int));
   mem_er((trmatn == NULL) ? 0 : 1, tnum_c * sizeof(int));
   for(i=0; i < tnum_c; i++){
       trmatn[i] = 1;
       trmat[i] = (struct tot_tr *)malloc_initl(sizeof(struct tot_tr));
       mem_er((trmat[i] == NULL) ? 0 : 1, sizeof(struct tot_tr));
       trtmp = trmat[i];
       trtmp->num = (trcntl + i)->num;
       trtmp->time = (trcntl + i)->time;
       trtmp->trid = 0;
       trtmp->awt = 0;
       trtmp->trpt = (struct fet_pt_tr *)calloc(trtmp->num, sizeof(struct fet_pt_tr));
       mem_er((trtmp->trpt == NULL) ? 0 : 1, trtmp->num * sizeof(struct fet_pt_tr));
       memcpy(trtmp->trpt, (trcntl + i)->trpt, trtmp->num * sizeof(struct fet_pt_tr));

       if(nff){
         for(j=0; j < trtmp->num; j++){
             fp1 = trtmp->trpt + j;
             fp2 = (trcntl + i)->trpt + j;
             fp1->add_fld = (float *)calloc(nfld, sizeof(float));
             mem_er((fp1->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
             for(k=0; k < nfld; k++) fp1->add_fld[k] = fp2->add_fld[k];
         }
      }
   }

   if(!(trcntl->time)){
      printf("***ERROR***, incorrect temporal identity\n\n");
      exit(1);
   }

   iper_num_t = iper_num;
   sum_per_t = sum_per;

   printf("How many ensemble members are there?\n\n");
   scanf("%d", &nens);

   if(nens < 1) {
      printf("****ERROR****, Not enough ensemble members.\n\n");
      exit(1);
   }


   printf("When matching tracks, do you want,           \r\n"
          "Whole track matching,            input '0'   \r\n"
          "First portion of track matching, input '1'   \r\n"
          "Genesis matching,                input '2'   \n\n");
   scanf("%d", &ty_match);

   if(ty_match < 0 || ty_match > 2){
      printf("****ERROR****, incorrect matching identifier.\n\n");
      exit(1);
   }

   if(ty_match == 1){
      printf("How many track time points should matching be done over from genesis of the ensemble member track?\n\n");
      scanf("%d", &ntrpt);
   }
   else if (ty_match == 2) ntrpt = 1;

   if(ntrpt < 1) {
      printf("****ERROR****, number of points to match over must be >= 1.\n\n");
      exit(1);
   }

   printf("Specify the number of time points from start of forecast within which genesis must occur.\n\n");
   scanf("%d", &ngg);

   if(ngg < 1){
      printf("****ERROR****, number of time points must be positive for genesis threshold.\n\n");
      exit(1);
   }

   printf("What is the seperation distance and temporal overlap thresholds?\n\n");
   scanf("%f %f", &disth, &timth);

   printf("Do you want mean or min seperation matching, '0' or '1'\n\n");
   scanf("%d", &isep);

   if(nsamp < 1 || nhisd < 1 || nhiss < 1) {
      printf("****ERROR****, number of time sampling points must be >= 1\n\n");
      exit(1);
   }

   printf("Do you want to write matching and no matching tracks to file, 'y' or 'n'\n\n");
   scanf("\n");
   if(getchar() == 'y'){
      ifw = 1;
      printf("What is the minimum number of matches for which tracks should be \r\n"
             "written to file.                                                 \n\n");
      scanf("%d", &nmth);

   }

   dt = (the - ths) / nsamp;

   printf("Time sampling is every %ld hours.\n\n", dt);

   ++nsamp;

   tsamp = (int *)calloc(nsamp, sizeof(int));
   mem_er((tsamp == NULL) ? 0 : 1, nsamp * sizeof(int));

   for(i=0; i < nsamp; i++) *(tsamp + i) = ths + i * dt;

   nspd = (int *)calloc(nsamp, sizeof(int));
   mem_er((nspd == NULL) ? 0 : 1, nsamp * sizeof(int));

   nsps = (int *)calloc(nsamp, sizeof(int));
   mem_er((nsps == NULL) ? 0 : 1, nsamp * sizeof(int));

   dmean = (float *)calloc(nsamp, sizeof(float));
   mem_er((dmean == NULL) ? 0 : 1, nsamp * sizeof(int));

   smean = (float *)calloc(nsamp, sizeof(float));
   mem_er((smean == NULL) ? 0 : 1, nsamp * sizeof(int));

   absmean = (float *)calloc(nsamp, sizeof(float));
   mem_er((absmean == NULL) ? 0 : 1, nsamp * sizeof(int));

   dstd = (float *)calloc(nsamp, sizeof(float));
   mem_er((dstd == NULL) ? 0 : 1, nsamp * sizeof(int));

   sstd = (float *)calloc(nsamp, sizeof(float));
   mem_er((sstd == NULL) ? 0 : 1, nsamp * sizeof(int));

   ddist = (double **)calloc(nsamp, sizeof(double *));
   mem_er((ddist == NULL) ? 0 : 1, nsamp * sizeof(double *));

   dstr = (double **)calloc(nsamp, sizeof(double *));
   mem_er((dstr == NULL) ? 0 : 1, nsamp * sizeof(double *));

   bind = (float *)calloc(nhisd, sizeof(float));
   mem_er((bind == NULL) ? 0 : 1, nhisd * sizeof(float));

   bins = (float *)calloc(nhiss, sizeof(float));
   mem_er((bins == NULL) ? 0 : 1, nhiss * sizeof(float));

/* calculate bins */

   dtd = (hisde - hisds) / nhisd;
   printf("Seperation distance bin widths are %f\n\n", dtd);
   for(i=0; i < nhisd; i++) *(bind + i) = hisds + 0.5 * dtd + i * dtd;
   dts = (hisse - hisss) / nhiss;
   printf("Intensity difference bin widths are %f\n\n", dts);
   for(i=0; i < nhiss; i++) *(bins + i) = hisss + 0.5 * dts + i * dts;

   for(i=0; i < nsamp; i++){
       *(ddist + i) = (double *)calloc(nhisd, sizeof(double));
       mem_er((*(ddist + i) == NULL) ? 0 : 1, nhisd * sizeof(double));
       *(dstr + i) = (double *)calloc(nhiss, sizeof(double));
       mem_er((*(dstr + i) == NULL) ? 0 : 1, nhiss * sizeof(double));
   }
   
/*   iskip = 0; */

   for(i=0; i < nens; i++){

       printf("What is the next ensemble member track file and start date?\r\n"
              "Date should be in the format YYYYMMDDHH.                   \n\n");
       scanf("%s %ld", ensmfil, &stdate1);

       iskip = 0; 

       if(!asty) stdate = stdate1;
       else {
          newt = stdate;
          for(j=0; j <ntims; j++) newt = new_time(newt, timint);
          sprintf(timchr, "%ld", newt);
          if(newt != stdate1 || !strstr(ensmfil, timchr)){
             printf("****WARNING****, dates %ld and %ld differ by more than the requested number of time steps %d,\r\n"
                    "                    or %ld not found in filename, do you want to continue, 'y' or 'n'.       \n\n", 
                                         stdate, stdate1, ntims, newt);
             if(!niskp){
                scanf("\n");
                if(getchar() == 'n') exit(1);
                else iskip = 1;
                niskp = 1;
             }
             else iskip = 1;
          }
       }

       fin = fopen(ensmfil, "r");
       if(!fin){
          printf("Can't open file %s for read.\n\n", ensmfil);
          exit(1);
       }

       trens = read_tracks(fin, &tnum_e, &gpr2, &ipr2, 's', &alat2, &alng2, NULL, NULL);

       fclose(fin);

       if(isepty) convert_track(trens, tnum_e, 0, 0);

       if((gpr1 != gpr2 || ipr1 != ipr2) ||
          (fabs(alat1 - alat2) > 1.0e-4 || fabs(alng1 - alng2) > 1.0e-4)){
          printf("***ERROR***, projection parameters do not match\n\n");
          exit(1);
       }

       if(!(trens->time)){
          printf("***ERROR***, Incorrect temporal identity.\n\n");
          exit(1);
       }

       if(itcheck){
          if(iper_num != iper_num_t || fabs(sum_per - sum_per_t) > 1.0e-6){
             printf("****ERROR****, number of time periods do not match\n\n");
             exit(1);
          }
       }

       if(nff != nfft || nfld != nfldt){
          printf("***ERROR***, additional field parameters do not match for file %s.\n\n", ensmfil);
          exit(1);
       }
       for(j=0; j < nff; j++){
           if(*(nfwpost + j) != *(nfwpos + j)){
             printf("***ERROR***, additional field location parameters do not match for additional field %d\n\n", j+1);
             exit(1);
           }
       }

/*       if(iskip && asty) {stdate = stdate1; continue;} */
       if(iskip && asty) goto SKIP;

       for(j=0; j < tnum_e; j++){

          atens = trens + j;

          ngtr = 0;
          ctim = stdate;

          while(atens->trpt->time != ctim) {++ngtr; ctim = new_time(ctim, timint);}

          if(ireg){
             nreg = 0;
             for(k=0; k < atens->num; k++){
                 fp1 = atens->trpt + k;
                 if((rlng2 - fp1->xf) * (fp1->xf - rlng1) >= 0.0 &&
                    (rlat2 - fp1->yf) * (fp1->yf - rlat1) >= 0.0    ) ++nreg;
             }

             if((float)nreg / (float)(atens->num) < 0.6) {
                if(!asty) {
                   if(nff){
                     for(k=0; k < atens->num; k++) {
                        fpp = atens->trpt + k;
                        free(fpp->add_fld);
                     }
                   }
                   atens->num = 0;
                   free(atens->trpt);
                }
                continue;
             }
             ++nr2;
             if(ngtr <= ngg ) ++nrs;
          }

          else {
            if(ngtr <= ngg ) ++nrs;
          }

          jtrmin = -1;
          distm = DISTLARGE;
          numpm = 0;
          frat = 0.0;

          if(atens->num){

             nmat = 0;

             for(k=0; k < tnum_c; k++){

                 atcntl = trcntl + k;

                 if(ireg){
                    nreg = 0;

                    for(l=0; l < atcntl->num; l++){
                        fp2 = atcntl->trpt + l;
                        if((rlng2 - fp2->xf) * (fp2->xf - rlng1) >= 0.0 &&
                           (rlat2 - fp2->yf) * (fp2->yf - rlat1) >= 0.0    ) ++nreg;
                    }
                    if((float)nreg / (float)(atcntl->num) < 0.6) continue;
                 }

                 if(atcntl->num){

/* check for overlap in time */

                    if(toverlap(atens, atcntl, &is1, &is2)){

                       dist = trdist_eps(atens, atcntl, disth, is1, is2, &nump, &it1, &it2, isep, ty_match, ntrpt, isepty);
                       ffrat = 2.0 * nump / (float)(atens->num + atcntl->num);

/* find number of time points for genesis from start of forecast */

                       if(dist <= disth && ffrat >= timth && ngtr <= ngg){

                          ++nmat;

                          if(jtrmin < 0 ){
                             jtrmin = k;
                             distm = dist;
                             frat = ffrat;
                             numpm = nump;
                             in1 = it1;
                             in2 = it2;
                          }

                          else if(jtrmin >= 0 && dist < distm) {
                             jtrmin = k;
                             distm = dist;
                             frat = ffrat;
                             numpm = nump;
                             in1 = it1;
                             in2 = it2;
                          }

                       }

                    }

                 }

             }

          }

          if(jtrmin >= 0){

             ++numtr;

/* add track to relevent control ensemble for output */

             ++(trmatn[jtrmin]);
             trmat[jtrmin] = (struct tot_tr *)realloc_n(trmat[jtrmin], trmatn[jtrmin] * sizeof(struct tot_tr));
             mem_er((trmat[jtrmin] == NULL) ? 0 : 1,  trmatn[jtrmin] * sizeof(struct tot_tr));
             trtmp = trmat[jtrmin] + trmatn[jtrmin] - 1;
             trtmp->num = atens->num;
             trtmp->time = atens->time;
             trtmp->trid = i + 1;
             trtmp->awt = 0;
             trtmp->trpt = (struct fet_pt_tr *)calloc(trtmp->num, sizeof(struct fet_pt_tr));
             mem_er((trtmp->trpt == NULL) ? 0 : 1, trtmp->num * sizeof(struct fet_pt_tr));
             memcpy(trtmp->trpt, atens->trpt, trtmp->num * sizeof(struct fet_pt_tr));

             if(nff){
                for(k=0; k < trtmp->num; k++){ 
                    fp1 = trtmp->trpt + k;
                    fp2 = atens->trpt + k;
                    fp1->add_fld = (float *)calloc(nfld, sizeof(float));
                    mem_er((fp1->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                    for(l=0; l < nfld; l++) fp1->add_fld[l] = fp2->add_fld[l];
                }
             }

             atcntl = trcntl + jtrmin;
             ctim = stdate;
             sumhr = 0;
             for(k=0; k < numpm; k++){
                 fp1 = atens->trpt + in1 + k;
                 fp2 = atcntl->trpt + in2 + k;
                 if(fp1->time != fp2->time){
                    printf("****WARNING****, times for matched tracks do not correspond.\n\n");
                 }

                 while(ctim != fp1->time){
                      sumhr += timint;
                      ctim = new_time(ctim, timint);
                 }

                 if(sumhr / dt >= nsamp) break;

/* for differences, ensemble member - control */

                 if(!(sumhr % dt)){
                    insamp = sumhr / dt;
                    (fpts1.x).xy = fp1->xf;
                    (fpts1.y).xy = fp1->yf;

                    (fpts2.x).xy = fp2->xf;
                    (fpts2.y).xy = fp2->yf;
                    dist = measure(&fpts1, &fpts2) / FP_PI;

                    sdif = fp1->zf - fp2->zf;
                    if(dist < 0.0 || dist > hisde){
                       printf("****WARNING****, distance difference value %f out of range.\n\n", dist);
                    }
                    else {
                       ++(*(nspd + insamp));
                       *(ddist[insamp] + (int)((dist - hisds) / dtd)) += 1.0;
                       *(dmean + insamp) += dist;
                       *(dstd + insamp) += dist * dist;
                    }
                    if(sdif < hisss || sdif > hisse){
                       printf("****WARNING****, intensity difference value %f out of range.\n\n", sdif);
                    }
                    else {
                       ++(*(nsps + insamp));
                       *(dstr[insamp] + (int)((sdif - hisss) / dts)) += 1.0;
                       *(absmean + insamp) += fabs(sdif);
                       *(smean + insamp) += sdif;
                       *(sstd + insamp) += sdif * sdif; 
                    }   

                 }

             }

             if(!asty){
                if(nff){
                  for(k=0; k < atens->num; k++) {
                     fpp = atens->trpt + k;
                     free(fpp->add_fld);
                  }
                }
                atens->num = 0;
                free(atens->trpt);
             }


          }

       }

SKIP:

       if(!asty){
          if(ifw){
             sprintf(trout, "trnomatch_ens%04d", i+1);
             fout=fopen(trout, "w");
             if(!fout){
                printf("****ERROR****, can't open file %s for write.\n\n", trout);
                exit(1);
             }
             meantrd(fout, trens, tnum_e, 's', gpr1, ipr1, alat1, alng1);
             fclose(fout);
          }
          for(j=0; j < tnum_e; j++) {
             if(nff){
               for(k=0; k < (trens + j)->num; k++) {
                  fpp = (trens + j)->trpt + k;
                  free(fpp->add_fld);
               }
             }
             free((trens + j)->trpt);
          }
          free(trens);
       }
       else {

          for(j=0; j < tnum_c; j++) {
             if(ifw && trmatn[j] >= nmth){
                sprintf(trout, "trmatch_ens_%10ld_tr%04d", stdate, j+1);
                fout=fopen(trout, "w");
                if(!fout){
                   printf("****ERROR****, can't open file %s for write.\n\n", trout);
                   exit(1);
                }
                meantrd(fout, *(trmat + j), trmatn[j], 's', gpr1, ipr1, alat1, alng1);
                fclose(fout);
             }
             free((trcntl + j)->trpt);
             for(k=0; k < trmatn[j]; k++){
                 trtmp = *(trmat + j) + k;
                 if(nff){
                    for(l=0; l < trtmp->num; l++) {
                       fpp = trtmp->trpt + l;
                      free(fpp->add_fld);
                    }
                 }
                 free(trtmp->trpt);
             }
             free(*(trmat + j));

          } 
          free(trcntl);
          trcntl = trens;
          tnum_c = tnum_e;
          stdate = stdate1;
          trmat = (struct tot_tr **)realloc_n(trmat, tnum_c * sizeof(struct tot_tr *));
          mem_er((trmat == NULL) ? 0 : 1, tnum_c * sizeof(struct tot_tr *));
          trmatn = (int *)realloc_n(trmatn, tnum_c * sizeof(int));
          mem_er((trmatn == NULL) ? 0 : 1, tnum_c * sizeof(int));
          for(j=0; j < tnum_c; j++){
             trmatn[j] = 1;
             trmat[j] = (struct tot_tr *)malloc_initl(sizeof(struct tot_tr));
             mem_er((trmat[j] == NULL) ? 0 : 1, sizeof(struct tot_tr));
             trtmp = trmat[j];
             trtmp->num = (trcntl + j)->num;
             trtmp->time = (trcntl + j)->time;
             trtmp->trid = i + 1;
             trtmp->awt = 0;
             trtmp->trpt = (struct fet_pt_tr *)calloc(trtmp->num, sizeof(struct fet_pt_tr));
             mem_er((trtmp->trpt == NULL) ? 0 : 1, trtmp->num * sizeof(struct fet_pt_tr));
             memcpy(trtmp->trpt, (trcntl + j)->trpt, trtmp->num * sizeof(struct fet_pt_tr));

             if(nff){
                for(k=0; k < trtmp->num; k++){
                    fp1 = trtmp->trpt + k;
                    fp2 = (trcntl + j)->trpt + k;
                    fp1->add_fld = (float *)calloc(nfld, sizeof(float));
                    mem_er((fp1->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                    for(l=0; l < nfld; l++) fp1->add_fld[l] = fp2->add_fld[l];
                }
             }
          }
       }

   }

   if(!asty){

      printf("Do you want to produce a mean track? \n\n");
      scanf("\n");
      if(getchar() == 'y'){
         printf("What is the minimum number of track points to average over?\n\n");
         scanf("%d", &ntpta);
      }

      for(i=0; i < tnum_c; i++){
        if(ifw && trmatn[i] >= nmth){
           sprintf(trout, "trmatch_cntl_tr%04d", i+1);
           fout=fopen(trout, "w");
           if(!fout){
              printf("****ERROR****, can't open file %s for write.\n\n", trout);
              exit(1);
           }
           meantrd(fout, *(trmat + i), trmatn[i], 's', gpr1, ipr1, alat1, alng1);
           fclose(fout);

           if(trmatn[i]-1 >= ntpta){

              strcat(trout, "_mean");
              nmaxpt = 0;
              tmin=(trmat[i] + 1)->trpt->time;
              for(j=1; j < trmatn[i]; j++) {
                  if((trmat[i] + j)->num > nmaxpt) nmaxpt = (trmat[i] + j)->num;
                  if((trmat[i] + j)->trpt->time < tmin) tmin = (trmat[i] + j)->trpt->time;
              }
         
              nmaxal = 2 * nmaxpt + (nmaxpt/2);
              mtrack = (struct tot_tr *)malloc_initl(sizeof(struct tot_tr));
              mem_er((mtrack == NULL) ? 0 : 1, sizeof(struct tot_tr));
              mtrack->trpt = (struct fet_pt_tr *)calloc(nmaxal, sizeof(struct fet_pt_tr));
              mem_er((mtrack->trpt == NULL) ? 0 : 1, nmaxal * sizeof(struct fet_pt_tr));
              napt = (int **)calloc(nmaxal, sizeof(int *));
              mem_er((napt == NULL) ? 0 : 1, nmaxal * sizeof(int *));
              if(nff){
                 napt = (int **)calloc(nmaxal, sizeof(int *));
                 mem_er((napt == NULL) ? 0 : 1, nmaxal * sizeof(int *));
                 for(j=0; j < nmaxal; j++){
                     fp1 = mtrack->trpt + j;
                     fp1->add_fld = (float *)calloc(nfld, sizeof(float));
                     mem_er((fp1->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                     *(napt + j) = (int *)calloc(nfld, sizeof(int)); 
                     mem_er((*(napt + j) == NULL) ? 0 : 1, nfld * sizeof(int));
                 }
              }
              mtrack->num = 0;
              for(j=1; j < trmatn[i]; j++) {
                  trtmp = trmat[i] + j;
                  ipt=0;
                  ctim = tmin;
                  while(ctim != trtmp->trpt->time) {++ipt; ctim = new_time(ctim, timint);}
                  if(j == 1){
                     memcpy(mtrack->trpt + ipt, trtmp->trpt, trtmp->num * sizeof(struct fet_pt_tr));

                     if(nff){
                        for(k=0; k < trtmp->num; k++){
                            fp1 = mtrack->trpt + ipt + k;
                            fp2 = trtmp->trpt + k;
                            for(l=0; l < nfld; l++) {
                                if(fp2->add_fld[l] > ADD_CHECK) continue; 
                                fp1->add_fld[l] = fp2->add_fld[l];
                                *(*(napt + ipt + k) + l) = 1;
                            }
                        }
                     }

                     for(k=0; k < trtmp->num; k++){
                         fp1 = mtrack->trpt + ipt + k;
                         fp1->npt = 1;
                     }
                  }
                  else {
                     for(k=0; k < trtmp->num; k++){
                         fp1 = mtrack->trpt + ipt + k;
                         fp2 = trtmp->trpt + k;
                         if(fp1->time && fp1->time != fp2->time){
                            printf("****ERROR****, times do not match for track averaging, %ld %ld.\n\n", fp1->time, fp2->time);
                            exit(1);
                         }
                         lngtmp = fp1->xf / fp1->npt;

                         ++(fp1->npt);
                         if(fabs(lngtmp - fp2->xf) > 180.0){
                            if(fabs(360.0 - lngtmp) < 180.0) fp1->xf += (fp2->xf + 360.0);
                            else fp1->xf += (fp2->xf - 360.0);
                         }
                         else fp1->xf += fp2->xf;
                         fp1->yf += fp2->yf;
                         fp1->zf += fp2->zf;
                         if(!(fp1->time)) fp1->time = fp2->time;

                         if(nff){
                            for(l=0; l < nfld; l++){
                                if(fp2->add_fld[l] > ADD_CHECK) continue;
                                fp1->add_fld[l] += fp2->add_fld[l];
                                ++(*(*(napt + ipt + k) + l));
                            }
                         }
                     }
                  }

              }

              mtrack->time = tmin; 

              for(j=0; j < nmaxal; j++){
                  fp1 = mtrack->trpt + j;
                  if(fp1->npt >= ntpta) {
                     fp1->xf /= (float) (fp1->npt);
                     if(fp1->xf > 360.0) fp1->xf -= 360.0;
                     else if(fp1->xf < 0.0) fp1->xf += 360.0;
                     fp1->yf /= (float) (fp1->npt);
                     fp1->zf /= (float) (fp1->npt);

                     if(nff){
                       for(k=0; k < nfld; k++){
                          if(*(*(napt + j) + k) > 0) fp1->add_fld[k] /= (float)(*(*(napt + j) + k));
                          else fp1->add_fld[k] = ADD_UNDEF;
                       }
                     } 

                     ++(mtrack->num);
                  }
                  else {
                     fp1->xf = 0.0;
                     fp1->yf = 0.0;
                     fp1->zf = 0.0;
                     fp1->time = 0;
                  }

              }

              ipt = 0;
              for(j=0; j < nmaxal; j++) { fp1 = mtrack->trpt + j; if(!fp1->time) ++ipt; else break;}

              fpp = mtrack->trpt;
              mtrack->trpt = mtrack->trpt + ipt;

              if(mtrack->num){
                 fout=fopen(trout, "w");
                 if(!fout){
                    printf("****ERROR****, can't open file %s for write.\n\n", trout);
                    exit(1);
                 }
                 meantrd(fout, mtrack, 1, 's', gpr1, ipr1, alat1, alng1);
                 fclose(fout);
              }
              mtrack->trpt = fpp;

             if(nff){
                for(j=0; j < nmaxal; j++){
                    fpp = mtrack->trpt + j;
                    free(fpp->add_fld);
                    free(*(napt + j));
                }
                free(napt);
             }
             free(mtrack->trpt);
             free(mtrack);

           }

        }

     }
   }

   for(i=0; i < tnum_c; i++) {
      if(nff){
        for(j=0; j < (trcntl + i)->num; j++) {
           fpp = (trcntl + i)->trpt + j;
           free(fpp->add_fld);
        }
      }
      free((trcntl + i)->trpt);
   }
   free(trcntl);

   printf("Do you want frequency distributions or pdf's, '0' or '1'\n\n");
   scanf("%d", &idty);

   printf("Normalize distributions as 1D, or 2D, '0' or '1'\n\n");
   scanf("%d", &inorm);

   if(inorm){
      for(i=0; i < nsamp; i++){
         ndt += nspd[i];
         nst += nsps[i];
      }
   }

   for(i=0; i < nsamp; i++){
       if(!idty){
          if(!inorm){
             for(j=0; j < nhisd; j++)
                *(ddist[i] + j) /= (float) nspd[i];

             for(j=0; j < nhiss; j++)
                *(dstr[i] + j) /= (float) nsps[i];
          }
          else {
             for(j=0; j < nhisd; j++)
                *(ddist[i] + j) /= (float) ndt;

             for(j=0; j < nhiss; j++)
                *(dstr[i] + j) /= (float) nst;
          }
       }
       else {
          if(!inorm){
             for(j=0; j < nhisd; j++) 
                *(ddist[i] + j) /= (dtd * (float) nspd[i]);

             for(j=0; j < nhiss; j++)
                *(dstr[i] + j) /= (dts * (float) nsps[i]);
          }
          else {
             for(j=0; j < nhisd; j++) 
                *(ddist[i] + j) /= (dt * dtd * (float) ndt);

             for(j=0; j < nhiss; j++)
                *(dstr[i] + j) /= (dt * dts * (float) nst);
          }
       }

       if(nspd[i] > 0){
          dmean[i] /= (float) nspd[i];
          dstd[i] = sqrt((dstd[i] / (float) nspd[i]) - dmean[i] * dmean[i]);
       }
       else {dstd[i] = 0.0; dmean[i] = 0.0;}

       if(nsps[i] > 0){
          absmean[i] /= (float) nsps[i];
          smean[i] /= (float) nsps[i];
          sstd[i] = sqrt((sstd[i] / (float) nsps[i]) - smean[i] * smean[i]);
       }
       else {sstd[i] = 0.0; smean[i] = 0.0;}


   }

/* write output to file */

   fout = fopen("ensemble.dat", "w");
   if(!fout){
      printf("****ERROR****, can't open file ensemble.dat for write.\n\n");
      exit(1);
   }

   fprintf(fout, "                Seperation distance distribution\n\n");
   fprintf(fout, "         ");
   for(i=0; i < nsamp; i++){
       fprintf(fout, "%8d  ", *(tsamp + i));
   }
   fprintf(fout, "\n");
   for(i=0; i < nhisd; i++){
       fprintf(fout, "%9.4f ", *(bind + i));
       if(!inorm){
          for(j=0; j < nsamp; j++)
             fprintf(fout, "%9.4f ", *(ddist[j] + i));
       }
       else {
          for(j=0; j < nsamp; j++)
             fprintf(fout, "%e ", *(ddist[j] + i));
       }
       fprintf(fout, "\n");
   }
   fprintf(fout, "\n\n");

   fprintf(fout, "                Intensity difference distribution\n\n");
   fprintf(fout, "         ");
   for(i=0; i < nsamp; i++){
       fprintf(fout, "%8d  ", *(tsamp + i));
   }
   fprintf(fout, "\n");
   for(i=0; i < nhiss; i++){
       fprintf(fout, "%9.4f ", *(bins + i));
       if(!inorm){
          for(j=0; j < nsamp; j++)
             fprintf(fout, "%9.4f ", *(dstr[j] + i));
       }
       else {
          for(j=0; j < nsamp; j++)
             fprintf(fout, "%e ", *(dstr[j] + i));
       }
       fprintf(fout, "\n");
   }
   fprintf(fout, "\n\n");

   fprintf(fout, "Distribution of Means and Std's\n\n");
   fprintf(fout, "             Distance           Intensity\n\n");
   fprintf(fout, "Time       Mean     STD       Mean(Abs)   STD\n\n");
   for(i=0; i < nsamp; i++){
       fprintf(fout, "%4d %9.4f %9.4f %9.4f %9.4f\n\n", tsamp[i], dmean[i], dstd[i], absmean[i], sstd[i]);
   }

   fclose(fout);


   printf("Total number of tracks in region in forecast files is %d\n", nr2);
   printf("Total number of tracks in region and within genesis threshold in forecast file is %d\n", nrs);

   printf("Total number of tracks that match are %d\n", numtr);

   for(i=0; i < tnum_c; i++){
       trtmp = *(trmat + i);
       for(j=0; j < trmatn[i]; j++) free((trtmp +j)->trpt);
       free(trtmp);
   }

   free(trmat);
   free(trmatn);

   free(tsamp);
   for(i=0; i < nsamp; i++){
       free(*(ddist + i));
       free(*(dstr + i));
   }
   free(ddist);
   free(dstr);
   free(bind);
   free(bins);
   free(dmean);
   free(smean);
   free(absmean);
   free(dstd);
   free(sstd);
   free(nspd);
   free(nsps);

   return 0;

}
