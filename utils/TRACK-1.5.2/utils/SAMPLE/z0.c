#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "statistic.h"
#include "mem_er.h"

#define  MAXCHR  500
#define  NSTAT   13

struct tot_stat *read_stats(FILE * );
void statdmp(FILE * , struct tot_stat * );
long double cnorminv(long double );

int prgr, prty;

int main(void)

{

    int i, j;
    int nsamp=0;
    int ptnum=0;
    int ibias=0;

    char statfil[MAXCHR];
    char stub[MAXCHR], infile[MAXCHR];
    char cnum[10];
    char stfilo[]="stats_z0.dat";
    char sflerr[]="stats_err.dat";

    FILE *fin=NULL, *fout=NULL;
    FILE *fstat=NULL;

    float nub=0.0;
    float tlev=0.0, tlim=0.0;
    float zfac=0.0;
    float stat=0.0;

    struct tot_stat *serr=NULL;
    struct tot_stat *smean=NULL;
    struct tot_stat *stato=NULL;
    struct tot_stat *sout=NULL;
    struct tot_stat *sin=NULL;
    struct tot_stat *sbias=NULL;
    struct tot_stat *s_lo=NULL, *s_hi=NULL;
    struct tot_stat *slen=NULL, *sshp=NULL;
    struct pt_stat *stmp1=NULL, *stmp2=NULL, *stmp3=NULL, *stmp4=NULL, *stmp5=NULL;
    struct pt_stat *sbb=NULL;

    
    printf("What is the file with the original statistics?\n\n");
    scanf("%s", statfil);

    fstat = fopen(statfil, "r");
    if(!fstat){
       printf("****ERROR****, unable to open file %s for read.\n\n", statfil);
       exit(1);
    }

    stato = read_stats(fstat);

    fclose(fstat);

    sout = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((sout == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(sout, stato, sizeof(struct tot_stat));

    sout->ptnum = stato->ptnum;
    ptnum = sout->ptnum;

    sout->ptst = (struct pt_stat *)calloc(ptnum, sizeof(struct pt_stat));
    mem_er((sout->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat));

    smean = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((smean == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(smean, stato, sizeof(struct tot_stat));

    smean->ptnum = stato->ptnum;

    smean->ptst = (struct pt_stat *)calloc(ptnum, sizeof(struct pt_stat));
    mem_er((smean->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat));

    serr = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((serr == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(serr, stato, sizeof(struct tot_stat));

    serr->ptnum = stato->ptnum;

    serr->ptst = (struct pt_stat *)calloc(ptnum, sizeof(struct pt_stat));
    mem_er((serr->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat));

    sbias = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((sbias == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(sbias, stato, sizeof(struct tot_stat));

    sbias->ptnum = stato->ptnum;

    sbias->ptst = (struct pt_stat *)calloc(ptnum, sizeof(struct pt_stat));
    mem_er((sbias->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat));


    for(i=0; i < ptnum; i++){
       stmp1 = stato->ptst + i;
       stmp2 = sout->ptst + i;
       stmp3 = smean->ptst + i;
       stmp4 = serr->ptst + i;
       stmp5 = sbias->ptst + i;
       stmp2->xs = stmp3->xs = stmp4->xs = stmp5->xs = stmp1->xs;
       stmp2->ys = stmp3->ys = stmp4->ys = stmp5->ys = stmp1->ys;
    }

    printf("How many samples are there to read?\n\n");
    scanf("%d", &nsamp);

    printf("Input file stub for input files.\n\n");
    scanf("%s", stub);

    nub = (float) nsamp / (float) (nsamp - 1);

    for(i=0; i < nsamp; i++){

       if(i < 10) sprintf(cnum, "00000%d", i);
       else if(i < 100) sprintf(cnum, "0000%d", i);
       else if(i < 1000) sprintf(cnum, "000%d", i);
       else if(i < 10000) sprintf(cnum, "00%d", i);
       strcpy(infile, stub);
       strcat(infile, cnum);

       printf("%s\n", infile);

       fin = fopen(infile, "r");
       if(fin == NULL){
          printf("****ERROR****, cant open file \r\n"
                 "               %s\r\n"
                 "               for read\n\n", infile);
          exit(1);
       }

       sin = read_stats(fin);

       fclose(fin);

       if(sin->ptnum != stato->ptnum){
          printf("****ERROR****, grids in stats files do not match.\n\n");
          exit(1);
       }

       for(j=0; j < ptnum; j++){
           stmp1 = stato->ptst + j;
           stmp2 = sin->ptst + j;
           stmp3 = sout->ptst + j;
           stmp4 = smean->ptst + j;
           stmp5 = serr->ptst + j;
           if((stmp2->stat1).mean < (stmp1->stat1).mean ) (stmp3->stat1).mean += 1.0;
           (stmp4->stat1).mean += (stmp2->stat1).mean;
           (stmp5->stat1).mean += (stmp2->stat1).mean * (stmp2->stat1).mean;
           if((stmp2->stat1).var < (stmp1->stat1).var ) (stmp3->stat1).var += 1.0;
           (stmp4->stat1).var += (stmp2->stat1).var;
           (stmp5->stat1).var += (stmp2->stat1).var * (stmp2->stat1).var;
           if((stmp2->stat2).mean < (stmp1->stat2).mean ) (stmp3->stat2).mean += 1.0;
           (stmp4->stat2).mean += (stmp2->stat2).mean;
           (stmp5->stat2).mean += (stmp2->stat2).mean * (stmp2->stat2).mean;
           if((stmp2->stat2).var < (stmp1->stat2).var ) (stmp3->stat2).var += 1.0;
           (stmp4->stat2).var += (stmp2->stat2).var;
           (stmp5->stat2).var += (stmp2->stat2).var * (stmp2->stat2).var;
           if(stmp2->stat3 < stmp1->stat3 ) stmp3->stat3 += 1.0;
           stmp4->stat3 += stmp2->stat3;
           stmp5->stat3 += stmp2->stat3 * stmp2->stat3;
           if(stmp2->stat4 < stmp1->stat4 ) stmp3->stat4 += 1.0;
           stmp4->stat4 += stmp2->stat4;
           stmp5->stat4 += stmp2->stat4 * stmp2->stat4;
           if(stmp2->stat5 < stmp1->stat5 ) stmp3->stat5 += 1.0;
           stmp4->stat5 += stmp2->stat5;
           stmp5->stat5 += stmp2->stat5 * stmp2->stat5;
           if(stmp2->stat6 < stmp1->stat6 ) stmp3->stat6 += 1.0;
           stmp4->stat6 += stmp2->stat6;
           stmp5->stat6 += stmp2->stat6 * stmp2->stat6;
           if(stmp2->stat8 < stmp1->stat8 ) stmp3->stat8 += 1.0;
           stmp4->stat8 += stmp2->stat8;
           stmp5->stat8 += stmp2->stat8 * stmp2->stat8;
           if(stmp2->stat9 < stmp1->stat9 ) stmp3->stat9 += 1.0;
           stmp4->stat9 += stmp2->stat9;
           stmp5->stat9 += stmp2->stat9 * stmp2->stat9;
           if(stmp2->stat10 < stmp1->stat10 ) stmp3->stat10 += 1.0;
           stmp4->stat10 += stmp2->stat10;
           stmp5->stat10 += stmp2->stat10 * stmp2->stat10;
           if(stmp2->stat12 < stmp1->stat12 ) stmp3->stat12 += 1.0;
           stmp4->stat12 += stmp2->stat12;
           stmp5->stat12 += stmp2->stat12 * stmp2->stat12;
           if(stmp2->stat13 < stmp1->stat13 ) stmp3->stat13 += 1.0;
           stmp4->stat13 += stmp2->stat13;
           stmp5->stat13 += stmp2->stat13 * stmp2->stat13;

       }

       free(sin->ptst);
       free(sin);

    }

    for(i=0; i < ptnum; i++){
        stmp1 = stato->ptst + i;
        stmp3 = sout->ptst + i;
        stmp4 = smean->ptst + i;
        stmp5 = serr->ptst + i;
        sbb = sbias->ptst + i;

        (stmp3->stat1).mean /= (float) nsamp;
        (stmp4->stat1).mean /= (float) nsamp;
        (stmp5->stat1).mean /= (float) nsamp;
        (stmp5->stat1).mean -= (stmp4->stat1).mean * (stmp4->stat1).mean;
        if((stmp5->stat1).mean > 0) (stmp5->stat1).mean = sqrt(nub * (stmp5->stat1).mean);
        else (stmp5->stat1).mean = 0.0;
        (sbb->stat1).mean = (stmp4->stat1).mean - (stmp1->stat1).mean;
        (stmp3->stat1).mean = (float)cnorminv((long double) ((stmp3->stat1).mean));

        (stmp3->stat1).var /= (float) nsamp;
        (stmp4->stat1).var /= (float) nsamp;
        (stmp5->stat1).var /= (float) nsamp;
        (stmp5->stat1).var -= (stmp4->stat1).var * (stmp4->stat1).var;
        if((stmp5->stat1).var > 0) (stmp5->stat1).var = sqrt(nub * (stmp5->stat1).var);
        else (stmp5->stat1).var = 0.0;
        (sbb->stat1).var = (stmp4->stat1).var - (stmp1->stat1).var;
        (stmp3->stat1).var = (float)cnorminv((long double) ((stmp3->stat1).var));

        (stmp3->stat2).mean /= (float) nsamp;
        (stmp4->stat2).mean /= (float) nsamp;
        (stmp5->stat2).mean /= (float) nsamp;
        (stmp5->stat2).mean -= (stmp4->stat2).mean * (stmp4->stat2).mean;
        if((stmp5->stat2).mean > 0) (stmp5->stat2).mean = sqrt(nub * (stmp5->stat2).mean);
        else (stmp5->stat2).mean = 0.0;
        (sbb->stat2).mean = (stmp4->stat2).mean - (stmp1->stat2).mean;
        (stmp3->stat2).mean = (float)cnorminv((long double) ((stmp3->stat2).mean));

        (stmp3->stat2).var /= (float) nsamp;
        (stmp4->stat2).var /= (float) nsamp;
        (stmp5->stat2).var /= (float) nsamp;
        (stmp5->stat2).var -= (stmp4->stat2).var * (stmp4->stat2).var;
        if((stmp5->stat2).var > 0) (stmp5->stat2).var = sqrt(nub * (stmp5->stat2).var);
        else (stmp5->stat2).var = 0.0;
        (sbb->stat2).var = (stmp4->stat2).var - (stmp1->stat2).var;
        (stmp3->stat2).var = (float)cnorminv((long double) ((stmp3->stat2).var));

        stmp3->stat3 /= (float) nsamp;
        stmp4->stat3 /= (float) nsamp;
        stmp5->stat3 /= (float) nsamp;
        stmp5->stat3 -= stmp4->stat3 * stmp4->stat3;
        if(stmp5->stat3 > 0) stmp5->stat3 = sqrt(nub * stmp5->stat3);
        else stmp5->stat3 = 0.0;
        sbb->stat3 = stmp4->stat3 - stmp1->stat3;
        stmp3->stat3 = (float)cnorminv((long double) (stmp3->stat3));

        stmp3->stat4 /= (float) nsamp;
        stmp4->stat4 /= (float) nsamp;
        stmp5->stat4 /= (float) nsamp;
        stmp5->stat4 -= stmp4->stat4 * stmp4->stat4;
        if(stmp5->stat4 > 0) stmp5->stat4 = sqrt(nub * stmp5->stat4);
        else stmp5->stat4 = 0.0;
        sbb->stat4 = stmp4->stat4 - stmp1->stat4;
        stmp3->stat4 = (float)cnorminv((long double) (stmp3->stat4));

        stmp3->stat5 /= (float) nsamp;
        stmp4->stat5 /= (float) nsamp;
        stmp5->stat5 /= (float) nsamp;
        stmp5->stat5 -= stmp4->stat5 * stmp4->stat5;
        if(stmp5->stat5 > 0) stmp5->stat5 = sqrt(nub * stmp5->stat5);
        else stmp5->stat5 = 0.0;
        sbb->stat5 = stmp4->stat5 - stmp1->stat5;
        stmp3->stat5 = (float)cnorminv((long double) (stmp3->stat5));

        stmp3->stat6 /= (float) nsamp;
        stmp4->stat6 /= (float) nsamp;
        stmp5->stat6 /= (float) nsamp;
        stmp5->stat6 -= stmp4->stat6 * stmp4->stat6;
        if(stmp5->stat6 > 0) stmp5->stat6 = sqrt(nub * stmp5->stat6);
        else stmp5->stat6 = 0.0;
        sbb->stat6 = stmp4->stat6 - stmp1->stat6;
        stmp3->stat6 = (float)cnorminv((long double) (stmp3->stat6));

        stmp3->stat8 /= (float) nsamp;
        stmp4->stat8 /= (float) nsamp;
        stmp5->stat8 /= (float) nsamp;
        stmp5->stat8 -= stmp4->stat8 * stmp4->stat8;
        if(stmp5->stat8 > 0) stmp5->stat8 = sqrt(nub * stmp5->stat8);
        else stmp5->stat8 = 0.0;
        sbb->stat8 = stmp4->stat8 - stmp1->stat8;
        stmp3->stat8 = (float)cnorminv((long double) (stmp3->stat8));

        stmp3->stat9 /= (float) nsamp;
        stmp4->stat9 /= (float) nsamp;
        stmp5->stat9 /= (float) nsamp;
        stmp5->stat9 -= stmp4->stat9 * stmp4->stat9;
        if(stmp5->stat9 > 0) stmp5->stat9 = sqrt(nub * stmp5->stat9);
        else stmp5->stat9 = 0.0;
        sbb->stat9 = stmp4->stat9 - stmp1->stat9;
        stmp3->stat9 = (float)cnorminv((long double) (stmp3->stat9));

        stmp3->stat10 /= (float) nsamp;
        stmp4->stat10 /= (float) nsamp;
        stmp5->stat10 /= (float) nsamp;
        stmp5->stat10 -= stmp4->stat10 * stmp4->stat10;
        if(stmp5->stat10 > 0) stmp5->stat10 = sqrt(nub * stmp5->stat10);
        else stmp5->stat10 = 0.0;
        sbb->stat10 = stmp4->stat10 - stmp1->stat10;
        stmp3->stat10 = (float)cnorminv((long double) (stmp3->stat10));

        stmp3->stat12 /= (float) nsamp;
        stmp4->stat12 /= (float) nsamp;
        stmp5->stat12 /= (float) nsamp;
        stmp5->stat12 -= stmp4->stat12 * stmp4->stat12;
        if(stmp5->stat12 > 0) stmp5->stat12 = sqrt(nub * stmp5->stat12);
        else stmp5->stat12 = 0.0;
        sbb->stat12 = stmp4->stat12 - stmp1->stat12;
        stmp3->stat12 = (float)cnorminv((long double) (stmp3->stat12));

        stmp3->stat13 /= (float) nsamp;
        stmp4->stat13 /= (float) nsamp;
        stmp5->stat13 /= (float) nsamp;
        stmp5->stat13 -= stmp4->stat13 * stmp4->stat13;
        if(stmp5->stat13 > 0) stmp5->stat13 = sqrt(nub * stmp5->stat13);
        else stmp5->stat13 = 0.0;
        sbb->stat13 = stmp4->stat13 - stmp1->stat13;
        stmp3->stat13 = (float)cnorminv((long double) (stmp3->stat13));

    }

    fout = fopen(stfilo, "w");
    if(!fout){
       printf("****ERROR****, unable to open file %s for write.\n\n", stfilo);
       exit(1);
    }
    statdmp(fout, sout);
    fclose(fout);

    fout = fopen(sflerr, "w");
    if(!fout){
       printf("****ERROR****, unable to open file %s for write.\n\n", sflerr);
       exit(1);
    }
    statdmp(fout, serr);
    fclose(fout);

    printf("Standard confidence intervals, 'y' or 'n'.\n\n");
    scanf("\n");
    if(getchar() == 'y'){
       printf("What is the confidence interval required, e.g. 95%%.\n\n");
       scanf("%f", &tlev);
       printf("Use bias correction, '0' for no and '1' for yes.\n\n");
       scanf("%d", &ibias);
       tlim = 1.0 - tlev * 0.01;
       tlim *= 0.5;
       tlim = 1.0 - tlim;
       zfac = (float)cnorminv((long double) tlim);
/* re-use sout and smean */
       s_lo = sout;
       s_hi = smean;
       for(i=0; i < ptnum; i++){
           stmp1 = stato->ptst + i;
           stmp2 = sbias->ptst + i;
           stmp3 = serr->ptst + i;
           stmp4 = s_lo->ptst + i;
           stmp5 = s_hi->ptst + i;

           stat = (ibias) ? (stmp1->stat1).mean - (stmp2->stat1).mean : (stmp1->stat1).mean;
           (stmp4->stat1).mean = stat - zfac * (stmp3->stat1).mean;
           (stmp5->stat1).mean = stat + zfac * (stmp3->stat1).mean;

           stat = (ibias) ? (stmp1->stat1).var - (stmp2->stat1).var : (stmp1->stat1).var;
           (stmp4->stat1).var = stat - zfac * (stmp3->stat1).var;
           (stmp5->stat1).var = stat + zfac * (stmp3->stat1).var;

           stat = (ibias) ? (stmp1->stat2).mean - (stmp2->stat2).mean : (stmp1->stat2).mean;
           (stmp4->stat2).mean = stat - zfac * (stmp3->stat2).mean;
           (stmp5->stat2).mean = stat + zfac * (stmp3->stat2).mean;

           stat = (ibias) ? (stmp1->stat2).var - (stmp2->stat2).var : (stmp1->stat2).var;
           (stmp4->stat2).var = stat - zfac * (stmp3->stat2).var;
           (stmp5->stat2).var = stat + zfac * (stmp3->stat2).var;

           stat = (ibias) ? stmp1->stat3 - stmp2->stat3 : stmp1->stat3;
           stmp4->stat3 = stat - zfac * stmp3->stat3;
           stmp5->stat3 = stat + zfac * stmp3->stat3;

           stat = (ibias) ? stmp1->stat4 - stmp2->stat4 : stmp1->stat4;
           stmp4->stat4 = stat - zfac * stmp3->stat4;
           stmp5->stat4 = stat + zfac * stmp3->stat4;

           stat = (ibias) ? stmp1->stat5 - stmp2->stat5 : stmp1->stat5;
           stmp4->stat5 = stat - zfac * stmp3->stat5;
           stmp5->stat5 = stat + zfac * stmp3->stat5;

           stat = (ibias) ? stmp1->stat6 - stmp2->stat6 : stmp1->stat6;
           stmp4->stat6 = stat - zfac * stmp3->stat6;
           stmp5->stat6 = stat + zfac * stmp3->stat6;

           stat = (ibias) ? stmp1->stat8 - stmp2->stat8 : stmp1->stat8;
           stmp4->stat8 = stat - zfac * stmp3->stat8;
           stmp5->stat8 = stat + zfac * stmp3->stat8;

           stat = (ibias) ? stmp1->stat9 - stmp2->stat9 : stmp1->stat9;
           stmp4->stat9 = stat - zfac * stmp3->stat9;
           stmp5->stat9 = stat + zfac * stmp3->stat9;

           stat = (ibias) ? stmp1->stat10 - stmp2->stat10 : stmp1->stat10;
           stmp4->stat10 = stat - zfac * stmp3->stat10;
           stmp5->stat10 = stat + zfac * stmp3->stat10;

           stat = (ibias) ? stmp1->stat12 - stmp2->stat12 : stmp1->stat12;
           stmp4->stat12 = stat - zfac * stmp3->stat12;
           stmp5->stat12 = stat + zfac * stmp3->stat12;

           stat = (ibias) ? stmp1->stat13 - stmp2->stat13 : stmp1->stat13;
           stmp4->stat13 = stat - zfac * stmp3->stat13;
           stmp5->stat13 = stat + zfac * stmp3->stat13;

       }

       fout = fopen("low_std.dat", "w");
       if(!fout){
          printf("****ERROR****, unable to open file %s for write.\n\n", "low_std.dat");
          exit(1);
       }
       statdmp(fout, s_lo);
       fclose(fout);

       fout = fopen("high_std.dat", "w");
       if(!fout){
          printf("****ERROR****, unable to open file %s for write.\n\n", "high_std.dat");
          exit(1);
       }
       statdmp(fout, s_hi);
       fclose(fout);

       printf("Do you want the length and shape diagnostics, 'y' or 'n'.\n\n");
       scanf("\n");
       if(getchar() == 'y'){

/* re-use serr and sbias */

          slen = sbias;
          sshp = serr;

          for(i=0; i < ptnum; i++){
             stmp1 = stato->ptst + i;
             stmp2 = s_lo->ptst + i;
             stmp3 = s_hi->ptst + i;
             stmp4 = slen->ptst + i;
             stmp5 = sshp->ptst + i;

             (stmp4->stat1).mean = (stmp3->stat1).mean - (stmp2->stat1).mean;
             if((stmp1->stat1).mean > 0.0) 
                (stmp5->stat1).mean = ((stmp3->stat1).mean - (stmp1->stat1).mean) / ((stmp1->stat1).mean - (stmp2->stat1).mean);
             else (stmp5->stat1).mean = 0.0;

             (stmp4->stat1).var = (stmp3->stat1).var - (stmp2->stat1).var;
             if((stmp1->stat1).var > 0.0) 
                (stmp5->stat1).var = ((stmp3->stat1).var - (stmp1->stat1).var) / ((stmp1->stat1).var - (stmp2->stat1).var);
             else (stmp5->stat1).var = 0.0;

             (stmp4->stat2).mean = (stmp3->stat2).mean - (stmp2->stat2).mean;
             if((stmp1->stat2).mean > 0.0) 
                (stmp5->stat2).mean = ((stmp3->stat2).mean - (stmp1->stat2).mean) / ((stmp1->stat2).mean - (stmp2->stat2).mean);
             else (stmp5->stat2).mean = 0.0;

             (stmp4->stat2).var = (stmp3->stat2).var - (stmp2->stat2).var;
             if((stmp1->stat2).var > 0.0) 
                (stmp5->stat2).var = ((stmp3->stat2).var - (stmp1->stat2).var) / ((stmp1->stat2).var - (stmp2->stat2).var);
             else (stmp5->stat2).var = 0.0;

             stmp4->stat3 = stmp3->stat3 - stmp2->stat3;
             if(stmp1->stat3 > 0.0) 
                stmp5->stat3 = (stmp3->stat3 - stmp1->stat3) / (stmp1->stat3 - stmp2->stat3);
             else stmp5->stat3 = 0.0;

             stmp4->stat4 = stmp3->stat4 - stmp2->stat4;
             if(stmp1->stat4 > 0.0) 
                stmp5->stat4 = (stmp3->stat4 - stmp1->stat4) / (stmp1->stat4 - stmp2->stat4);
             else stmp5->stat4 = 0.0;

             stmp4->stat5 = stmp3->stat5 - stmp2->stat5;
             if(stmp1->stat5 > 0.0) 
                stmp5->stat5 = (stmp3->stat5 - stmp1->stat5) / (stmp1->stat5 - stmp2->stat5);
             else stmp5->stat5 = 0.0;

             stmp4->stat6 = stmp3->stat6 - stmp2->stat6;
             if(stmp1->stat6 > 0.0) 
                stmp5->stat6 = (stmp3->stat6 - stmp1->stat6) / (stmp1->stat6 - stmp2->stat6);
             else stmp5->stat6 = 0.0;

             stmp4->stat8 = stmp3->stat8 - stmp2->stat8;
             if(stmp1->stat8 > 0.0) 
                stmp5->stat8 = (stmp3->stat8 - stmp1->stat8) / (stmp1->stat8 - stmp2->stat8);
             else stmp5->stat8 = 0.0;

             stmp4->stat9 = stmp3->stat9 - stmp2->stat9;
             if(stmp1->stat9 > 0.0) 
                stmp5->stat9 = (stmp3->stat9 - stmp1->stat9) / (stmp1->stat9 - stmp2->stat9);
             else stmp5->stat9 = 0.0;

             stmp4->stat10 = stmp3->stat10 - stmp2->stat10;
             if(stmp1->stat10 > 0.0) 
                stmp5->stat10 = (stmp3->stat10 - stmp1->stat10) / (stmp1->stat10 - stmp2->stat10);
             else stmp5->stat10 = 0.0;

             stmp4->stat12 = stmp3->stat12 - stmp2->stat12;
             if(stmp1->stat12 > 0.0) 
                stmp5->stat12 = (stmp3->stat12 - stmp1->stat12) / (stmp1->stat12 - stmp2->stat12);
             else stmp5->stat12 = 0.0;

             stmp4->stat13 = stmp3->stat13 - stmp2->stat13;
             if(stmp1->stat13 > 0.0) 
                stmp5->stat13 = (stmp3->stat13 - stmp1->stat13) / (stmp1->stat13 - stmp2->stat13);
             else stmp5->stat13 = 0.0;

          }

          fout = fopen("len_std.dat", "w");
          if(!fout){
             printf("****ERROR****, unable to open file %s for write.\n\n", "len_std.dat");
             exit(1);
          }
          statdmp(fout, slen);
          fclose(fout);

          fout = fopen("shp_std.dat", "w");
          if(!fout){
             printf("****ERROR****, unable to open file %s for write.\n\n", "shp_std.dat");
             exit(1);
          }
          statdmp(fout, sshp);
          fclose(fout);
         
       }

    }

    free(sout->ptst);
    free(sout);

    free(stato->ptst);
    free(stato);

    free(smean->ptst);
    free(smean);

    free(serr->ptst);
    free(serr);

    free(sbias->ptst);
    free(sbias);

    return 0;

}
