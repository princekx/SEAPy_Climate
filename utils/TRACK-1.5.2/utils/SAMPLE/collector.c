#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistic.h"
#include "mem_er.h"

#define  MAXCHR   100
#define  NSTAT    13

struct tot_stat *read_stats(FILE * );
void statdmp(FILE * , struct tot_stat * );
void sumprb(float , float * , float , float , float , int * , int , int , int );

int prgr, prty;

int main(void )
{
    int i, j, k;
    int nsamp=0, nnsmp=0;
    int inc=0;
    int ptnum, ptnn;
    int nstn;
    int nbin;
    int ii=0;
    int itest=0;
    int iside=1;
    int ibin=0;

    int **distb[NSTAT];

    char stub[MAXCHR], infile[MAXCHR], mxmn[MAXCHR];
    char prbf[MAXCHR], outfil[MAXCHR];
    char cnum[10], dum[20];
    char pv='t';

    struct tot_stat *stt=NULL;
    struct tot_stat *prb=NULL;
    struct tot_stat *sout1=NULL;
    struct pt_stat *stmp=NULL, *stmpo=NULL;

    float *smin[NSTAT], *smax[NSTAT];
    float *binw[NSTAT];
    float tlev=0.0;

    FILE *fin=NULL, *fout=NULL, *fprbf=NULL;

    printf("Input file stub for input files.\n\n");
    scanf("%s", stub);

    printf("How many samples are there to collect statistics on?\n\n");
    scanf("%d", &nsamp);

    nnsmp = nsamp - 1;

    printf("What is the maxmin file to read?\n\n");
    scanf("%s", mxmn);

    fin = fopen(mxmn, "r");
    if(fin == NULL){
       printf("****ERROR****, cant open file \r\n"
              "               %s\r\n"
              "               for read\n\n", mxmn);
       exit(1);
    }

    fscanf(fin, "%d %d", &ptnn, &nstn);

    if(nstn != NSTAT){
       printf("****ERROR****, the maxmin file has the wrong number of min-max\r\n"
              "               ranges for the defined number of statistics.   \n\n");
       exit(1);
    }

    for(i=0; i < NSTAT; i++){
       smin[i] = (float *)calloc(ptnn, sizeof(float));
       mem_er((smin[i] == NULL) ? 0 : 1, ptnn * sizeof(float));
       smax[i] = (float *)calloc(ptnn, sizeof(float));
       mem_er((smax[i] == NULL) ? 0 : 1, ptnn * sizeof(float));
       binw[i] = (float *)calloc(ptnn, sizeof(float));
       mem_er((binw[i] == NULL) ? 0 : 1, ptnn * sizeof(float));
    }

    for(i=0; i < ptnn; i++){
        for(j=0; j < NSTAT; j++)fscanf(fin, "%f %f", smin[j] + i, smax[j] + i);
        fgets(dum, 20, fin);
    }

    fclose(fin);

    printf("How many bins are required for the distributions?\n\n");
    scanf("%d", &nbin);

    for(i=0; i < NSTAT; i++){
       distb[i] = (int **)calloc(ptnn, sizeof(int *));
       mem_er((distb[i] == NULL) ? 0 : 1, ptnn * sizeof(int *));
       for(j=0; j < ptnn; j++){
           *(distb[i] + j) = (int *)calloc(nbin+1, sizeof(int));
           mem_er((distb[i] + j == NULL) ? 0 : 1, (nbin+1) * sizeof(int));
       }
    }

    for(i=0; i<ptnn; i++){

        for(j=0; j < NSTAT; j++){

            *(binw[j] + i) = (*(smax[j] + i) - *(smin[j] + i)) / nbin;
            *(smin[j] + i) -= *(binw[j] + i) / 2.0;
            *(smax[j] + i) += *(binw[j] + i) / 2.0;

        }


    }

    ++nbin;

    printf("What is the statistics file to test?\n\n");
    scanf("%s", prbf);

    fprbf = fopen(prbf, "r");
    if(fprbf == NULL){
       printf("****ERROR****, cant open file \r\n"
              "               %s\r\n"
              "               for read\n\n", prbf);
       exit(1);
    }

    prb = read_stats(fprbf);

    fclose(fprbf);


    printf("Is this a confidence interval estimation, '0' or a significance test, '1'\n\n");
    scanf("%d", &itest);

    if(!itest){

       printf("Specify the level for the confidence intervals, e.g. 95%%\n\n");
       scanf("%f", &tlev);

    }

    else{
       printf("Do you want a single sided test, '0' or a double sided test '1'\n\n");
       scanf("%d", &iside);

       printf("Do you want to compute p-values 'p' or just a simple test 't'\n\n");
       scanf("\n");
       if((pv=getchar()) == 't'){
          printf("Specify the level for the test, e.g. 95%%\n\n");
          scanf("%f", &tlev);
       }

    }
    
    while(inc <= nnsmp){

       sprintf(cnum, ".%d", inc);
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

       stt = read_stats(fin);

       fclose(fin);

       if(!inc){
          ptnum = stt->ptnum;
          if(ptnum != ptnn){
             printf("****ERROR****, the number of points in the maxmin and statistics\r\n"
                    "               files are inconsistent.                          \n\n");
             exit(1);
          }
       }


       for(i=0; i < ptnum; i++){
           stmp = stt->ptst + i;

           if(*(binw[0] + i) > 0.){
             ii = (int)(((stmp->stat1).mean - *(smin[0] + i)) / *(binw[0] + i));
             ++(*(*(distb[0] + i) + ii));
           }
           if(*(binw[1] + i) > 0.){
              ii = (int)(((stmp->stat1).var - *(smin[1] + i)) / *(binw[1] + i));
              ++(*(*(distb[1] + i) + ii));
           }
           if(*(binw[2] + i) > 0.){
              ii = (int)(((stmp->stat2).mean - *(smin[2] + i)) / *(binw[2] + i));
              ++(*(*(distb[2] + i) + ii));
           }
           if(*(binw[3] + i) > 0.){
              ii = (int)(((stmp->stat2).var - *(smin[3] + i)) / *(binw[3] + i));
              ++(*(*(distb[3] + i) + ii));
           }
           if(*(binw[4] + i) > 0.){
              ii = (int)((stmp->stat3 - *(smin[4] + i)) / *(binw[4] + i));
              ++(*(*(distb[4] + i) + ii));
           }
           if(*(binw[5] + i) > 0.){
              ii = (int)((stmp->stat4 - *(smin[5] + i)) / *(binw[5] + i));
              ++(*(*(distb[5] + i) + ii));
           }
           if(*(binw[6] + i) > 0.){
              ii = (int)((stmp->stat5 - *(smin[6] + i)) / *(binw[6] + i));
              ++(*(*(distb[6] + i) + ii));
           }
           if(*(binw[7] + i) > 0.){
              ii = (int)((stmp->stat6 - *(smin[7] + i)) / *(binw[7] + i));
              ++(*(*(distb[7] + i) + ii));
           }
           if(*(binw[8] + i) > 0.){
              ii = (int)((stmp->stat8 - *(smin[8] + i)) / *(binw[8] + i));
              ++(*(*(distb[8] + i) + ii));
           }
           if(*(binw[9] + i) > 0.){
              ii = (int)((stmp->stat9 - *(smin[9] + i)) / *(binw[9] + i));
              ++(*(*(distb[9] + i) + ii));
           }
           if(*(binw[10] + i) > 0.){
              ii = (int)((stmp->stat10 - *(smin[10] + i)) / *(binw[10] + i));
              ++(*(*(distb[10] + i) + ii));
           }
           if(*(binw[11] + i) > 0.){
              ii = (int)((stmp->stat12 - *(smin[11] + i)) / *(binw[11] + i));
              ++(*(*(distb[11] + i) + ii));
           }
           if(*(binw[12] + i) > 0.){
              ii = (int)((stmp->stat13 - *(smin[12] + i)) / *(binw[12] + i));
              ++(*(*(distb[12] + i) + ii));
           }
       }


       free(stt->ptst);
       free(stt);

       ++inc;

    }


/*for(k=0; k < NSTAT; k++){
for(i=0; i < nbin; i++){
    printf("%7.4f ", *(smin[k]) + ((float)i * *(binw[k]) + (*(binw[k]) / (float)nbin)));
    for(j=0; j < *(*(distb[k]) + i); j++) printf("#");
    printf("\n");

}
printf("\n\n");
} */

/* assign memory for test output */

    sout1 = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((sout1 == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(sout1, prb, sizeof(struct tot_stat));

    sout1->ptnum = ptnum;
    sout1->datnm[8] = 0;
    sout1->datnm[12] = 0;

    sout1->ptst = (struct pt_stat * )calloc(ptnum, sizeof(struct pt_stat));
    mem_er((sout1->ptst == NULL) ? 0 : 1, ptnum * sizeof(struct pt_stat));


    if(!itest){

    }
    else {

        for(i=0; i < ptnum; i++){
            stmp = prb->ptst + i;
            stmpo = sout1->ptst + i;

            stmpo->ptyp = stmp->ptyp;
            stmpo->xs = stmp->xs;
            stmpo->ys = stmp->ys;

            if(*(binw[0] + i) > 0.) 
               sumprb((stmp->stat1).mean, &((stmpo->stat1).mean), *(smin[0] + i), *(smax[0] + i), *(binw[0] + i), *(distb[0] + i), nbin, iside, nsamp);
            else (stmpo->stat1).mean = 1.0;

            if(*(binw[1] + i) > 0.)            
               sumprb((stmp->stat1).var, &((stmpo->stat1).var), *(smin[1] + i), *(smax[1] + i), *(binw[1] + i), *(distb[1] + i), nbin, iside, nsamp);
            else (stmpo->stat1).var = 1.0;

            if(*(binw[2] + i) > 0.) 
               sumprb((stmp->stat2).mean, &((stmpo->stat2).mean), *(smin[2] + i), *(smax[2] + i), *(binw[2] + i), *(distb[2] + i), nbin, iside, nsamp);
            else (stmpo->stat2).mean = 1.0;

            if(*(binw[3] + i) > 0.)            
               sumprb((stmp->stat2).var, &((stmpo->stat2).var), *(smin[3] + i), *(smax[3] + i), *(binw[3] + i), *(distb[3] + i), nbin, iside, nsamp);
            else (stmpo->stat2).var = 1.0;

            if(*(binw[4] + i) > 0.)            
               sumprb(stmp->stat3, &(stmpo->stat3), *(smin[4] + i), *(smax[4] + i), *(binw[4] + i), *(distb[4] + i), nbin, iside, nsamp);
            else stmpo->stat3 = 1.0;

            if(*(binw[5] + i) > 0.)            
               sumprb(stmp->stat4, &(stmpo->stat4), *(smin[5] + i), *(smax[5] + i),  *(binw[5] + i), *(distb[5] + i), nbin, iside, nsamp);
            else stmpo->stat4 = 1.0;

            if(*(binw[6] + i) > 0.)            
               sumprb(stmp->stat5, &(stmpo->stat5), *(smin[6] + i), *(smax[6] + i), *(binw[6] + i), *(distb[6] + i), nbin, iside, nsamp);
            else stmpo->stat5 = 1.0;

            if(*(binw[7] + i) > 0.)            
               sumprb(stmp->stat6, &(stmpo->stat6), *(smin[7] + i), *(smax[7] + i), *(binw[7] + i), *(distb[7] + i), nbin, iside, nsamp);
            else stmpo->stat6 = 1.0;

/* printf("%f %f %f %f\n",  *(smin[7] + i), *(smax[7] + i), stmp->stat6, stmpo->stat6); */

            if(*(binw[8] + i) > 0.)            
               sumprb(stmp->stat8, &(stmpo->stat8), *(smin[8] + i), *(smax[8] + i), *(binw[8] + i), *(distb[8] + i), nbin, iside, nsamp);
            else stmpo->stat8 = 1.0;

            if(*(binw[9] + i) > 0.)            
               sumprb(stmp->stat9, &(stmpo->stat9), *(smin[9] + i), *(smax[9] + i), *(binw[9] + i), *(distb[9] + i), nbin, iside, nsamp);
            else stmpo->stat9 = 1.0;

            if(*(binw[10] + i) > 0.)            
               sumprb(stmp->stat10, &(stmpo->stat10), *(smin[10] + i), *(smax[10] + i), *(binw[10] + i), *(distb[10] + i), nbin, iside, nsamp);
            else stmpo->stat10 = 1.0;

            if(*(binw[11] + i) > 0.)            
               sumprb(stmp->stat12, &(stmpo->stat12), *(smin[11] + i), *(smax[11] + i), *(binw[11] + i), *(distb[11] + i), nbin, iside, nsamp);
            else stmpo->stat12 = 1.0;

            if(*(binw[12] + i) > 0.)            
               sumprb(stmp->stat13, &(stmpo->stat13), *(smin[12] + i), *(smax[12] + i), *(binw[12] + i), *(distb[12] + i), nbin, iside, nsamp);
            else stmpo->stat13 = 1.0;

        }

    }


/* write statistics to file */


    printf("What output file do you want to write the test information?\n\n");
    scanf("%s", outfil);

    fout = fopen(outfil, "w");
    statdmp(fout, sout1);
    fclose(fout);


    free(sout1->ptst);
    free(sout1);

    for(i=0; i < NSTAT; i++){
        for(j=0; j < ptnn; j++) free(*(distb[i] + j));
        free(distb[i]);
    }

    free(prb->ptst);
    free(prb);

    for(i=0; i < NSTAT; i++){
        free(smin[i]);
        free(smax[i]);
        free(binw[i]);
    }

    return 0;

}

void sumprb(float stat, float *sout, float smin, float smax, float binw, int *distb, int nbin, int iside, int nsamp)
{

    int i;
    int ibin=0;
    int ipp=0;
    int ifval=0;

    float sump=0.0;

    if(stat < smin || stat > smax) {
       *sout = 0.0;
       return;
    }

    ibin = (int)((stat - smin) / binw);
    sump = 0.0;

    if(stat < 0.0){

       for(i=0; i <= ibin; i++) sump += *(distb + i);

    }

    else {

       for(i=ibin; i < nbin; i++) sump += *(distb + i);
       ipp = 1;
    }

    *sout = sump / nsamp;


    if(iside){

       ifval = *(distb + ibin);

       if(!ipp){

          for(i=nbin - 1; i > ibin; i--) {
              if(*(distb + i) >= ifval) {ibin = i; break;}
          }

          sump = 0.0;

          for(i=ibin; i < nbin; i++) sump += *(distb + i);
          *sout += sump / nsamp;

       }

       else {

          for(i=0; i < ibin; i++){
             if(*(distb + i) >= ifval) {ibin = i; break;}
          }

          sump = 0.0;
          for(i=0; i < ibin; i++) sump += *(distb + i);
          *sout += sump / nsamp;

       }


    }

/*    else{

       if(stat < 0.0) *sout *= -1.0;

    } */

    return;

}
