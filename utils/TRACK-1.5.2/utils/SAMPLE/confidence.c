#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistic.h"
#include "mem_er.h"

#define  MAXCHR   100
#define  NSTAT    13

struct tot_stat *read_stats(FILE * );
void statdmp(FILE * , struct tot_stat * );
void clim(int * , int * , int * , int , float );
float ran3(int * );

int prgr, prty;

int main(void )
{
    int i, j;
    int nsamp=0, nnsmp=0;
    int inc=0;
    int ptnum, ptnn;
    int nstn;
    int nbin;
    int ii=0;
    int low=0, high=0;
    int iord=0, idum=0, nfile=0;
    int nf=0;

    int **distb[NSTAT];
    int *ifil=NULL;

    char stub[MAXCHR], infile[MAXCHR], mxmn[MAXCHR], stato[MAXCHR];
    char outfil1[MAXCHR], outfil2[MAXCHR];
    char cnum[10], dum[20];

    struct tot_stat *stt=NULL, *prb=NULL, *sto=NULL;
    struct tot_stat *sout1=NULL, *sout2=NULL;
    struct pt_stat *stmp=NULL, *stmp1=NULL, *stmp2=NULL;

    float *smin[NSTAT], *smax[NSTAT];
    float *binw[NSTAT];
    float tlev=0.0;
    float tlim=0.0;

    float stlen=0.0, stshp=0.0;

    FILE *fin=NULL, *fout=NULL;

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
            *(smin[j] + i) -= (*(binw[j] + i) * 0.5);
            *(smax[j] + i) += (*(binw[j] + i) * 0.5);

        }

    }

    ++nbin;


    printf("Specify the level for the confidence intervals, e.g. 95%%\n\n");
    scanf("%f", &tlev);

    tlim = 1.0 - tlev * 0.01;
    tlim *= 0.5 * nsamp;

    prb = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((prb == NULL) ? 0 : 1, sizeof(struct tot_stat));

    printf("Do you want sequentially chosen files or randomly chosen files, \r\n"
           "input '0' for sequential and '1' for random.                    \n\n");
    scanf("%d", &iord);

    if(iord){
       printf("Input a seed value, must be negative\n\n");
       scanf("%d", &idum);
       printf("How many files are there to sample from?\n\n");
       scanf("%d", &nfile);
       ifil = (int *)calloc(nfile, sizeof(int));
       mem_er((ifil == NULL) ? 0 : 1, nfile*sizeof(int));
       for(i=0; i < nfile; i++) *(ifil + i) = 1;
    }
    
    while(inc <= nnsmp){

       if(iord) {
          nf = (int)(ran3(&idum) * (float) (nfile - 1) + 0.5);
          if( ! (*(ifil + nf))) continue;
          else *(ifil + nf) = 0;
       }
       else nf = inc;

       if(nf < 10) sprintf(cnum, "00000%d", nf);
       else if(nf < 100) sprintf(cnum, "0000%d", nf);
       else if(nf < 1000) sprintf(cnum, "000%d", nf);
       else if(nf < 10000) sprintf(cnum, "00%d", nf);
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
          memcpy(prb, stt, sizeof(struct tot_stat));

          prb->ptst = (struct pt_stat * )calloc(ptnum, sizeof(struct pt_stat));
          mem_er((prb->ptst == NULL) ? 0 : 1, ptnum * sizeof(struct pt_stat));
          
          for(i=0; i<ptnum; i++) memcpy(prb->ptst + i, stt->ptst + i, sizeof(struct pt_stat));

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

    free(ifil);

/*for(k=0; k < NSTAT; k++){
for(i=0; i < nbin; i++){
    printf("%7.4f ", *(smin[k]) + ((float)i * *(binw[k]) + (*(binw[k]) / (float)nbin)));
    for(j=0; j < *(*(distb[k]) + i); j++) printf("#");
    printf("\n");

}
printf("\n\n");
} */

/* assign memory for output */

    sout1 = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((sout1 == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(sout1, prb, sizeof(struct tot_stat));

    sout1->ptnum = ptnum;
    sout1->datnm[8] = 0;
    sout1->datnm[12] = 0;

    sout1->ptst = (struct pt_stat * )calloc(ptnum, sizeof(struct pt_stat));
    mem_er((sout1->ptst == NULL) ? 0 : 1, ptnum * sizeof(struct pt_stat));

    sout2 = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((sout2 == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(sout2, prb, sizeof(struct tot_stat));

    sout2->ptnum = ptnum;
    sout2->datnm[8] = 0;
    sout2->datnm[12] = 0;

    sout2->ptst = (struct pt_stat * )calloc(ptnum, sizeof(struct pt_stat));
    mem_er((sout2->ptst == NULL) ? 0 : 1, ptnum * sizeof(struct pt_stat));


    for(i=0; i < ptnum; i++){
        stmp = prb->ptst + i;
        stmp1 = sout1->ptst + i;

        stmp1->ptyp = stmp->ptyp;
        stmp1->xs = stmp->xs;
        stmp1->ys = stmp->ys;

        stmp2 = sout2->ptst + i;

        stmp2->ptyp = stmp->ptyp;
        stmp2->xs = stmp->xs;
        stmp2->ys = stmp->ys;

        if(*(binw[0] + i) > 0.) {
           clim(&low, &high , *(distb[0] + i), nbin, tlim);
           (stmp1->stat1).mean = *(smin[0] + i) + ((float)low + 0.5) * *(binw[0] + i);
           (stmp2->stat1).mean = *(smin[0] + i) + ((float)high + 0.5) * *(binw[0] + i);

        }

        if(*(binw[1] + i) > 0.) {
           clim(&low, &high , *(distb[1] + i), nbin, tlim);
           (stmp1->stat1).var = *(smin[1] + i) + ((float)low + 0.5) * *(binw[1] + i);
           (stmp2->stat1).var = *(smin[1] + i) + ((float)high + 0.5) * *(binw[1] + i);
        }

        if(*(binw[2] + i) > 0.) {
           clim(&low, &high , *(distb[2] + i), nbin, tlim);
           (stmp1->stat2).mean = *(smin[2] + i) + ((float)low + 0.5) * *(binw[2] + i);
           (stmp2->stat2).mean = *(smin[2] + i) + ((float)high + 0.5) * *(binw[2] + i);
        }

        if(*(binw[3] + i) > 0.) {
           clim(&low, &high , *(distb[3] + i), nbin, tlim);
           (stmp1->stat2).var = *(smin[3] + i) + ((float)low + 0.5) * *(binw[3] + i);
           (stmp2->stat2).var = *(smin[3] + i) + ((float)high + 0.5) * *(binw[3] + i);
        }

        if(*(binw[4] + i) > 0.) {
           clim(&low, &high , *(distb[4] + i), nbin, tlim);
           stmp1->stat3 = *(smin[4] + i) + ((float)low + 0.5) * *(binw[4] + i);
           stmp2->stat3 = *(smin[4] + i) + ((float)high + 0.5) * *(binw[4] + i);
        }

        if(*(binw[5] + i) > 0.) {
           clim(&low, &high , *(distb[5] + i), nbin, tlim);
           stmp1->stat4 = *(smin[5] + i) + ((float)low + 0.5) * *(binw[5] + i);
           stmp2->stat4 = *(smin[5] + i) + ((float)high + 0.5) * *(binw[5] + i);
        }

        if(*(binw[6] + i) > 0.) {
           clim(&low, &high , *(distb[6] + i), nbin, tlim);
           stmp1->stat5 = *(smin[6] + i) + ((float)low + 0.5) * *(binw[6] + i);
           stmp2->stat5 = *(smin[6] + i) + ((float)high + 0.5) * *(binw[6] + i);
        }

        if(*(binw[7] + i) > 0.) {
           clim(&low, &high , *(distb[7] + i), nbin, tlim);
           stmp1->stat6 = *(smin[7] + i) + ((float)low + 0.5) * *(binw[7] + i);
           stmp2->stat6 = *(smin[7] + i) + ((float)high + 0.5) * *(binw[7] + i);
        }

        if(*(binw[8] + i) > 0.) {
           clim(&low, &high , *(distb[8] + i), nbin, tlim);
           stmp1->stat8 = *(smin[8] + i) + ((float)low + 0.5) * *(binw[8] + i);
           stmp2->stat8 = *(smin[8] + i) + ((float)high + 0.5) * *(binw[8] + i);
        }

        if(*(binw[9] + i) > 0.) {
           clim(&low, &high , *(distb[9] + i), nbin, tlim);
           stmp1->stat9 = *(smin[9] + i) + ((float)low + 0.5) * *(binw[9] + i);
           stmp2->stat9 = *(smin[9] + i) + ((float)high + 0.5) * *(binw[9] + i);
        }

        if(*(binw[10] + i) > 0.) {
           clim(&low, &high , *(distb[10] + i), nbin, tlim);
           stmp1->stat10 = *(smin[10] + i) + ((float)low + 0.5) * *(binw[10] + i);
           stmp2->stat10 = *(smin[10] + i) + ((float)high + 0.5) * *(binw[10] + i);
        }

        if(*(binw[11] + i) > 0.) {
           clim(&low, &high , *(distb[11] + i), nbin, tlim);
           stmp1->stat12 = *(smin[11] + i) + ((float)low + 0.5) * *(binw[11] + i);
           stmp2->stat12 = *(smin[11] + i) + ((float)high + 0.5) * *(binw[11] + i);
        }

        if(*(binw[12] + i) > 0.) {
           clim(&low, &high , *(distb[12] + i), nbin, tlim);
           stmp1->stat13 = *(smin[12] + i) + ((float)low + 0.5) * *(binw[12] + i);
           stmp2->stat13 = *(smin[12] + i) + ((float)high + 0.5) * *(binw[12] + i);
        }


    }


/* write statistics to file */


    printf("What two output files do you want to write the confidence interval information?\n\n");
    scanf("%s %s", outfil1, outfil2);

    fout = fopen(outfil1, "w");
    statdmp(fout, sout1);
    fclose(fout);

    fout = fopen(outfil2, "w");
    statdmp(fout, sout2);
    fclose(fout);

/* compute the length and shape of the confidence intervals */

    printf("Do you want to compute length and shape pf the confidence intervals? 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y'){

       printf("What is the original statistics file.\n\n");
       scanf("%s", stato);

       fin = fopen(stato, "r");
       if(fin == NULL){
          printf("****ERROR****, cant open file \r\n"
                 "               %s\r\n"
                 "               for read\n\n", stato);
          exit(1);
       }

       sto = read_stats(fin);

       fclose(fin);
       
       if(sto->ptnum != sout1->ptnum){
          printf("****ERROR****, stats number of points not consistent.\n\n");
          exit(1);
       }

/* reuse sout1 and sout2 */

       for(i=0; i < sout1->ptnum; i++){

          stmp = sto->ptst + i;
          stmp1 = sout1->ptst + i;
          stmp2 = sout2->ptst + i;

          stlen = (stmp2->stat1).mean - (stmp1->stat1).mean;
          if((stmp->stat1).mean > 0.){
             stshp = ((stmp2->stat1).mean - (stmp->stat1).mean) / ((stmp->stat1).mean - (stmp1->stat1).mean);
          }
          (stmp1->stat1).mean = stlen; 
          (stmp2->stat1).mean = stshp; 

          stlen = (stmp2->stat1).var - (stmp1->stat1).var;
          if((stmp->stat1).var > 0.){
             stshp = ((stmp2->stat1).var - (stmp->stat1).var) / ((stmp->stat1).var - (stmp1->stat1).var);
          }
          (stmp1->stat1).var = stlen;
          (stmp2->stat1).var = stshp;

          stlen = (stmp2->stat2).mean - (stmp1->stat2).mean;
          if((stmp->stat2).mean > 0.){
             stshp = ((stmp2->stat2).mean - (stmp->stat2).mean) / ((stmp->stat2).mean - (stmp1->stat2).mean);
          }
          (stmp1->stat2).mean = stlen;
          (stmp2->stat2).mean = stshp;

          stlen = (stmp2->stat2).var - (stmp1->stat2).var;
          if((stmp->stat2).var > 0.){
             stshp = ((stmp2->stat2).var - (stmp->stat2).var) / ((stmp->stat2).var - (stmp1->stat2).var);
          }
          (stmp1->stat2).var = stlen;
          (stmp2->stat2).var = stshp;

          stlen = stmp2->stat3 - stmp1->stat3;
          if(stmp->stat3 > 0.){
             stshp = (stmp2->stat3 - stmp->stat3) / (stmp->stat3 - stmp1->stat3);
          }
          stmp1->stat3 = stlen;
          stmp2->stat3 = stshp;

          stlen = stmp2->stat4 - stmp1->stat4;
          if(stmp->stat4 > 0.){
             stshp = (stmp2->stat4 - stmp->stat4) / (stmp->stat4 - stmp1->stat4);
          }
          stmp1->stat4 = stlen;
          stmp2->stat4 = stshp;

          stlen = stmp2->stat5 - stmp1->stat5;
          if(stmp->stat5 > 0.){
             stshp = (stmp2->stat5 - stmp->stat5) / (stmp->stat5 - stmp1->stat5);
          }
          stmp1->stat5 = stlen;
          stmp2->stat5 = stshp;

          stlen = stmp2->stat6 - stmp1->stat6;
          if(stmp->stat6 > 0.){
             stshp = (stmp2->stat6 - stmp->stat6) / (stmp->stat6 - stmp1->stat6);
          }
          stmp1->stat6 = stlen;
          stmp2->stat6 = stshp;

          stlen = stmp2->stat8 - stmp1->stat8;
          if(stmp->stat8 > 0.){
             stshp = (stmp2->stat8 - stmp->stat8) / (stmp->stat8 - stmp1->stat8);
          }
          stmp1->stat8 = stlen;
          stmp2->stat8 = stshp;

          stlen = stmp2->stat9 - stmp1->stat9;
          if(stmp->stat9 > 0.){
             stshp = (stmp2->stat9 - stmp->stat9) / (stmp->stat9 - stmp1->stat9);
          }
          stmp1->stat9 = stlen;
          stmp2->stat9 = stshp;

          stlen = stmp2->stat10 - stmp1->stat10;
          if(stmp->stat10 > 0.){
             stshp = (stmp2->stat10 - stmp->stat10) / (stmp->stat10 - stmp1->stat10);
          }
          stmp1->stat10 = stlen;
          stmp2->stat10 = stshp;

          stlen = stmp2->stat12 - stmp1->stat12;
          if(stmp->stat12 > 0.){
             stshp = (stmp2->stat12 - stmp->stat12) / (stmp->stat12 - stmp1->stat12);
          }
          stmp1->stat12 = stlen;
          stmp2->stat12 = stshp;

          stlen = stmp2->stat13 - stmp1->stat13;
          if(stmp->stat13 > 0.){
             stshp = (stmp2->stat13 - stmp->stat13) / (stmp->stat13 - stmp1->stat13);
          }
          stmp1->stat13 = stlen;
          stmp2->stat13 = stshp;


       }

       fout = fopen("conf_len.dat", "w");
       statdmp(fout, sout1);
       fclose(fout);

       fout = fopen("conf_shape.dat", "w");
       statdmp(fout, sout2);
       fclose(fout);

    }

    free(sout1->ptst);
    free(sout1);

    free(sout2->ptst);
    free(sout2);
    
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


void clim(int *low, int *high, int *dist, int nbin, float tlim)
{

     int i;

     float stat1, stat2;

     stat1 = *dist;
     *low = 0;
     if(stat1 < tlim){
       for(i=1; i < nbin; i++){
           if(stat1 + (float)(*(dist + i)) <= tlim){
              stat1 += (float)(*(dist + i));
              *low = i;
           }
           else break;

       }
     }

     stat2 = *(dist + nbin - 1);
     *high = nbin - 1;
     if(stat2 < tlim){
       for(i=1; i < nbin; i++){
           if(stat2 + (float)(*(dist + nbin - i - 1)) <= tlim){
           stat2 += (float)(*(dist + nbin - i - 1));
           *high = nbin - i - 1;
           }
           else break;
        }
     }

     return;
}
