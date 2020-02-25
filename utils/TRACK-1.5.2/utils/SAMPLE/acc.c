#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "statistic.h"
#include "statistic_double.h"
#include "mem_er.h"

#define  MAXCHR  200
#define  LINE    40000
#define  NSTAT   13
#define  TOLBOT  1.0e-6

struct tot_stat *read_stats(FILE * );
void statdmp(FILE * , struct tot_stat * );

int prgr, prty;

int main(void)
{


    int i, j;
    int ptnum=0;
    int ntrack=0;
    int isamp=0;
    int nb=0, nbb=0;
    long int nbtot=0;

    int iret=0;

    double normb=0.0;
    double wi=0;

    double s1, s2, s3;
    double top, bot;

    char stub[MAXCHR], infile[MAXCHR];
    char cnum[10];
    char sflacc[]="stats_acc.dat";
    char missf[]="conv_missing.dat";
    char line[LINE], *ll=NULL;

    FILE *fin=NULL, *fout=NULL;
    FILE *fstat=NULL, *fmiss=NULL;

    struct tot_stat_d *sb=NULL;
    struct tot_stat *sacc=NULL;
    struct tot_stat *sin=NULL;
    struct tot_stat_d *ssum1=NULL, *ssum2=NULL, *ssum3=NULL;
    struct pt_stat_d *pts1=NULL, *pts2=NULL, *pts3=NULL;
    struct pt_stat_d *sbb=NULL;
    struct pt_stat *sac=NULL, *psin=NULL;

    fmiss = fopen(missf, "r");
    if(!fmiss){
       printf("****ERROR****, opening file %s for read.\n\n", missf);
       exit(1);
    }
    fgets(line, LINE, fmiss);
    sscanf(line, "%d %ld", &ntrack, &nbtot);

    normb = (float)nbtot / (float)ntrack;

    printf("normb %f\n", normb);

    printf("Input file stub for input files.\n\n");
    scanf("%s", stub);

    sacc = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((sacc == NULL) ? 0 : 1, sizeof(struct tot_stat));

    ssum1 = (struct tot_stat_d * )malloc_initl(sizeof(struct tot_stat_d));
    mem_er((ssum1 == NULL) ? 0 : 1, sizeof(struct tot_stat_d));

    ssum2 = (struct tot_stat_d * )malloc_initl(sizeof(struct tot_stat_d));
    mem_er((ssum2 == NULL) ? 0 : 1, sizeof(struct tot_stat_d));

    ssum3 = (struct tot_stat_d * )malloc_initl(sizeof(struct tot_stat_d));
    mem_er((ssum3 == NULL) ? 0 : 1, sizeof(struct tot_stat_d)); 

ntrack=5;

    for(i=0; i < ntrack; i++){
        nbb = 0;

        printf("track %d\n", i);

        fgets(line, LINE, fmiss);
        ll = line;
        sscanf(ll, "%d", &nb);

        ll += 7;

        while((iret=sscanf(ll, "%d", &isamp)) != EOF){

             ++nbb;
             if(isamp < 10) sprintf(cnum, "00000%d", isamp);
             else if(isamp < 100) sprintf(cnum, "0000%d", isamp);
             else if(isamp < 1000) sprintf(cnum, "000%d", isamp);
             else if(isamp < 10000) sprintf(cnum, "00%d", isamp);
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

             if(! sacc->ptst){
                memcpy(sacc, sin, sizeof(struct tot_stat));
                sacc->ptnum = sin->ptnum;
             }

             if(!sb){

                sb = (struct tot_stat_d * )malloc_initl(sizeof(struct tot_stat_d));
                mem_er((sb == NULL) ? 0 : 1, sizeof(struct tot_stat_d));

                ptnum = sb->ptnum = sin->ptnum;

                sb->ptst = (struct pt_stat_d *)calloc(ptnum, sizeof(struct pt_stat_d));
                mem_er((sb->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat_d));
                for(j=0; j < ptnum; j++){
                    psin = sin->ptst + j;
                    sbb = sb->ptst + j;
                    sbb->xs = psin->xs;
                    sbb->ys = psin->ys;
                }
             }

             for(j=0; j < ptnum; j++){
                psin = sin->ptst + j;
                sbb = sb->ptst + j;
                (sbb->stat1).mean += (psin->stat1).mean;
                (sbb->stat1).var += (psin->stat1).var;
                (sbb->stat2).mean += (psin->stat2).mean;
                (sbb->stat2).var += (psin->stat2).var;
                sbb->stat3 += psin->stat3;
                sbb->stat4 += psin->stat4;
                sbb->stat5 += psin->stat5;
                sbb->stat6 += psin->stat6;
                sbb->stat8 += psin->stat8;
                sbb->stat9 += psin->stat9;
                sbb->stat10 += psin->stat10;
                sbb->stat12 += psin->stat12;
                sbb->stat13 += psin->stat13;
             }

             free(sin->ptst);
             free(sin);

             ll = strchr(ll, ' ') + 1;

        }

        for(j=0; j < ptnum; j++){
            sbb = sb->ptst + j;
            (sbb->stat1).mean /= nb;
            (sbb->stat1).var /= nb;
            (sbb->stat2).mean /= nb;
            (sbb->stat2).var /= nb;
            sbb->stat3 /= nb;
            sbb->stat4 /= nb;
            sbb->stat5 /= nb;
            sbb->stat6 /= nb;
            sbb->stat8 /= nb;
            sbb->stat9 /= nb;
            sbb->stat10 /= nb;
            sbb->stat12 /= nb;
            sbb->stat13 /= nb;
        }

        if(!ssum1->ptst){
           ssum1->ptnum = ssum2->ptnum = ssum3->ptnum = ptnum;

           ssum1->ptst = (struct pt_stat_d *)calloc(ptnum, sizeof(struct pt_stat_d));
           mem_er((ssum1->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat_d));

           ssum2->ptst = (struct pt_stat_d *)calloc(ptnum, sizeof(struct pt_stat_d));
           mem_er((ssum2->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat_d));

           ssum3->ptst = (struct pt_stat_d *)calloc(ptnum, sizeof(struct pt_stat_d));
           mem_er((ssum3->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat_d));

           sacc->ptst = (struct pt_stat *)calloc(ptnum, sizeof(struct pt_stat));
           mem_er((sacc->ptst == NULL) ? 0 : 1, ptnum*sizeof(struct pt_stat));

           for(j=0; j < ptnum; j++){
               sbb = sb->ptst + j;
               sac = sacc->ptst + j;
               sac->xs = sbb->xs;
               sac->ys = sbb->ys;
           }           

        }

        wi = (float) nb / normb;

        for(j=0; j < ptnum; j++){
            sbb = sb->ptst + j;
            pts1 = ssum1->ptst + j;
            pts2 = ssum2->ptst + j;
            pts3 = ssum3->ptst + j;
            s1 = wi * (sbb->stat1).mean;
            s2 = s1 * s1;
            s3 = s2 * s1;
            (pts1->stat1).mean += s1;
            (pts2->stat1).mean += s2;
            (pts3->stat1).mean += s3;
            s1 = wi * (sbb->stat1).var;
            s2 = s1 * s1;
            s3 = s2 * s1;
            (pts1->stat1).var += s1;
            (pts2->stat1).var += s2;
            (pts2->stat1).var += s2;
            s1 = wi * (sbb->stat2).mean;
            s2 = s1 * s1;
            s3 = s2 * s1;
            (pts1->stat2).mean += s1;
            (pts2->stat2).mean += s2;
            (pts3->stat2).mean += s3;
            s1 = wi * (sbb->stat2).var;
            s2 = s1 * s1;
            s3 = s2 * s1;
            (pts1->stat2).var += s1;
            (pts2->stat2).var += s2;
            (pts2->stat2).var += s2;
            s1 = wi * sbb->stat3;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat3 += s1;
            pts2->stat3 += s2;
            pts3->stat3 += s3;
            s1 = wi * sbb->stat4;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat4 += s1;
            pts2->stat4 += s2;
            pts3->stat4 += s3;
            s1 = wi * sbb->stat5;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat5 += s1;
            pts2->stat5 += s2;
            pts3->stat5 += s3;
            s1 = wi * sbb->stat6;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat6 += s1;
            pts2->stat6 += s2;
            pts3->stat6 += s3;
            s1 = wi * sbb->stat8;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat8 += s1;
            pts2->stat8 += s2;
            pts3->stat8 += s3;
            s1 = wi * sbb->stat9;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat9 += s1;
            pts2->stat9 += s2;
            pts3->stat9 += s3;
            s1 = wi * sbb->stat10;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat10 += s1;
            pts2->stat10 += s2;
            pts3->stat10 += s3;
            s1 = wi * sbb->stat12;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat12 += s1;
            pts2->stat12 += s2;
            pts3->stat12 += s3;
            s1 = wi * sbb->stat13;
            s2 = s1 * s1;
            s3 = s2 * s1;
            pts1->stat13 += s1;
            pts2->stat13 += s2;
            pts3->stat13 += s3;
        }

        printf("%d %d\n", nb, nbb);
        if(nb != nbb){
           printf("****ERROR****, number of samples for track %d is not consitent.\n\n", i);
           exit(1);
        }

        free(sb->ptst);
        free(sb);

    }

    fclose(fmiss);

    for(i=0; i < ptnum; i++){
        sac = sacc->ptst + i;
        pts1 = ssum1->ptst + i;
        pts2 = ssum2->ptst + i;
        pts3 = ssum3->ptst + i;

        s1 = (pts1->stat1).mean / ntrack;
        s2 = (pts2->stat1).mean / ntrack;
        s3 = (pts3->stat1).mean / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) (sac->stat1).mean = top / (6.0 * bot);

        s1 = (pts1->stat1).var / ntrack;
        s2 = (pts2->stat1).var / ntrack;
        s3 = (pts3->stat1).var / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) (sac->stat1).var = top / (6.0 * bot);

        s1 = (pts1->stat2).mean / ntrack;
        s2 = (pts2->stat2).mean / ntrack;
        s3 = (pts3->stat2).mean / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) (sac->stat2).mean = top / (6.0 * bot);

        s1 = (pts1->stat2).var / ntrack;
        s2 = (pts2->stat2).var / ntrack;
        s3 = (pts3->stat2).var / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) (sac->stat2).var = top / (6.0 * bot);

        s1 = pts1->stat3 / ntrack;
        s2 = pts2->stat3 / ntrack;
        s3 = pts3->stat3 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat3 = top / (6.0 * bot);

        s1 = pts1->stat4 / ntrack;
        s2 = pts2->stat4 / ntrack;
        s3 = pts3->stat4 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat4 = top / (6.0 * bot);

        s1 = pts1->stat5 / ntrack;
        s2 = pts2->stat5 / ntrack;
        s3 = pts3->stat5 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat5 = top / (6.0 * bot);

        s1 = pts1->stat6 / ntrack;
        s2 = pts2->stat6 / ntrack;
        s3 = pts3->stat6 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat6 = top / (6.0 * bot);

        s1 = pts1->stat8 / ntrack;
        s2 = pts2->stat8 / ntrack;
        s3 = pts3->stat8 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat8 = top / (6.0 * bot);

        s1 = pts1->stat9 / ntrack;
        s2 = pts2->stat9 / ntrack;
        s3 = pts3->stat9 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat9 = top / (6.0 * bot);

        s1 = pts1->stat10 / ntrack;
        s2 = pts2->stat10 / ntrack;
        s3 = pts3->stat10 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat10 = top / (6.0 * bot);

        s1 = pts1->stat12 / ntrack;
        s2 = pts2->stat12 / ntrack;
        s3 = pts3->stat12 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat12 = top / (6.0 * bot);

        s1 = pts1->stat13 / ntrack;
        s2 = pts2->stat13 / ntrack;
        s3 = pts3->stat13 / ntrack;
        top = 3.0 * s1 * s2 - 2.0 * pow(s1, 3.0) - s3;
        bot = pow(s2 - s1 * s1, 1.5);
        if(bot > TOLBOT) sac->stat13 = top / (6.0 * bot);
    }

    fout = fopen("stats_acc.dat", "w");
    if(!fout){
       printf("****ERROR****, unable to open file %s for write.\n\n", "stats_acc.dat");
       exit(1);
    }
    statdmp(fout, sacc);
    fclose(fout);

    free(sacc->ptst);
    free(sacc);
    free(ssum1->ptst);
    free(ssum1);
    free(ssum2->ptst);
    free(ssum2);
    free(ssum3->ptst);
    free(ssum3);

    return 0;

}
