#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "statistic.h"
#include "mem_er.h"
#include "grid.h"

#define  MAXCHR   200
#define  NSTAT    13

struct tot_stat *read_stats(FILE * );
void statdmp(FILE * , struct tot_stat * );
void read_gstat(struct tot_stat * , GRID * , char * , int * );
void netcdf_write_stats(struct tot_stat * , GRID * , char * , char * , int , int );

int prgr, prty;

int main(void )
{
    int i;
    int nsamp=0, nnsmp=0;
    int inc=0;
    int ptnum;
    int snum=0;
    int idd=0;

    float nns=0.0;

    char stub[MAXCHR], infile[MAXCHR];
    char prbf[MAXCHR], outfil[MAXCHR];
    char netf[MAXCHR], trfil[MAXCHR];
    char cnum[10];

    struct tot_stat *stt=NULL;
    struct tot_stat *prb=NULL;
    struct tot_stat *sout1=NULL, *sout2=NULL;
    struct pt_stat *stmp=NULL, *stmpt=NULL, *stmpo=NULL, *stmpot=NULL;

    GRID gstat;

    FILE *fin=NULL, *fout=NULL, *fprbf=NULL;
    
    gstat.xgrid = NULL;

    printf("Input file stub for input files.\n\n");
    scanf("%s", stub);

    printf("How many samples are there to collect statistics on?\n\n");
    scanf("%d", &nsamp);

    nnsmp = nsamp - 1;
    nns = (float) nsamp;
    
    printf("Are samples diff files, '0' for no and '1' for yes.\n\n");
    scanf("%d", &idd);


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

/* assign output statistics data structure */

    sout1 = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((sout1 == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(sout1, prb, sizeof(struct tot_stat));

    sout1->ptnum = prb->ptnum;
    ptnum = prb->ptnum;

    sout1->ptst = (struct pt_stat * )calloc(ptnum, sizeof(struct pt_stat));
    mem_er((sout1->ptst == NULL) ? 0 : 1, ptnum * sizeof(struct pt_stat));

    sout2 = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((sout2 == NULL) ? 0 : 1, sizeof(struct tot_stat));

    memcpy(sout2, prb, sizeof(struct tot_stat));

    sout2->ptnum = prb->ptnum;

    sout2->ptst = (struct pt_stat * )calloc(ptnum, sizeof(struct pt_stat));
    mem_er((sout2->ptst == NULL) ? 0 : 1, ptnum * sizeof(struct pt_stat));

    for(i=0; i < ptnum; i++){
       stmp = prb->ptst + i;
       stmpo = sout1->ptst + i;

       stmpo->ptyp = stmp->ptyp;
       stmpo->xs = stmp->xs;
       stmpo->ys = stmp->ys;
    }

    while(inc <= nnsmp){

       if(inc < 10) sprintf(cnum, "00000%d", inc);
       else if(inc < 100) sprintf(cnum, "0000%d", inc);
       else if(inc < 1000) sprintf(cnum, "000%d", inc);
       else if(inc < 10000) sprintf(cnum, "00%d", inc);
       strcpy(infile, stub);
       strcat(infile, cnum);
       if(idd) strcat(infile, "_diff"); 

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

       if(stt->ptnum != prb->ptnum){
          printf("****ERROR****, grids in stats files do not match.\n\n");
          exit(1);
       }

       for(i=0; i < ptnum; i++){
           stmp = stt->ptst + i;
           stmpt = prb->ptst + i;
           stmpo = sout1->ptst + i;
           stmpot = sout2->ptst + i;
           if((stmp->stat1).mean > (stmpt->stat1).mean) (stmpo->stat1).mean += 1.0;
           else (stmpot->stat1).mean += 1.0;

           if((stmp->stat1).var > (stmpt->stat1).var) (stmpo->stat1).var += 1.0;
           else (stmpot->stat1).var += 1.0;

           if((stmp->stat2).mean > (stmpt->stat2).mean) (stmpo->stat2).mean += 1.0;
           else (stmpot->stat2).mean += 1.0;

           if((stmp->stat2).var > (stmpt->stat2).var) (stmpo->stat2).var += 1.0;
           else (stmpot->stat2).var += 1.0;

           if(stmp->stat3 > stmpt->stat3) stmpo->stat3 += 1.0;
           else stmpot->stat3 += 1.0;

           if(stmp->stat4 > stmpt->stat4) stmpo->stat4 += 1.0;
           else stmpot->stat4 += 1.0;

           if(stmp->stat5 > stmpt->stat5) stmpo->stat5 += 1.0;
           else stmpot->stat5 += 1.0;

           if(stmp->stat6 > stmpt->stat6) stmpo->stat6 += 1.0;
           else stmpot->stat6 += 1.0;

           if(stmp->stat8 > stmpt->stat8) stmpo->stat8 += 1.0;
           else stmpot->stat8 += 1.0;

           if(stmp->stat9 > stmpt->stat9) stmpo->stat9 += 1.0;
           else stmpot->stat9 += 1.0;

           if(stmp->stat10 > stmpt->stat10) stmpo->stat10 += 1.0;
           else stmpot->stat10 += 1.0;

           if(stmp->stat12 > stmpt->stat12) stmpo->stat12 += 1.0;
           else stmpot->stat12 += 1.0;

           if(stmp->stat13 > stmpt->stat13) stmpo->stat13 += 1.0;
           else stmpot->stat13 += 1.0;


       }

       free(stt->ptst);
       free(stt);

       ++inc;

    }

    for(i=0; i < ptnum; i++){
       stmpo = sout1->ptst + i;
       stmpot = sout2->ptst + i;
       (stmpo->stat1).mean = ((stmpo->stat1).mean < (stmpot->stat1).mean) ? (stmpo->stat1).mean : (stmpot->stat1).mean;
       (stmpo->stat1).mean *= 2.0 / nns;
       (stmpo->stat1).var = ((stmpo->stat1).var < (stmpot->stat1).var) ? (stmpo->stat1).var : (stmpot->stat1).var;
       (stmpo->stat1).var *= 2.0 / nns;
       (stmpo->stat2).mean = ((stmpo->stat2).mean < (stmpot->stat2).mean) ? (stmpo->stat2).mean : (stmpot->stat2).mean;
       (stmpo->stat2).mean *= 2.0 / nns;
       (stmpo->stat2).var = ((stmpo->stat2).var < (stmpot->stat2).var) ? (stmpo->stat2).var : (stmpot->stat2).var;
       (stmpo->stat2).var *= 2.0 / nns;
       stmpo->stat3 = (stmpo->stat3 < stmpot->stat3) ? stmpo->stat3 : stmpot->stat3;
       stmpo->stat3 *= 2.0 / nns;
       stmpo->stat4 = (stmpo->stat4 < stmpot->stat4) ? stmpo->stat4 : stmpot->stat4;
       stmpo->stat4 *= 2.0 / nns;
       stmpo->stat5 = (stmpo->stat5 < stmpot->stat5) ? stmpo->stat5 : stmpot->stat5;
       stmpo->stat5 *= 2.0 / nns;
       stmpo->stat6 = (stmpo->stat6 < stmpot->stat6) ? stmpo->stat6 : stmpot->stat6;
       stmpo->stat6 *= 2.0 / nns;
       stmpo->stat8 = (stmpo->stat8 < stmpot->stat8) ? stmpo->stat8 : stmpot->stat8;
       stmpo->stat8 *= 2.0 / nns;
       stmpo->stat9 = (stmpo->stat9 < stmpot->stat9) ? stmpo->stat9 : stmpot->stat9;
       stmpo->stat9 *= 2.0 / nns;
       stmpo->stat10 = (stmpo->stat10 < stmpot->stat10) ? stmpo->stat10 : stmpot->stat10;
       stmpo->stat10 *= 2.0 / nns;
       stmpo->stat12 = (stmpo->stat12 < stmpot->stat12) ? stmpo->stat12 : stmpot->stat12;
       stmpo->stat12 *= 2.0 / nns;
       stmpo->stat13 = (stmpo->stat13 < stmpot->stat13) ? stmpo->stat13 : stmpot->stat13;
       stmpo->stat13 *= 2.0 / nns;
    }


    printf("What output file do you want to write the test information?\n\n");
    scanf("%s", outfil);

    fout = fopen(outfil, "w");
    statdmp(fout, sout1);
    fclose(fout);

    printf("Do you want netcdf output?, 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y'){
       printf("What output netcdf file do you want to write to?\n\n");
       scanf("%s", netf);
       read_gstat(sout1, &gstat, trfil, &snum);
       netcdf_write_stats(sout1, &gstat, netf, trfil, snum, 0);
    }


    if(gstat.xgrid){
       free(gstat.xgrid);
       free(gstat.ygrid);
    }

    free(sout1->ptst);
    free(sout1)

    free(sout2->ptst);
    free(sout2)

    free(prb->ptst);
    free(prb);

    return 0;

}

void read_gstat(struct tot_stat *stat, GRID *gs, char *trfil, int *snum)
{

    int i;

    char gf[MAXCHR];
    char line[MAXCHR];

    struct pt_stat *pst=NULL;

    FILE *fg=NULL;

    printf("What statistics grid file do you want to read?\n\n");
    scanf("%s", gf);

    fg = fopen(gf, "r");
    if(fg == NULL){
       printf("****ERROR****, cant open file \r\n"
              "               %s\r\n"
              "               for read\n\n", gf);
       exit(1);
    }

    fgets(line, MAXCHR, fg);
    sscanf(line, "%s", trfil);
    fgets(line, MAXCHR, fg);
    sscanf(line, "%d", snum);
    fgets(line, MAXCHR, fg);
    sscanf(line, "%d %d", &(gs->ix), &(gs->iy));
    fgets(line, MAXCHR, fg);
    sscanf(line, "%d %d", &(gs->prgr), &(gs->prty));
    fgets(line, MAXCHR, fg);
    sscanf(line, "%f %f", &(gs->alng), &(gs->alat));
    fgets(line, MAXCHR, fg);
    sscanf(line, "%s", gs->prgr_nm);
    fgets(line, MAXCHR, fg);
    sscanf(line, "%s", gs->prty_nm);
    gs->xgrid = (float * )calloc(gs->ix, sizeof(float));
    mem_er((gs->xgrid == NULL) ? 0 :1, gs->ix * sizeof(float));
    gs->ygrid = (float * )calloc(gs->iy, sizeof(float));
    mem_er((gs->ygrid == NULL) ? 0 :1, gs->iy * sizeof(float));
    fread(gs->xgrid, gs->ix * sizeof(float), 1, fg);
    fgets(line, MAXCHR, fg);
    fread(gs->ygrid, gs->iy * sizeof(float), 1, fg);
    fgets(line, MAXCHR, fg);

    pst = stat->ptst;
    for(i=0; i < *snum; i++) fscanf(fg, "%d", &((pst + i)->apos));

    fclose(fg);

    return;
}
