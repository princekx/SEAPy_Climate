#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "statistic.h"
#include "mem_er.h"

#define  MAXCHR   1000
#define  TOLFLD   1.0e-6


struct tot_stat *read_stats(FILE * );
void statdmp(FILE * , struct tot_stat * );

int prgr, prty;


int main(int argc, char **argv)

{

     int i, j;
     int pg, pt;
     int ns=0;
     int ipt=0;

     float pth=0.0;

     struct tot_stat *base=NULL, *stmn=NULL, *stsd=NULL, *stat=NULL;
     struct pt_stat *pts=NULL, *ptsb=NULL, *ptss=NULL;

      
     FILE *sbase=NULL, *stdo=NULL, *statin=NULL;

     char fbase[MAXCHR], prstf[MAXCHR], statf[MAXCHR];

     printf("What is the base statistic file to read?\n");
     scanf("%s", fbase);

     sbase = fopen(fbase, "r");
     if(sbase == NULL){
        printf("****ERROR****, cant open file \r\n"
               "               %s\r\n"
               "               for read\n\n", fbase);
        exit(1);
     }

     base = read_stats(sbase);

     pg = prgr;
     pt = prty;

     fclose(sbase);

     stmn = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
     mem_er((stmn == NULL) ? 0 : 1, sizeof(struct tot_stat));

     stsd = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
     mem_er((stsd == NULL) ? 0 : 1, sizeof(struct tot_stat)); 

     stsd->ptnum = stmn->ptnum = base->ptnum;

     stsd->xa1 = stmn->xa1 = base->xa1;
     stsd->xa2 = stmn->xa2 = base->xa2;
     stsd->ya1 = stmn->ya1 = base->ya1;
     stsd->ya2 = stmn->ya2 = base->ya2;

     for(i=0; i < STNM; i++){

         stsd->kern[i] = stmn->kern[i] = base->kern[i];
         stsd->sm[i] = stmn->sm[i] = base->sm[i];
         stsd->scden = stmn->scden = base->scden;
         stsd->add_den_sc = stmn->add_den_sc = base->add_den_sc;

     }

     stmn->ptst = (struct pt_stat * )calloc(stmn->ptnum, sizeof(struct pt_stat));
     mem_er((stmn->ptst == NULL) ? 0 : 1, stmn->ptnum * sizeof(struct pt_stat));

     stsd->ptst = (struct pt_stat * )calloc(stsd->ptnum, sizeof(struct pt_stat));
     mem_er((stsd->ptst == NULL) ? 0 : 1, stsd->ptnum * sizeof(struct pt_stat));

     for(i=0; i < base->ptnum; i++){
         pts = stmn->ptst + i;
         ptss = stsd->ptst + i;
         pts->xs = (base->ptst + i)->xs;
         pts->ys = (base->ptst + i)->ys;
         ptss->xs = (base->ptst + i)->xs;
         ptss->ys = (base->ptst + i)->ys;
        (ptss->stat1).mean = (pts->stat1).mean = 0.0;
        (ptss->stat1).var = (pts->stat1).var = 0.0;
        (ptss->stat2).mean = (pts->stat2).mean = 0.0;
        (ptss->stat2).var = (pts->stat2).var = 0.0;
         ptss->stat3 = pts->stat3 = 0.0;
         ptss->stat4 = pts->stat4 = 0.0;
         ptss->stat5 = pts->stat5 = 0.0;
         ptss->stat6 = pts->stat6 = 0.0;
        (ptss->stat7).xcomp = (pts->stat7).xcomp = 0.0; 
        (ptss->stat7).ycomp = (pts->stat7).ycomp = 0.0;
         ptss->stat8 = pts->stat8 = 0.0;
         ptss->stat9 = pts->stat9 = 0.0;
         ptss->stat10 = pts->stat10 = 0.0;
        (ptss->stat11).xcomp = (pts->stat11).xcomp = 0.0;
        (ptss->stat11).ycomp = (pts->stat11).ycomp = 0.0;
         ptss->stat12 = pts->stat12 = 0.0;
         ptss->stat13 = pts->stat13 = 0.0;
     }


     printf("How many statistics files should be used.\n\n");
     scanf("%d", &ns);

     for(i=0; i < ns; i++){
        printf("Input next statistic filename.\n");
        scanf("%s", statf);

        statin = fopen(statf, "r");

        if(statin == NULL){
           printf("****ERROR****, cant open file \r\n"
                  "               %s\r\n"
                  "               for read\n\n", statf);
           exit(1);
        }

         
        stat = read_stats(statin);

        if(pg != prgr || pt != prty || stat->ptnum != base->ptnum ||
           stat->scden != base->scden || stat->add_den_sc != base->add_den_sc){

           printf("****WARNING****, stats files maybe incompatable.\n");


        }

        fclose(statin);
        statin = NULL;

        for(j=0; j < base->ptnum; j++){ 

           pts = stmn->ptst + j;
           ptss = stsd->ptst + j;
           ptsb = stat->ptst + j;
   
          (pts->stat1).mean += (ptsb->stat1).mean;
          (ptss->stat1).mean += (ptsb->stat1).mean * (ptsb->stat1).mean;
          (pts->stat1).var += (ptsb->stat1).var;
          (pts->stat2).mean += (ptsb->stat2).mean;
          (ptss->stat2).mean += (ptsb->stat2).mean * (ptsb->stat2).mean;
          (pts->stat2).var += (ptsb->stat2).var;
           pts->stat3 += ptsb->stat3;
           ptss->stat3 += ptsb->stat3 * ptsb->stat3;
           pts->stat4 += ptsb->stat4;
           ptss->stat4 += ptsb->stat4 * ptsb->stat4;
           pts->stat5 += ptsb->stat5;
           ptss->stat5 += ptsb->stat5 * ptsb->stat5;
           pts->stat6 += ptsb->stat6;
           ptss->stat6 += ptsb->stat6 * ptsb->stat6;
          (pts->stat7).xcomp += (ptsb->stat7).xcomp;
          (ptss->stat7).xcomp += (ptsb->stat7).xcomp * (ptsb->stat7).xcomp;
          (pts->stat7).ycomp += (ptsb->stat7).ycomp;
          (ptss->stat7).ycomp += (ptsb->stat7).ycomp * (ptsb->stat7).ycomp;
           pts->stat8 += ptsb->stat8;
           ptss->stat8 += ptsb->stat8 * ptsb->stat8; 
           pts->stat9 += ptsb->stat9;
           ptss->stat9 += ptsb->stat9 * ptsb->stat9;
           pts->stat10 += ptsb->stat10;
           ptss->stat10 += ptsb->stat10 * ptsb->stat10;
          (pts->stat11).xcomp += (ptsb->stat11).xcomp;
          (ptss->stat11).xcomp += (ptsb->stat11).xcomp * (ptsb->stat11).xcomp;
          (pts->stat11).ycomp += (ptsb->stat11).ycomp;              
          (ptss->stat11).ycomp += (ptsb->stat11).ycomp * (ptsb->stat11).ycomp;
           pts->stat12 += ptsb->stat12;
           ptss->stat12 += ptsb->stat12 * ptsb->stat12;
           pts->stat13 += ptsb->stat13;
           ptss->stat13 += ptsb->stat13 * ptsb->stat13;

        }

        free(stat->ptst); 
        free(stat);

     }

/* normalise */

     for(i=0; i <= base->ptnum; i++){
         pts = stmn->ptst + i; 
         ptss = stsd->ptst + i;
        (pts->stat1).mean /= ns;
        (ptss->stat1).mean /= ns;
        (ptss->stat1).mean = (ptss->stat1).mean - (pts->stat1).mean * (pts->stat1).mean;
        if((ptss->stat1).mean > TOLFLD) (ptss->stat1).mean = sqrt((ptss->stat1).mean);
        else (ptss->stat1).mean = 0.0;
        (pts->stat1).var /= ns;
        (pts->stat2).mean /= ns;
        (ptss->stat2).mean /= ns;
        (ptss->stat2).mean = (ptss->stat2).mean - (pts->stat2).mean * (pts->stat2).mean;
        if((ptss->stat2).mean > TOLFLD) (ptss->stat2).mean = sqrt((ptss->stat2).mean);
        else (ptss->stat2).mean = 0.0;
        (pts->stat2).var /= ns;
         pts->stat3 /= ns;
         ptss->stat3 /= ns;
         ptss->stat3 = ptss->stat3 - pts->stat3 * pts->stat3;
         if(ptss->stat3 > TOLFLD) ptss->stat3 = sqrt(ptss->stat3);
         else ptss->stat3 = 0.0;
         pts->stat4 /= ns;
         ptss->stat4 /= ns;
         ptss->stat4 = ptss->stat4 - pts->stat4 * pts->stat4;
         if(ptss->stat4 > TOLFLD) ptss->stat4 = sqrt(ptss->stat4);
         else ptss->stat4 = 0.0;
         pts->stat5 /= ns;
         ptss->stat5 /= ns;
         ptss->stat5 = ptss->stat5 - pts->stat5 * pts->stat5;
         if(ptss->stat5 > TOLFLD) ptss->stat5 = sqrt(ptss->stat5);
         else ptss->stat5 = 0.0;
         pts->stat6 /= ns;
         ptss->stat6 /= ns; 
         ptss->stat6 = ptss->stat6 - pts->stat6 * pts->stat6;
         if(ptss->stat6 > TOLFLD) ptss->stat6 = sqrt(ptss->stat6);
         else ptss->stat6 = 0.0;
        (pts->stat7).xcomp /= ns;
        (ptss->stat7).xcomp /= ns;
        (ptss->stat7).xcomp = (ptss->stat7).xcomp - (pts->stat7).xcomp * (pts->stat7).xcomp;
         if((ptss->stat7).xcomp > TOLFLD) (ptss->stat7).xcomp = sqrt((ptss->stat7).xcomp);
         else (ptss->stat7).xcomp = 0.0;
        (pts->stat7).ycomp /= ns;
        (ptss->stat7).ycomp /= ns;
        (ptss->stat7).ycomp = (ptss->stat7).ycomp - (pts->stat7).ycomp * (pts->stat7).ycomp;
         if((ptss->stat7).ycomp > TOLFLD) (ptss->stat7).ycomp = sqrt((ptss->stat7).ycomp);
         else (ptss->stat7).ycomp = 0.0;
         pts->stat8 /= ns;
         ptss->stat8 /= ns;
         ptss->stat8 = ptss->stat8 - pts->stat8 * pts->stat8;
         if(ptss->stat8 > TOLFLD) ptss->stat8 = sqrt(ptss->stat8);
         else ptss->stat8 = 0.0;
         pts->stat9 /= ns;
         ptss->stat9 /= ns;
         ptss->stat9 = ptss->stat9 - pts->stat9 * pts->stat9;
         if(ptss->stat9 > TOLFLD) ptss->stat9 = sqrt(ptss->stat9);
         else ptss->stat9 = 0.0;
         pts->stat10 /= ns;
         ptss->stat10 /= ns;
         ptss->stat10 = ptss->stat10 - pts->stat10 * pts->stat10;
         if(ptss->stat10 > TOLFLD) ptss->stat10 = sqrt(ptss->stat10);
         else ptss->stat10 = 0.0;
        (pts->stat11).xcomp /= ns;
        (ptss->stat11).xcomp /= ns;
        (ptss->stat11).xcomp = (ptss->stat11).xcomp - (pts->stat11).xcomp * (pts->stat11).xcomp;
         if((ptss->stat11).xcomp > TOLFLD) (ptss->stat11).xcomp = sqrt((ptss->stat11).xcomp);
         else (ptss->stat11).xcomp= 0.0;
        (pts->stat11).ycomp /= ns;
        (ptss->stat11).ycomp /= ns;
        (ptss->stat11).ycomp = (ptss->stat11).ycomp - (pts->stat11).ycomp * (pts->stat11).ycomp;
         if((ptss->stat11).ycomp > TOLFLD) (ptss->stat11).ycomp = sqrt((ptss->stat11).ycomp);
         else (ptss->stat11).ycomp = 0.0;
         pts->stat12 /= ns;
         ptss->stat12 /= ns;
         ptss->stat12 = ptss->stat12 - pts->stat12 * pts->stat12;
         if(ptss->stat12 > TOLFLD) ptss->stat12 = sqrt(ptss->stat12);
         else ptss->stat12 = 0.0;
         pts->stat13 /= ns;
         ptss->stat13 /= ns;
         ptss->stat13 = ptss->stat13 - pts->stat13 * pts->stat13;
         if(ptss->stat13 > TOLFLD) ptss->stat13 = sqrt(ptss->stat13);
         else ptss->stat13 = 0.0;
     }

     printf("Input output filename for mean.\n");
     scanf("%s", prstf);

     stdo = fopen(prstf, "w");

     if(stdo == NULL){
        printf("****ERROR****, cant open file \r\n"
               "               %s\r\n"
               "               for write\n\n", prstf);
        exit(1);
     }

     statdmp(stdo, stmn);

     fclose(stdo);

     printf("Input output filename for standard deviation.\n");
     scanf("%s", prstf);

     stdo = fopen(prstf, "w");

     if(stdo == NULL){
        printf("****ERROR****, cant open file \r\n"
               "               %s\r\n"
               "               for write\n\n", prstf);
        exit(1);
     }

     statdmp(stdo, stsd);

     fclose(stdo);

     return 0;

}

