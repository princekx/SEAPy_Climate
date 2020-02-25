#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "statistic.h"
#include "mem_er.h"

#define  MAXCHR   1000
#define  TOLFLD   1.0e-6


struct tot_stat *read_stats(FILE * );
void statdmp(FILE * , struct tot_stat * );
float samesign(float , float );
float extreme(float , float );

int prgr, prty;


int main(int argc, char **argv)

{

     int i, j;
     int pg, pt;
     int ns=0;
     int ipt=0;

     float pth=0.0;

     struct tot_stat *base=NULL, *stat=NULL, *prst=NULL;
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

     prst = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
     mem_er((prst == NULL) ? 0 : 1, sizeof(struct tot_stat));

     prst->ptnum = base->ptnum;

     prst->xa1 = base->xa1;
     prst->xa2 = base->xa2;
     prst->ya1 = base->ya1;
     prst->ya2 = base->ya2;

     for(i=0; i < STNM; i++){

         prst->kern[i] = base->kern[i];
         prst->sm[i] = base->sm[i];
         prst->scden = base->scden;
         prst->add_den_sc = base->add_den_sc;

     }

     prst->ptst = (struct pt_stat * )calloc(prst->ptnum, sizeof(struct pt_stat));
     mem_er((prst->ptst == NULL) ? 0 : 1, prst->ptnum * sizeof(struct pt_stat));

     for(i=0; i < prst->ptnum; i++){
         pts = prst->ptst + i;
         pts->xs = (base->ptst + i)->xs;
         pts->ys = (base->ptst + i)->ys;
        (pts->stat1).mean = 0.0;
        (pts->stat1).var = 0.0;
        (pts->stat2).mean = 0.0;
        (pts->stat2).var = 0.0;
         pts->stat3 = 0.0;
         pts->stat4 = 0.0;
         pts->stat5 = 0.0;
         pts->stat6 = 0.0;
        (pts->stat7).xcomp = 0.0; 
        (pts->stat7).ycomp = 0.0;
         pts->stat8 = 0.0;
         pts->stat9 = 0.0;
         pts->stat10 = 0.0;
        (pts->stat11).xcomp = 0.0;
        (pts->stat11).ycomp = 0.0;
         pts->stat12 = 0.0;
         pts->stat13 = 0.0;
     }

     printf("For the persistence statistic, do you want:     \r\n"
            "   '0' number with same sign,                   \r\n"
            "   '1' number with more extreme value.          \n\n");
     scanf("%d", &ipt); 

     if(ipt < 0 || ipt > 1){
        printf("****ERROR****, incorrect input.\n\n");
        exit(1);
     }

     printf("How many other statistics file should be compared with the base statistics.\n\n");
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

        if(!ipt){
           for(j=0; j < base->ptnum; j++){ 

              pts = prst->ptst + j;
              ptsb = base->ptst + j;
              ptss = stat->ptst + j;
   
             (pts->stat1).mean += samesign((ptsb->stat1).mean, (ptss->stat1).mean);
             (pts->stat1).var += samesign((ptsb->stat1).var, (ptss->stat1).var);
             (pts->stat2).mean += samesign((ptsb->stat2).mean, (ptss->stat2).mean);
             (pts->stat2).var += samesign((ptsb->stat2).var, (ptss->stat2).var);
              pts->stat3 += samesign(ptsb->stat3, ptss->stat3);
              pts->stat4 += samesign(ptsb->stat4, ptss->stat4);
              pts->stat5 += samesign(ptsb->stat5, ptss->stat5);
              pts->stat6 += samesign(ptsb->stat6, ptss->stat6);
             (pts->stat7).xcomp += samesign((ptsb->stat7).xcomp, (ptss->stat7).xcomp);
             (pts->stat7).ycomp += samesign((ptsb->stat7).ycomp, (ptss->stat7).ycomp);
              pts->stat8 += samesign(ptsb->stat8, ptss->stat8);
              pts->stat9 += samesign(ptsb->stat9, ptss->stat9);
              pts->stat10 += samesign(ptsb->stat10, ptss->stat10);
             (pts->stat11).xcomp += samesign((ptsb->stat11).xcomp, (ptss->stat11).xcomp);
             (pts->stat11).ycomp += samesign((ptsb->stat11).ycomp, (ptss->stat11).ycomp);              
              pts->stat12 += samesign(ptsb->stat12, ptss->stat12);
              pts->stat13 += samesign(ptsb->stat13, ptss->stat13);

           }
        }
        else if(ipt == 1){
           for(j=0; j < base->ptnum; j++){
              pts = prst->ptst + j;
              ptsb = base->ptst + j;
              ptss = stat->ptst + j;
             (pts->stat1).mean += extreme((ptsb->stat1).mean, (ptss->stat1).mean);
             (pts->stat1).var += extreme((ptsb->stat1).var, (ptss->stat1).var);
             (pts->stat2).mean += extreme((ptsb->stat2).mean, (ptss->stat2).mean);
             (pts->stat2).var += extreme((ptsb->stat2).var, (ptss->stat2).var);
              pts->stat3 += extreme(ptsb->stat3, ptss->stat3);
              pts->stat4 += extreme(ptsb->stat4, ptss->stat4);
              pts->stat5 += extreme(ptsb->stat5, ptss->stat5);
              pts->stat6 += extreme(ptsb->stat6, ptss->stat6);
             (pts->stat7).xcomp += extreme((ptsb->stat7).xcomp, (ptss->stat7).xcomp);
             (pts->stat7).ycomp += extreme((ptsb->stat7).ycomp, (ptss->stat7).ycomp);
              pts->stat8 += extreme(ptsb->stat8, ptss->stat8);
              pts->stat9 += extreme(ptsb->stat9, ptss->stat9);
              pts->stat10 += extreme(ptsb->stat10, ptss->stat10);
             (pts->stat11).xcomp += extreme((ptsb->stat11).xcomp, (ptss->stat11).xcomp);
             (pts->stat11).ycomp += extreme((ptsb->stat11).ycomp, (ptss->stat11).ycomp);
              pts->stat12 += extreme(ptsb->stat12, ptss->stat12);
              pts->stat13 += extreme(ptsb->stat13, ptss->stat13);

           }
        }
        else {
          printf("****ERROR****, unknown option.\n\n");
          exit(1);
        }

        free(stat->ptst); 
        free(stat);

     }

/* normalise */

     printf("What threshold is required?\n\n");
     scanf("%f", &pth);

     for(i=0; i <= base->ptnum; i++){
         pts = prst->ptst + i; 
        (pts->stat1).mean /= ns;
         if((pts->stat1).mean < pth) (pts->stat1).mean = 0.0;
        (pts->stat1).var /= ns;
         if((pts->stat1).var < pth) (pts->stat1).var = 0.0;
        (pts->stat2).mean /= ns;
         if((pts->stat2).mean < pth) (pts->stat2).mean = 0.0;
        (pts->stat2).var /= ns;
         if((pts->stat2).var < pth) (pts->stat2).var = 0.0;
         pts->stat3 /= ns;
         if(pts->stat3 < pth) pts->stat3 = 0.0;
         pts->stat4 /= ns;
         if(pts->stat4 < pth) pts->stat4 = 0.0;
         pts->stat5 /= ns;
         if(pts->stat5 < pth) pts->stat5 = 0.0;
         pts->stat6 /= ns;
         if(pts->stat6 < pth) pts->stat6 = 0.0;
        (pts->stat7).xcomp /= ns;
         if((pts->stat7).xcomp < pth) (pts->stat7).xcomp = 0.0;
        (pts->stat7).ycomp /= ns;
         if((pts->stat7).ycomp < pth) (pts->stat7).ycomp = 0.0;
         pts->stat8 /= ns;
         if(pts->stat8 < pth) pts->stat8 = 0.0;
         pts->stat9 /= ns;
         if(pts->stat9 < pth) pts->stat9 = 0.0;
         pts->stat10 /= ns;
         if(pts->stat10 < pth) pts->stat10 = 0.0;
        (pts->stat11).xcomp /= ns;
         if((pts->stat11).xcomp < pth) (pts->stat11).xcomp = 0.0;
        (pts->stat11).ycomp /= ns;
         if((pts->stat11).ycomp < pth) (pts->stat11).ycomp = 0.0;
         pts->stat12 /= ns;
         if(pts->stat12 < pth) pts->stat12 = 0.0;
         pts->stat13 /= ns;
         if(pts->stat13 < pth) pts->stat13 = 0.0;
     }

     printf("Input output filename.\n");
     scanf("%s", prstf);

     stdo = fopen(prstf, "w");

     if(stdo == NULL){
        printf("****ERROR****, cant open file \r\n"
               "               %s\r\n"
               "               for write\n\n", prstf);
        exit(1);
     }

     statdmp(stdo, prst);

     fclose(stdo);

     return 0;

}

float samesign(float sb, float ss)
{
    int ib=0, is=0;

    float fss=0.0;

    ib = (sb > 0.0) - (sb < 0.0);
    is = (ss > 0.0) - (ss < 0.0);
    if(ib && (ib == is)) fss = 1.0; 

    return fss;
}

float extreme(float sb, float ss)
{
    int ib=0;

    float fss=0.0;

    ib = (sb > 0.0) - (sb < 0.0);
    if(ib < 0 && ss <= sb) fss = 1.0;
    else if(ib > 0 && ss >= sb) fss = 1.0;

    return fss;
}
