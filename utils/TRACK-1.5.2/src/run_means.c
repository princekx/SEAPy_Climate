#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include "statistic.h"
#include "file_handle.h"
#include "file_cat_out.h"
#include "mem_er.h"

#define  TOLFD   1.0e-4
#define  TOLWT   1.0e-10
#define  IWS     19

/* routine to compute the combined statistics of seperate runs
   e.g. seasonal mean over a number of years.                     */

struct tot_stat *read_stats(FILE * );
void error_stats(int , int );
void statdmp(FILE * , struct tot_stat * );

extern char *fext;

extern int iext;

void run_means()

{

   int lpc=0, ns=0;
   int i;
   int isc=0;
   int iavg=0;

   float wt;

   char  statfil[MAXCHR];
   char  soutf[MAXCHR];

   double ws[IWS], mdif;
   double wa, wb, wc, wd;
   double sumwt=0.0;

   FILE *statin=NULL, *statout=NULL;

   struct tot_stat *comstat=NULL, *strun=NULL;
   struct pt_stat *st1=NULL, *st2=NULL;

   strncpy(soutf, STATCOM, MAXCHR);
   if(iext) strcpy(strstr(soutf, EXTENSION), fext);


   printf("How many sets of statistics are to be combined?\n");
   scanf("%d", &ns);


   while(lpc < ns){

      printf("What is the %d input file name or 'q' for quit??\n\n", lpc+1);

      scanf("%s", statfil);

      if(statfil[0] == 'q') break;

/* read in stats data and check compatability */

      statin = open_file(statfil, "r");

      if(!statin) {

         printf("****ERROR****, unable to open file %s for read.\n", statfil);
         break;
      }


      if(!lpc){

        comstat = read_stats(statin);

        if(comstat->scden == 1 || comstat->scden == 3) {
          isc = 1;

          printf("***WARNING***, weighted statistics not available for statistics    \r\n"
                 "               already scaled to number density. Use raw statistics\r\n\n");
          wt = 1.0;

          printf("Do you want accumulated statistics or average statistics, input '0' or '1'\n\n");
          scanf("%d", &iavg);

        }

        else {

           printf("****INFORMATION****, each set of statistics can now have an        \r\n"
                  "                     additional weight for the combination process.\r\n"
                  "                     Set these weights to unity if no weighting is \r\n"
                  "                     required.                                     \n\n");


           printf("What weight is required, choose 1 for no weighting. \n\n");
           scanf("%f", &wt);


           if(fabs(wt) < TOLWT){

              printf("****ERROR****, weight cannot be zero for first set of statistics.\r\n"
                  "               Load the sets in a different order.               \n\n");
              exit(1);
           }

           for(i=0; i<STNM; i++){

               comstat->datnm[i] *= wt;

           }

        }

        sumwt = wt;


      }

      else if(lpc > 0){

        strun = read_stats(statin);

        if(isc){

           wt = 1.0;

        }

        else {

           printf("What weight is required, choose 1 for no weighting. \n\n");
           scanf("%f", &wt);

           for(i=0; i<STNM; i++){

               strun->datnm[i] *= wt;

           }

        }

        sumwt += wt;

        if(comstat->ptnum != strun->ptnum) error_stats(1, __LINE__);

        else if(fabs(comstat->xa1 - strun->xa1) > TOLFD ||
                fabs(comstat->xa2 - strun->xa2) > TOLFD ||
                fabs(comstat->ya1 - strun->ya1) > TOLFD ||
                fabs(comstat->ya2 - strun->ya2) > TOLFD    ) error_stats(2, __LINE__);

        if(comstat->scden != strun->scden) error_stats(6, __LINE__);

        for(i=0; i < STNM; i++){

           if(comstat->kern[i] != strun->kern[i]) error_stats(3, __LINE__);
           if(fabs(comstat->sm[i] - strun->sm[i]) > TOLFD) error_stats(4, __LINE__);

        }

        for(i=0; i < comstat->ptnum; i++){

           st1 = (comstat->ptst) + i;
           st2 = (strun->ptst) + i;

           if(fabs(st1->xs - st2->xs) > TOLFD ||
              fabs(st1->ys - st2->ys) > TOLFD    ) error_stats(5, __LINE__);

        }


/* combine distributions */

        for(i=0; i < comstat->ptnum; i++){

          
            st1 = (comstat->ptst) + i;
            st2 = (strun->ptst) + i;

            if(isc){

/* ws[0] == feature density, sum(n*f) */

               ws[0] = st1->stat3 + st2->stat3;

               wa = ((double)(st1->stat3) / (double)(comstat->datnm[4]));
               wb = ((double)(st2->stat3) / (double)(strun->datnm[4]));

/* ws[1] == spead density, sum(n*f) */

               ws[1] = comstat->datnm[2] * wa + strun->datnm[2] * wb;

/* ws[2] == genesis density, sum(n*f) */

               ws[2] = st1->stat4 + st2->stat4;

/* ws[3] == lysis density, sum(n*f) */

               ws[3] = st1->stat5 + st2->stat5;
      
/* ws[4] == track density, sum(n*f) */

               ws[4] = st1->stat6 + st2->stat6;

               wc = ((double)(st1->stat6) / (double)(comstat->datnm[7]));
               wd = ((double)(st2->stat6) / (double)(strun->datnm[7]));

/* ws[5] == mean strength, sum(n*f*m) */

               ws[5] = comstat->datnm[0] * wa * (st1->stat1).mean +
                       strun->datnm[0] * wb * (st2->stat1).mean;

/* ws[6] == mean spead, sum(n*f*m) */

               ws[6] = comstat->datnm[2] * wa * (st1->stat2).mean +
                       strun->datnm[2] * wb * (st2->stat2).mean;

/* ws[7] == X-compenent of track flux, sum(n*f*m) */

               ws[7] = comstat->datnm[8] * wa * (st1->stat7).xcomp +
                       strun->datnm[8] * wb * (st2->stat7).xcomp;

/* ws[8] == Y-compenent of track flux, sum(n*f*m) */

               ws[8] = comstat->datnm[8] * wa * (st1->stat7).ycomp +
                       strun->datnm[8] * wb * (st2->stat7).ycomp;

/* ws[9] == mean lifetime, sum(n*f*m) */

               ws[9] = comstat->datnm[9] * wc * st1->stat8 + 
                       strun->datnm[9] * wd * st2->stat8;

/* ws[14] == mean growth/decay rate */

               ws[14] = comstat->datnm[10] * wa * st1->stat9 + 
                        strun->datnm[10] * wb * st2->stat9;

/* ws[15] == mean anisotropy */

               ws[15] = comstat->datnm[11] * wa * st1->stat10 + 
                        strun->datnm[11] * wb * st2->stat10;

/* ws[16] == X-compenent of orientation, sum(n*f*m) */

               ws[16] = comstat->datnm[12] * wa * (st1->stat11).xcomp +
                        strun->datnm[12] * wb * (st2->stat11).xcomp;

/* ws[17] == Y-compenent of orientation, sum(n*f*m) */

               ws[17] = comstat->datnm[12] * wa * (st1->stat11).ycomp +
                        strun->datnm[12] * wb * (st2->stat11).ycomp;

/* ws[10] == 1st term of varience of strength, sum(n*f*s*s) */

               ws[10] = comstat->datnm[0] * wa * (st1->stat1).var * (st1->stat1).var +
                        strun->datnm[0] * wb * (st2->stat1).var * (st2->stat1).var;

/* ws[11] == 2nd term of varience of strength, N*F*n*f*(M - m)^2 */

               mdif = (st1->stat1).mean - (st2->stat1).mean;

               ws[11] = comstat->datnm[0] * wa * strun->datnm[0] * wb *
                        mdif * mdif;

/* ws[12] == 1st term of varience of speed, sum(n*f*s*s) */

               ws[12] = comstat->datnm[2] * wa * (st1->stat2).var * (st1->stat2).var +
                        strun->datnm[2] * wb * (st2->stat2).var * (st2->stat2).var;

/* ws[13] == 2nd term of varience of spead, N*F*n*f*(M - m)^2 */

               mdif = (st1->stat2).mean - (st2->stat2).mean;

               ws[13] = comstat->datnm[2] * wa * strun->datnm[2] * wb *
                        mdif * mdif;

/* ws[18] == tendency */

               ws[18] = comstat->datnm[13] * wa * st1->stat12 + 
                        strun->datnm[13] * wb * st2->stat12;


            }

            else{

/* ws[0] == feature density, sum(n*f) */

               ws[0] = comstat->datnm[4] * st1->stat3 + strun->datnm[4] * st2->stat3;

/* ws[1] == spead density, sum(n*f) */

               ws[1] = comstat->datnm[2] * st1->stat3 + strun->datnm[2] * st2->stat3;

/* ws[2] == genesis density, sum(n*f) */

               ws[2] = comstat->datnm[5] * st1->stat4 + strun->datnm[5] * st2->stat4;

/* ws[3] == lysis density, sum(n*f) */

               ws[3] = comstat->datnm[6] * st1->stat5 + strun->datnm[6] * st2->stat5;
      
/* ws[4] == track density, sum(n*f) */

               ws[4] = comstat->datnm[7] * st1->stat6 + strun->datnm[7] * st2->stat6;

/* ws[5] == mean strength, sum(n*f*m) */

               ws[5] = comstat->datnm[0] * st1->stat3 * (st1->stat1).mean +
                       strun->datnm[0] * st2->stat3 * (st2->stat1).mean;

/* ws[6] == mean spead, sum(n*f*m) */

               ws[6] = comstat->datnm[2] * st1->stat3 * (st1->stat2).mean +
                       strun->datnm[2] * st2->stat3 * (st2->stat2).mean;

/* ws[7] == X-compenent of track flux, sum(n*f*m) */

               ws[7] = comstat->datnm[8] * st1->stat3 * (st1->stat7).xcomp +
                       strun->datnm[8] * st2->stat3 * (st2->stat7).xcomp;

/* ws[8] == Y-compenent of track flux, sum(n*f*m) */

               ws[8] = comstat->datnm[8] * st1->stat3 * (st1->stat7).ycomp +
                       strun->datnm[8] * st2->stat3 * (st2->stat7).ycomp;

/* ws[9] == mean lifetime, sum(n*f*m) */

               ws[9] = comstat->datnm[9] * st1->stat6 * st1->stat8 + 
                       strun->datnm[9] * st2->stat6 * st2->stat8;

/* ws[14] == mean growth/decay rate */

               ws[14] = comstat->datnm[10] * st1->stat3 * st1->stat9 + 
                        strun->datnm[10] * st2->stat3 * st2->stat9;

/* ws[15] == mean anisotropy */

               ws[15] = comstat->datnm[11] * st1->stat3 * st1->stat10 + 
                        strun->datnm[11] * st2->stat3 * st2->stat10;

/* ws[16] == X-compenent of orientation, sum(n*f*m) */

               ws[16] = comstat->datnm[12] * st1->stat3 * (st1->stat11).xcomp +
                        strun->datnm[12] * st2->stat3 * (st2->stat11).xcomp;

/* ws[17] == Y-compenent of orientation, sum(n*f*m) */

               ws[17] = comstat->datnm[12] * st1->stat3 * (st1->stat11).ycomp +
                        strun->datnm[12] * st2->stat3 * (st2->stat11).ycomp;

/* ws[10] == 1st term of varience of strength, sum(n*f*s*s) */

               ws[10] = comstat->datnm[0] * st1->stat3 * (st1->stat1).var * (st1->stat1).var +
                        strun->datnm[0] * st2->stat3 * (st2->stat1).var * (st2->stat1).var;

/* ws[11] == 2nd term of varience of strength, N*F*n*f*(M - m)^2 */

               mdif = (st1->stat1).mean - (st2->stat1).mean;

               ws[11] = comstat->datnm[0] * st1->stat3 * strun->datnm[0] * st2->stat3 *
                        mdif * mdif;

/* ws[12] == 1st term of varience of speed, sum(n*f*s*s) */

               ws[12] = comstat->datnm[2] * st1->stat3 * (st1->stat2).var * (st1->stat2).var +
                        strun->datnm[2] * st2->stat3 * (st2->stat2).var * (st2->stat2).var;

/* ws[13] == 2nd term of varience of spead, N*F*n*f*(M - m)^2 */

               mdif = (st1->stat2).mean - (st2->stat2).mean;

               ws[13] = comstat->datnm[2] * st1->stat3 * strun->datnm[2] * st2->stat3 *
                        mdif * mdif;

/* ws[18] == tendency */

               ws[18] = comstat->datnm[13] * st1->stat3 * st1->stat12 + 
                        strun->datnm[13] * st2->stat3 * st2->stat12;

            }

/* new values */


            if(ws[0] > DTOL){
              (st1->stat1).mean = ws[5] / ws[0];
              (st1->stat1).var = sqrt((ws[10] / ws[0]) + (ws[11] / (ws[0] * ws[0])));

              st1->stat10 = ws[15] / ws[0];

              (st1->stat11).xcomp = ws[16] / ws[0];
              (st1->stat11).ycomp = ws[17] / ws[0];

            }
            else {
              (st1->stat1).mean = (st1->stat1).var = 0.0;
              st1->stat10 = 0.0;
              (st1->stat11).xcomp = (st1->stat11).ycomp = 0.0;
            }

            if(ws[1] > DTOL){
              (st1->stat2).mean = ws[6] / ws[1];
              (st1->stat2).var = sqrt((ws[12] / ws[1]) + (ws[13] / (ws[1] * ws[1])));
              (st1->stat7).xcomp = ws[7] / ws[1];
              (st1->stat7).ycomp = ws[8] / ws[1];

              st1->stat9 = ws[14] / ws[1];
              st1->stat12 = ws[18] / ws[1];

            }
            else {
              (st1->stat2).mean = (st1->stat2).var = 0.0;
              (st1->stat7).xcomp = (st1->stat7).ycomp = 0.0;


              st1->stat9 = st1->stat10 = st1->stat12 = 0.0;

            }

            if(isc){

               st1->stat3 = ws[0];
               st1->stat4 = ws[2];
               st1->stat5 = ws[3];
               st1->stat6 = ws[4];
 
               if(ws[4] > DTOL) st1->stat8 = ws[9] / ws[4]; 
               else st1->stat8 = 0.0;

            }

            else{

               if(comstat->datnm[4] + strun->datnm[4] > 0.0)
                  st1->stat3 = ws[0] / (double)(comstat->datnm[4] + strun->datnm[4]);
               if(comstat->datnm[5] + strun->datnm[5] > 0.0)
                  st1->stat4 = ws[2] / (double)(comstat->datnm[5] + strun->datnm[5]);
               if(comstat->datnm[6] + strun->datnm[6] > 0.0)
                  st1->stat5 = ws[3] / (double)(comstat->datnm[6] + strun->datnm[6]);

               if(comstat->datnm[7] + strun->datnm[7] > 0)
                  st1->stat6 = ws[4] / (double)(comstat->datnm[7] + strun->datnm[7]);
 
               if(ws[4] > DTOL) st1->stat8 = ws[9] / ws[4]; 
               else st1->stat8 = 0.0;

            }

        }

        for(i=0; i < STNM; i++) comstat->datnm[i] += strun->datnm[i];

        free(strun->ptst);
        free(strun);

      }

      ++lpc;


   }

   if(statin) close_file(statin, statfil);


   if(isc && iavg) {

     for(i=0; i < STNM; i++) comstat->datnm[i] /= lpc;

     for(i=0; i < comstat->ptnum; i++){
        st1 = (comstat->ptst) + i;
        st1->stat3 /= sumwt;
        st1->stat4 /= sumwt;
        st1->stat5 /= sumwt;
        st1->stat6 /= sumwt;
     }

   }

   statout = open_file(soutf, "w");
   if(comstat) statdmp(statout, comstat);
   close_file(statout, soutf);

   printf("***INFORMATION***, the statistics written to file %s            \r\n"
          "                   can be displayed by choosing the options:    \r\n"
          "                   'display only existing analyses' and         \r\n"
          "                   'propog. speeds, and statistics'.            \n\n", soutf);

   return;

}
