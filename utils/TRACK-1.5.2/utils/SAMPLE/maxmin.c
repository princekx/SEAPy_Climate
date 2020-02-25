#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistic.h"
#include "mem_er.h"

#define  MAXCHR   100
#define  NSTAT    13

struct tot_stat *read_stats(FILE * );

int prgr, prty;

int main(void )
{

    int i,j;
    int nsamp=0;
    int inc=0;
    int ptnum=0;

    char stub[MAXCHR], infile[MAXCHR], outfile[MAXCHR];
    char cnum[10];
    char *statnm[] = {"Mean Intensity", "Intensity STD", "Mean Speed", "Speed STD",
                      "Feature Density", "Genesis Density", "Lysis Density", "Track Density",
                      "Mean Lifetime", "Growth/Decay", "Mean Anisotropy", "Tendency", "Mean Area"};

    struct tot_stat *stt=NULL;

    float *smin[NSTAT], *smax[NSTAT];
    float gsmin[NSTAT], gsmax[NSTAT];
    float val;

    FILE *fin=NULL, *fout=NULL;

    printf("****WARNING****, vector statistics, e.g. velocity are\r\n"
           "                 excluded for the moment.            \n\n");

    printf("Input file stub for input files.\n\n");
    scanf("%s", stub);

    printf("How many samples are there to collect statistics on?\n\n");
    scanf("%d", &nsamp);

    --nsamp;

    
    while(inc <= nsamp){

       if(inc < 10) sprintf(cnum, "00000%d", inc);
       else if(inc < 100) sprintf(cnum, "0000%d", inc);
       else if(inc < 1000) sprintf(cnum, "000%d", inc);
       else if(inc < 10000) sprintf(cnum, "00%d", inc);
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

          for(i=0; i < NSTAT; i++){
              smin[i] = (float *)calloc(stt->ptnum, sizeof(float));
              mem_er((smin[i] == NULL) ? 0 : 1, stt->ptnum * sizeof(float));
              smax[i] = (float *)calloc(stt->ptnum, sizeof(float));
              mem_er((smax[i] == NULL) ? 0 : 1, stt->ptnum * sizeof(float));
          }

          for(i=0; i < stt->ptnum; i++){

              *(smax[0] + i) = *(smin[0] + i) = ((stt->ptst + i)->stat1).mean;
              *(smax[1] + i) = *(smin[1] + i) = ((stt->ptst + i)->stat1).var;
              *(smax[2] + i) = *(smin[2] + i) = ((stt->ptst + i)->stat2).mean;
              *(smax[3] + i) = *(smin[3] + i) = ((stt->ptst + i)->stat2).var;
              *(smax[4] + i) = *(smin[4] + i) = (stt->ptst + i)->stat3;
              *(smax[5] + i) = *(smin[5] + i) = (stt->ptst + i)->stat4;
              *(smax[6] + i) = *(smin[6] + i) = (stt->ptst + i)->stat5;
              *(smax[7] + i) = *(smin[7] + i) = (stt->ptst + i)->stat6;
              *(smax[8] + i) = *(smin[8] + i) = (stt->ptst + i)->stat8;
              *(smax[9] + i) = *(smin[9] + i) = (stt->ptst + i)->stat9;
              *(smax[10] + i) = *(smin[10] + i) = (stt->ptst + i)->stat10;
              *(smax[11] + i) = *(smin[11] + i) = (stt->ptst + i)->stat12;
              *(smax[12] + i) = *(smin[12] + i) = (stt->ptst + i)->stat13;
          }

          ptnum = stt->ptnum;
       }

       else {

          for(i=0; i < stt->ptnum; i++){

              val = ((stt->ptst + i)->stat1).mean;
              if(val < *(smin[0] + i)) *(smin[0] + i) = val;
              else if (val > *(smax[0] + i)) *(smax[0] + i) = val;
              val = ((stt->ptst + i)->stat1).var;
              if(val < *(smin[1] + i)) *(smin[1] + i) = val;
              else if (val > *(smax[1] + i)) *(smax[1] + i) = val;
              val = ((stt->ptst + i)->stat2).mean;
              if(val < *(smin[2] + i)) *(smin[2] + i) = val;
              else if (val > *(smax[2] + i)) *(smax[2] + i) = val;
              val = ((stt->ptst + i)->stat2).var;
              if(val < *(smin[3] + i)) *(smin[3] + i) = val;
              else if (val > *(smax[3] + i)) *(smax[3] + i) = val;
              val = (stt->ptst + i)->stat3;
              if(val < *(smin[4] + i)) *(smin[4] + i) = val;
              else if (val > *(smax[4] + i)) *(smax[4] + i) = val;
              val = (stt->ptst + i)->stat4;
              if(val < *(smin[5] + i)) *(smin[5] + i) = val;
              else if (val > *(smax[5] + i)) *(smax[5] + i) = val;
              val = (stt->ptst + i)->stat5;
              if(val < *(smin[6] + i)) *(smin[6] + i) = val;
              else if (val > *(smax[6] + i)) *(smax[6] + i) = val;
              val = (stt->ptst + i)->stat6;
              if(val < *(smin[7] + i)) *(smin[7] + i) = val;
              else if (val > *(smax[7] + i)) *(smax[7] + i) = val;
              val = (stt->ptst + i)->stat8;
              if(val < *(smin[8] + i)) *(smin[8] + i) = val;
              else if (val > *(smax[8] + i)) *(smax[8] + i) = val;
              val = (stt->ptst + i)->stat9;
              if(val < *(smin[9] + i)) *(smin[9] + i) = val;
              else if (val > *(smax[9] + i)) *(smax[9] + i) = val;
              val = (stt->ptst + i)->stat10;
              if(val < *(smin[10] + i)) *(smin[10] + i) = val;
              else if (val > *(smax[10] + i)) *(smax[10] + i) = val;
              val = (stt->ptst + i)->stat12;
              if(val < *(smin[11] + i)) *(smin[11] + i) = val;
              else if (val > *(smax[11] + i)) *(smax[11] + i) = val;
              val = (stt->ptst + i)->stat13;
              if(val < *(smin[12] + i)) *(smin[12] + i) = val;
              else if (val > *(smax[12] + i)) *(smax[12] + i) = val;
          }  

       }
       
       free(stt->ptst);
       free(stt);

       ++inc;

    }

    printf("What is the ouput file for max/min info.\n\n");
    scanf("%s", outfile);

    fout = fopen(outfile, "w");
    if(fout == NULL){
       printf("****ERROR****, cant open file \r\n"
              "               %s\r\n"
              "               for read\n\n", outfile);
       exit(1);
    }

    fprintf(fout, "%d %d\n", ptnum, NSTAT);

    for(i=0; i < ptnum; i++){
       for(j=0; j < NSTAT; j++) fprintf(fout, "%f %f ", *(smin[j] + i), *(smax[j] + i));
       fprintf(fout, "\n");
    }


    fclose(fout);


    printf("GLOBAL MAX/MIN\n\n");

    for(i=0; i<NSTAT; i++){

        gsmin[i] = *smin[i];
        gsmax[i] = *smax[i];
        for(j=0; j < ptnum; j++){
            if(*(smin[i] + j) < gsmin[i]) gsmin[i] = *(smin[i] + j);
            if(*(smax[i] + j) > gsmax[i]) gsmax[i] = *(smax[i] + j);
        }

        printf("%20s %e %e\n", statnm[i], gsmin[i], gsmax[i]);

    }


    for(i=0; i < NSTAT; i++){
        free(smin[i]);
        free(smax[i]);
    }

    return 0;

}
