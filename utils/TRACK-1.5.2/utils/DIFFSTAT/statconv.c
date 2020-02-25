#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "statistic.h"
#include "mem_er.h"
#include "grid.h"

#define  MAXCHR   200
#define  TOLFLD   1.0e-6


struct tot_stat *read_stats(FILE * );
void read_gstat(struct tot_stat * , GRID * , char * , char * , int * );
void netcdf_write_stats(struct tot_stat * , GRID * , char * , char * , int , int );

int prgr, prty;

int main(int argc, char **argv)

{

     int snum=0;

     struct tot_stat *stat=NULL;

     GRID gstat;
      
     FILE *statin=NULL;


     char statfo[MAXCHR], statfn[MAXCHR];
     char fgstat[MAXCHR], trfil[MAXCHR];

     if(argc != 4) {
       printf("Usage:- %s [ASCII stats file] [ NETCDF stats file] [gstat file]\n", argv[0]);
       exit(1); 
     }

     sscanf(argv[1], "%s", statfn);
     sscanf(argv[2], "%s", statfo);
     sscanf(argv[3], "%s", fgstat);


     statin = fopen(statfn, "r");
     if(statin == NULL){
        printf("****ERROR****, cant open file \r\n"
               "               %s\r\n"
               "               for read\n\n", statfn);
        exit(1);
     }

     stat = read_stats(statin);

     fclose(statin);

     read_gstat(stat, &gstat, fgstat, trfil, &snum);
     netcdf_write_stats(stat, &gstat, statfo, trfil, snum, 0);

     return 0;

}

void read_gstat(struct tot_stat *stat, GRID *gs, char *fgs, char *trfil, int *snum)
{

    int i;

    char line[MAXCHR];

    struct pt_stat *pst=NULL;

    FILE *fg=NULL;

    fg = fopen(fgs, "r");
    if(fg == NULL){
       printf("****ERROR****, cant open file \r\n"
              "               %s\r\n"
              "               for read\n\n", fgs);
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
