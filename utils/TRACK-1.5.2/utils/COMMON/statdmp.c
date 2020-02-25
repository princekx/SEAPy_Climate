#include <stdio.h>
#include <stdlib.h>
#include "statistic.h"

/* function to write statistical data to file */

extern int prgr, prty;


void statdmp(FILE *fil, struct tot_stat *st)

{

    int i;

    struct pt_stat *pt;

    switch(st->scden){
         case 0:
            printf("****INFORMATION****, writing densities as pdf's except for track density\n\n");
            break;
         case 1:
            printf("****INFORMATION****, writing densities scaled to traditional number densities.\n\n");
            break;
         case 2:
            printf("****INFORMATION****, writing all densities as pdf's including track density. \n\n");
            break;
         case 3:
            printf("****INFORMATION****, writing all densities scaled to number densities from pdf's\n\n");
            break;
         default:
            printf("****INFORMATION****, unknown density types.\n\n");
    }

    fprintf(fil, "PROJ_ID  %d %d\n", prgr, prty);

    fprintf(fil, "POINT_NUM  %d SCALED_DENSITY %d ADDITIONAL_DENSITY_SCALING %.8f\n", st->ptnum, st->scden, st->add_den_sc);


    fprintf(fil, "DOMAIN %f %f %f %f\n", st->xa1, st->xa2, st->ya1, st->ya2);

    fprintf(fil, "SMOOTHING\n");

    for(i=0; i < STNM; i++) fprintf(fil, "%d ", st->kern[i]);
    fprintf(fil, "\n");
    for(i=0; i < STNM; i++) fprintf(fil, "%f ", st->sm[i]);
    fprintf(fil, "\n");
    for(i=0; i < STNM; i++) fprintf(fil, "%f ", st->datnm[i]);
    fprintf(fil, "\n");

    for(i=0; i < st->ptnum; i++){

       pt = (st->ptst) + i;

       fprintf(fil, "%f  %f  ", pt->xs, pt->ys);
       fprintf(fil, "%f %f ", (pt->stat1).mean, (pt->stat1).var);
       fprintf(fil, "%f %f ", (pt->stat2).mean, (pt->stat2).var);
       fprintf(fil, "%f %f %f %f ", pt->stat3, pt->stat4, pt->stat5, pt->stat6);
       fprintf(fil, "%f %f %f %f ", (pt->stat7).xcomp, (pt->stat7).ycomp, pt->stat8, pt->stat9);
       fprintf(fil, "%f %f %f ", pt->stat10, (pt->stat11).xcomp, (pt->stat11).ycomp);
       fprintf(fil, "%f %f\n", pt->stat12, pt->stat13);

    }

    return;

}
