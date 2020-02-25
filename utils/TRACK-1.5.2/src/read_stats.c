#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "statistic.h"
#include "mem_er.h"

/* function to read in data from a file containing statistical data
   from a previous analysis.                                         */

struct tot_stat *read_stats(FILE *fil)

{

    int i;
    int prj, gpr;

    struct tot_stat *trsav;
    struct pt_stat *tr;

    char tex[MXCHR];
    char *stok=NULL;

    fscanf(fil, "%s %d %d\n", tex, &gpr, &prj);

    if(gpr || prj) {

       printf("***ERROR***, incorrect projection specifiers in STATS file read by %s\n\n", __FILE__);
       exit(1);

    }

    trsav = (struct tot_stat * )malloc_initl(sizeof(struct tot_stat));
    mem_er((trsav == NULL) ? 0 : 1, sizeof(struct tot_stat));

    fgets(tex, MXCHR, fil);

    if(sscanf(tex, "%*s %d %*s %d %*s %e", &trsav->ptnum, &trsav->scden, &trsav->add_den_sc) != 3) {trsav->scden=-1; trsav->add_den_sc = 0.0;}

    switch(trsav->scden){
         case 0:
            printf("****INFORMATION****, densities are pdf's except for track density\n\n");
            break;
         case 1:
            printf("****INFORMATION****, densities scaled to traditional number densities.\n\n");
            break;
         case 2:
            printf("****INFORMATION****, all densities are pdf's including track density. \n\n");
            break;
         case 3:
            printf("****INFORMATION****, all densities scaled to number densities from pdf's\n\n");
            break;
         default:
            printf("****INFORMATION****, unknown density types.\n\n");
    }

    trsav->ptst = (struct pt_stat * )calloc(trsav->ptnum, sizeof(struct pt_stat));
    mem_er((trsav->ptst == NULL) ? 0 : 1, trsav->ptnum * sizeof(struct pt_stat));

    fgets(tex, MXCHR, fil);
    sscanf(tex, "%*s %f %f %f %f", &trsav->xa1, &trsav->xa2, &trsav->ya1, &trsav->ya2);

    fgets(tex, MXCHR, fil);

    fgets(tex, MXCHR, fil);
    stok = strtok(tex, " ");

    for(i=0; i < STNM; i++) {
       sscanf(stok, "%d ", &trsav->kern[i]); 

       stok = strtok(NULL, " ");
       if(!stok) break;
    }

    fgets(tex, MXCHR, fil);
    stok = strtok(tex, " ");

    for(i=0; i < STNM; i++) {
        sscanf(stok, "%f ", &trsav->sm[i]);
        stok = strtok(NULL, " ");
        if(!stok) break;
    }

    fgets(tex, MXCHR, fil);
    stok = strtok(tex, " ");

    for(i=0; i < STNM; i++) {
       sscanf(stok, "%f ", &trsav->datnm[i]);

       stok = strtok(NULL, " ");
       if(!stok) break;    
    }

    for(i=0; i < trsav->ptnum; i++) {

        tr = (trsav->ptst) + i;

        fscanf(fil, "%f %f", &tr->xs, &tr->ys);
        fscanf(fil, "%f %f", &(tr->stat1).mean, &(tr->stat1).var);
        fscanf(fil, "%f %f", &(tr->stat2).mean, &(tr->stat2).var);
        fscanf(fil, "%f %f %f %f", &tr->stat3, &tr->stat4, &tr->stat5, &tr->stat6);

        fscanf(fil, "%f %f %f", &(tr->stat7).xcomp, &(tr->stat7).ycomp, &tr->stat8);

        fgets(tex, MXCHR, fil);
        sscanf(tex, "%f %f %f %f %f %f", &tr->stat9, &tr->stat10, &(tr->stat11).xcomp, &(tr->stat11).ycomp, &tr->stat12, &tr->stat13);

    }


    return trsav;

}
