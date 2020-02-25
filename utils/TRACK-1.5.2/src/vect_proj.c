#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "mem_er.h"
#include "proj.h"
#include "grid.h"

#define  VECFRAC   0.1
#define  TOLVEC    1.0e-5

/* function to perform transformation of vectors to non-cylyndrical 
   projections.                                                        */


void vect_proj(float xss, float yss, float *xcmp, float *ycmp, PROJ proj, GRDP *g1, GRDP *g2)

{

    float x1, y1, x2, y2, xc, yc;
 
    double vlen1, vlen2, vrat;


    if(!g2->prgr){
       printf("***ERROR***, incorrect projection group type\n\n");
       exit(1);
    }

    xc = *xcmp;
    yc = *ycmp;
    vlen1 = sqrt(xc*xc + yc*yc);

    if(vlen1 <= TOLVEC) return;

    x1 = xss;
    y1 = yss;
    x2 = x1 + (xc * VECFRAC / vlen1);
    y2 = y1 + (yc * VECFRAC / vlen1);

    proj_point(proj, &x1, &y1, g1, g2);
    if(!proj_point(proj, &x2, &y2, g1, g2)) {x2 = x1; y2 = y1;}

    *xcmp = x2 - x1;
    *ycmp = y2 - y1;

    vlen2 = sqrt(*xcmp * *xcmp + *ycmp * *ycmp);

    if(vlen2 > 0.) vrat = vlen1 / vlen2;
    else vrat = 0.;

    *xcmp *= vrat;
    *ycmp *= vrat;

       
    return;

}
