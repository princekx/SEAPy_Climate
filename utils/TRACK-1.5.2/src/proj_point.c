#include <Stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "proj.h"


/* function to take point from one projection and project onto a second */


int proj_point(PROJ proj, float *xpp, float *ypp, GRDP *g1, GRDP *g2)

{
    int ift=0;

    int irt=1;

    switch(g1->prgr){
      case 0:
         if(proj.prj1 != NULL){
            (*proj.prj1)(xpp, 'x', 0);
            (*proj.prj1)(ypp, 'y', 0);
         }
         break;
      case 1:
	 azimuthal(ypp, xpp, g1->alat, g1->alng, 0, g1->prty);
         break;	 
      case 2:   
         irt=conic(ypp, xpp, g1->alat, g1->alng, g1->sp1, g1->sp2, 0, g1->prty, &ift);
	 break;
      case 3:
         irt=rotated_cyln(ypp, xpp, g1->alat, g1->alng, 0, &ift);
	 break;
      default:
         printf("***ERROR***, no such projection group for use in %s\n\n", __FILE__);
         exit(1);
     }

     switch(g2->prgr){
      case 0:
         if(g2->prty > 0){
            (*proj.prj2)(xpp, 'x', 1);
            (*proj.prj2)(ypp, 'y', 1);
         }
         break;
      case 1:
	 irt = azimuthal(ypp, xpp, g2->alat, g2->alng, 1, g2->prty);
         break;
      case 2: 
         irt=conic(ypp, xpp, g2->alat, g2->alng, g2->sp1, g2->sp2, 1, g2->prty, &ift);
	 break;
      case 3:
         irt=rotated_cyln(ypp, xpp, g2->alat, g2->alng, 1, &ift);
	 break;
      default:
         printf("***ERROR***, no such projection group for use in %s\n\n", __FILE__);
         exit(1);
     }

     return irt;

}
