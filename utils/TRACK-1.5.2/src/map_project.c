#include <Stdio.h>
#include <string.h>
#include <proj.h>
#include "grid.h"

/* perform a change of map projection */

extern int init;
extern GRID *gr;

extern FILE *finit;

PROJ map_project(int *pt, int cpy)

{

   int ptt;

   static char *prname[] = {"Plate Caree", "Equirectangular",
                            "Equal Area", "Urmaev", "Mercator",
                            "Stereographic", "Miller"};

   PROJ proj={NULL, NULL};
   PRFP prj1, prj2;

   proj_report(*pt, 0);

   ptt = *pt;

   prj1 = proj_assign(ptt);
   proj.prj1 = prj1;  

   if(cpy < 0){

     printf("The available cylindrical projection types are: \r\n"
          "Plate Caree             input '0'               \r\n"
          "Equirectangular,        input '1'               \r\n"
          "Equal Area,             input '2'               \r\n"
          "Urmaev,                 input '3'               \r\n"
          "Mercator,               input '4'               \r\n"
          "Stereographic,          input '5'               \r\n"
          "Miller,                 input '6'               \r\n");


     if(init) fscanf(finit, "%d", pt);
     else{
        scanf("%d", pt);
        if(finit)fprintf(finit, "%d\n", *pt);
     }

   }

   else *pt = cpy;

   strncpy(gr->prty_nm, prname[*pt], MXPRCH);

   prj2 = proj_assign(*pt);

/* transform to geo-coordinates first */

   if(ptt > 0) {

     if(prj2 == NULL) {*pt = 0; prj2 = prj1;}

     trans_proj(prj1, 0);

   }

/* perform required transformation */

   if(*pt != 0) trans_proj(prj2, 1);

   proj.prj2 = prj2;

   proj_report(*pt, 0);

   return proj;

}
