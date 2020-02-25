#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grid.h"
#include "proj.h"
#include "mem_er.h"


/* function to set the projection group type */


extern int init;

extern FILE *finit;

extern GRID *gr;

extern GRID *gr1, *gr2;
extern CNTRY *cm1, *cm2;

PROJ proj_group(int pgr, int cyfb)

{

     int pgold;

     PROJ proj={NULL, NULL};

     pgold = gr->prgr;


     if(pgr < 0){

        printf("what kind of projection would you like,                \r\n"
               "Cylindrical            input '0', Mercator etc.        \r\n"
               "Azimuthal              input '1', Stereographic etc    \r\n"
	       "Conic                  input '2', Lambert etc          \r\n"
	       "Rotated Cylindrical    input '3', Mercator etc         \n\n");

        if(init){
           fscanf(finit, "%d", &gr->prgr);
        }
        else{
           scanf("%d", &gr->prgr);

           if(finit) fprintf(finit, "%d\n", gr->prgr);

        }

     }

     else gr->prgr = pgr;

     switch(gr->prgr){
         case 0:
           if(pgold)gr->prty = 0;
           strncpy(gr->prgr_nm, "Cylindrical", MXPRCH);
           proj = map_project(&gr->prty, cyfb);
           break;
         case 1:
           if(!pgold)map_project(&gr->prty, 0);
           strncpy(gr->prgr_nm, "Azimuthal", MXPRCH);
           hemi_grid(cyfb);
           break;
         case 2:
           if(!pgold)map_project(&gr->prty, 0);
           strncpy(gr->prgr_nm, "Conic", MXPRCH);
           conic_grid(cyfb);
           break;
	 case 3:
	   if(!pgold)map_project(&gr->prty, 0);
           strncpy(gr->prgr_nm, "Rotated Cylindrical", MXPRCH);
           rotated_cyln_grid(cyfb);
	   break;	 
         default:
           printf("***ERROR***, projection group type not known\n");
           exit(1);
     }

     if(pgold && !(gr->prgr)){

       if(gr1){free(gr1->xgrid); free(gr1->ygrid); free(gr1);}
       if(cm1){free(cm1->cmxg); free(cm1->cmyg); free(cm1->cmi); free(cm1);}
       if(gr2){free(gr2->xgrid); free(gr2->ygrid); free(gr2);}
       if(cm2){free(cm2->cmxg); free(cm2->cmyg); free(cm2->cmi); free(cm2);}

     }    

     return proj;


}
