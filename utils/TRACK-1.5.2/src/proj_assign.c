#include <Stdio.h>
#include <proj.h>

/* function to assign projection */

PRFP proj_assign(int prj)

{

   PRFP proj=NULL;

   switch(prj){
         case 0:
            break;
         case 1:
            proj = eq_rect;
            break;
         case 2:
            proj = eq_area;
            break;
         case 3:
            proj = urmaev;
            break;
         case 4:
            proj = mercator;
            break;
         case 5:
            proj = stereog;
            break;
         case 6:
            proj = miller;
            break;
         default:
            printf("no projection available for chosen flag\n");
            break;
   }

   return proj;

}
