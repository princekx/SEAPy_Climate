#include <Stdio.h>


void proj_report(int pty, int pgr)

{


    switch(pgr){
         case 0:
           printf("\n***Cylindrical Projection***\n");
           switch(pty){
              case 0:
                 printf("Current world projection is Plate Caree\n\n");
                 break;
              case 1:
                 printf("Current world projection is Equirectangular\n\n");
                 break;
              case 2:
                 printf("Current world projection is Equal Area\n\n");
                 break;
              case 3:
                 printf("Current world projection is Urmaev\n\n");
                 break;
              case 4:
                 printf("Current world projection is Mercator\n\n");
                 break;
              case 5:
                 printf("Current world projection is Stereographic\n\n");
                 break;
              case 6:
                 printf("Current world projection is Miller\n\n");
                 break;
              default:
                 printf("***error***, no world projection defined for this identifier in this group\n\n");
                return;

            }
            break;

         case 1:
           printf("\n***Azimuthal Projection***\n");
           switch(pty){
              case 1:
                 printf("Current world projection is Orthographic\n\n");
                 break;
              case 2:
                 printf("Current world projection is Stereographic\n\n");
                 break;
              case 3:
                 printf("Current world projection is Gnomnic\n\n");
                 break;
              case 4:
                 printf("Current world projection is Lambert Equal Area\n\n");
                 break;
              case 5:
                 printf("Current world projection is Azimuthal Equidistant\n\n");
                 break;
              default:
                 printf("***error***, no world projection defined for this identifier in this group\n\n");
                return;

            }
            break;
         case 2:
           printf("\n***Conic Projection***\n");
           switch(pty){
              case 1:
                 printf("Current world projection is Albers Equal Area\n\n");
                 break;
              case 2:
                 printf("Current world projection is Lambert Conformal\n\n");
                 break;
              case 3:
                 printf("Current world projection is Equidistant Conic\n\n");
                 break;	   	 
               default:
                 printf("***error***, no world projection defined for this identifier in this group\n\n");
                return; 
            }
            break; 
         case 3:
           printf("\n***Rotated Cylindrical Projection***\n");
           switch(pty){
              case 0:
                 printf("Current world projection is Plate Caree\n\n");
                 break;
              case 1:
                 printf("Current world projection is Equirectangular\n\n");
                 break;
              case 2:
                 printf("Current world projection is Equal Area\n\n");
                 break;
              case 3:
                 printf("Current world projection is Urmaev\n\n");
                 break;
              case 4:
                 printf("Current world projection is Mercator\n\n");
                 break;
              case 5:
                 printf("Current world projection is Stereographic\n\n");
                 break;
              case 6:
                 printf("Current world projection is Miller\n\n");
                 break;
              default:
                 printf("***error***, no world projection defined for this identifier in this group\n\n");
                return;

            }
            break;	          
     }

     return;

}
