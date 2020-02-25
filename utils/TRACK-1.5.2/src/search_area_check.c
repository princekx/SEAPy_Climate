#include <Stdio.h>
#include <stdlib.h>

/* check user defined search area */

int search_area_check(int u1, int u2, int mdim)

{

     int snum, dim=u2-u1+1; 

     if((u1 >= u2) || (u1< 1) || (u2>mdim)) {

        printf("***ERROR***, chosen numbers are out of bounds\n");

        exit(1);

     }

     switch(dim){

        case 4:   case 8:   case 16:   case 32:   case 64:
        case 128: case 256: case 512:  case 1024: case 2048:

            printf("****INFORMATION****, chosen search area is consistent with data hierachy\n\n");

            snum=1;
            break;

         default:

            printf("****INFORMATION****, chosen search area is inconsistent with data \r\n "
                   "                     hierachy, but will be padded out to continue\n\n ");

            snum=0;
            break;

      }

      return snum;

}
