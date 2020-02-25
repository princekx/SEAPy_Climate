#include <Stdio.h>
#include <stdlib.h>

/* function to determine the sontype that shares an edge with the 
   specified sontype for the hierarchical data tree                */

int reflect(int dir, int ntype)

{

    int type=0;

    if(dir == 'N' || dir == 'S'){

          switch(ntype){
          case 1: case 2:
              type=ntype+2;
              break;
          case 3: case 4:
              type=ntype-2;
              break;
          }
     }

     else if(dir == 'E' || dir == 'W') {

          switch(ntype){
          case 1: case 3:
              type=ntype+1;
              break;
          case 2: case 4:
              type=ntype-1;
              break;
          }
     }

     else if(dir == 'V' || dir == 'U') {

          switch(ntype){
          case 1: 
              type=4;
              break;
          case 2:
              type=3;
              break;
          case 3:
              type=2;
              break;
          case 4:
              type=1;
              break;
          }
     }

     else {
       printf("****ERROR****, tag %c not known for neigbour finding.\n\n", dir);
       exit(1);
     }

     return type;

}
