#include <Stdio.h>
#include <stdlib.h>

int c_e(int dir, int ntype)

{

    int c='O';

    if(dir == 'V'){

       switch(ntype){
       case 1: case 4:
          c= 'O';
          break;
       case 2:
          c= 'N';
          break;
       case 3:
          c= 'W';
          break;
       }
     }

     else if(dir == 'U'){
       switch(ntype){
       case 1:
          c= 'N';
          break;
       case 2: case 3:
          c= 'O';
          break;
       case 4:
          c= 'E';
          break;
       }
    }

    else {
       printf("****ERROR****, tag %c not known for neigbour finding.\n\n", dir);
       exit(1);
    }

    return c;

}
