#include <Stdio.h>
#include <stdlib.h>
#include "sqt.h"

int sqt_reflect(int dir, int type)
{

    int typ=0;

    if(dir == 'L'){
       if(type == 1) typ = 4;
       else if(type == 2) typ = 1;
       else if(type == 3) typ = 2;
       else if(type == 4) typ = 3;
    }
    else if(dir == 'R'){
       if(type == 1) typ = 2;
       else if(type == 2) typ = 3;
       else if(type == 3) typ = 4;
       else if(type == 4) typ = 1;
    }
    else if(dir == 'V'){
       if(type == 1) typ = 3;
       else if(type == 2) typ = 2;
       else if(type == 3) typ = 1;
       else if(type == 4) typ = 4;
    }
    else {
       printf("****ERROR****, tag %c not known for neigbour finding.\n\n", dir);
       exit(1);
    }

    return typ;

}

int reflect_root(int dir, int type)
{
   int typ=0;

    if(dir == 'L'){
       if(type == 1) typ = 1;
       else if(type == 2) typ = 4;
    }
    else if(dir == 'R'){
       if(type == 1) typ = 1;
       else if(type == 4) typ = 2;
    }
    else if(dir == 'V'){
       if(type == 2) typ = 2;
       else if(type == 4) typ = 4;
    }
    else {
       printf("****ERROR****, tag %c not known for neigbour finding.\n\n", dir);
       exit(1);
    }   



   return typ;

}
