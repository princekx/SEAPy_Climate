#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"

/* function to determine the son-type of a node in the hierachical data tree 
   for north west quadrant sontype = 1
    "  north east     "        "   = 2
    "  south west     "        "   = 3
    "  south east     "        "   = 4                                        */

int sontype(struct image *ih)

{

    struct image *ip;
    int type=0;

    ip=ih->mother;
    if(ip->nw == ih) type=1;
    else if(ip->ne == ih) type=2;
    else if(ip->sw == ih) type=3;
    else if(ip->se == ih) type=4;
    else {
       printf("***ERROR***: file \"%s\", line %d\n", __FILE__, __LINE__);
       exit(1);
    }

    return type;

}
