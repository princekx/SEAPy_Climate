#include <Stdio.h>
#include <stdlib.h>
#include "sqt.h"

/* logical function to determine adjacency of nodes in a spherical quad tree */


int sqt_adj(int dir, int type)
{

    int lgcal = 0;

    switch(dir) {

        case 'L':
           lgcal = (type == 1 || type == 2) ? 1 : 0;
           break;
        case 'R':
           lgcal = (type == 1 || type == 4) ? 1 : 0;
           break;
        case 'V':
           lgcal = (type == 2 || type == 4) ? 1 : 0;
           break;
        default:
           printf("****ERROR****, tag %c not known for neigbour finding.\n\n", dir);
           exit(1);
    }

    return lgcal;

}
