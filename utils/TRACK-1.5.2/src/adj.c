#include <Stdio.h>
#include <stdlib.h>

/* logical function to determine the adjacency of nodes in the
   hierarchical data tree                                          */

int adj(int dir, int ntype)

{

      int logical = 0;

      switch(dir) {

         case 'N':

           logical = (ntype == 1 || ntype == 2) ? 1 : 0;
           break;
      
         case 'E':

           logical = (ntype == 2 || ntype == 4) ? 1 : 0;
           break;

         case 'S':

           logical = (ntype == 3 || ntype == 4) ? 1 : 0;
           break;

         case 'W':

           logical = (ntype == 1 || ntype == 3) ? 1 : 0;
           break;

         case 'V':

           logical = (ntype == 1) ? 1 : 0;
           break;

         case'U':

           logical = (ntype == 2) ? 1 : 0;
           break;

         default:
           printf("****ERROR****, tag %c not known for neigbour finding in quad tree.\n\n", dir);
           exit(1);

      }

      return logical;

}
