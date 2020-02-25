#include <Stdio.h>
#define  SUM4(a, b, c, d)  (int)(a+b+c+d)


/* function to assign node colour 1 for black 0 for white and -1 for grey */

int assign_node(int sw, int se, int nw, int ne)
{

     int nnode, sum;

     if(sw >= 0 && se >= 0 && nw >= 0 && ne >= 0) {

           sum = SUM4(sw, se, nw, ne);

           if(sum == 4 )

              nnode = 1;

           else if(sum == 0)

              nnode = 0;

           else

             nnode = -1;

      }

      else

         nnode = -1;


      return nnode;

}
