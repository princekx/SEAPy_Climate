#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"

/* a recursive pointer function to find an edge neigbour of a node in
   the hierarchical data tree                                         */

int sontype(struct image * );
int reflect(int , int );
int adj(int , int );

struct image *edge_neighbour(struct image *ih, int dir)

{

    struct image *iq;

    if(ih->mother != NULL && adj(dir, sontype(ih)) == 1)

       iq = edge_neighbour(ih->mother, dir);

    else

       iq = ih->mother;

    if(iq != NULL && iq->pval == -1) {

       switch(reflect(dir, sontype(ih))) {
       case 1:
         iq=iq->nw;
         break;
       case 2:
         iq=iq->ne;
         break;
       case 3:
         iq=iq->sw;
         break;
       case 4:
         iq=iq->se;
         break;
       }

       return iq;

    }

    else return iq;
}
