#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"

/* recursive pointer function to find a vertex neighbour of a node in 
   the hierarchical data structure.                                    */

int sontype(struct image * );
int adj(int , int );
int reflect(int , int );
int c_e(int , int );
struct image *edge_neighbour(struct image * , int );

struct image *vertex_neighbour(struct image *ip, int dir)

{

      struct image *iq;

       if(ip->mother == NULL) iq=NULL; 

       else if(adj(dir, sontype(ip)) == 1) 
            iq=vertex_neighbour(ip->mother, dir);

       else if(c_e(dir, sontype(ip)) != 'O')
            iq=edge_neighbour(ip->mother, c_e(dir, sontype(ip)));

       else iq=ip->mother;

       if(iq != NULL && iq->pval == -1) {

          switch(reflect(dir, sontype(ip))) {
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
