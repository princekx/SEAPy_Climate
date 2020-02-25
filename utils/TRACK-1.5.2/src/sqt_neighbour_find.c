#include <Stdio.h>
#include "sqt.h"

SQT *sqt_neighbour(SQT * , int , int , int , int , int , int * );

void sqt_neighbour_find(SQT **sqt, int nlev, int nlf, int iplat, int nfc)
{

   int i,j;
   int iroot=0;
   int irf = 0;

   LEAF *lf=NULL, *ll=NULL, *lnew=NULL;
   SQT *sq=NULL;

   lf = (LEAF *)sqt[nlev-1];

   for(i=0; i < nlf; i++){

       ll = lf + i;

/* find root trinagle ID */

       sq = ll->parent;
       iroot = 0;
       while(sq->parent) sq = sq->parent;

       for(j=0; j < nfc; j++) {
           if(sq == *sqt + j){iroot = j; break;}
       }       

       if(!(ll->lft)) {
          irf = 0;
          lnew = (LEAF *)sqt_neighbour((SQT *)ll, 'L', nlev, nlev, iplat, iroot, &irf);
          ll->lft = lnew;
          lnew->rgt = ll;   
 

/*printf("%p %p\n",  (void *)(ll->lft), (void *)(lnew->rgt));
printf("L %d %d %d -- %d %d %d\n", lnew->ivec[0],  lnew->ivec[1], lnew->ivec[2], ll->ivec[0], ll->ivec[1], ll->ivec[2]); */


       }
       if(!(ll->rgt)) {
          irf = 0;
          lnew = (LEAF *)sqt_neighbour((SQT *)ll, 'R', nlev, nlev, iplat, iroot, &irf);
          ll->rgt = lnew;
          lnew->lft = ll;
 

/*printf("%p %p\n",  (void *)(ll->rgt), (void *)(lnew->lft));
printf("R %d %d %d -- %d %d %d\n", lnew->ivec[0],  lnew->ivec[1], lnew->ivec[2], ll->ivec[0], ll->ivec[1], ll->ivec[2]); */


       }
       if(!(ll->udw)) {
          irf = 0;
          lnew = (LEAF *)sqt_neighbour((SQT *)ll, 'V', nlev, nlev, iplat, iroot, &irf);
          ll->udw = lnew;
          lnew->udw = ll;
 

/*printf("%p %p\n",  (void *)(ll->udw), (void *)(lnew->udw));
printf("V %d %d %d -- %d %d %d\n", lnew->ivec[0],  lnew->ivec[1], lnew->ivec[2], ll->ivec[0], ll->ivec[1], ll->ivec[2]); */


       }

   }

   return;

}
