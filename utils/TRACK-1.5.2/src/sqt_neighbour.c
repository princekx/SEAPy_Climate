#include <Stdio.h>
#include "sqt.h"

/* function to find leaf neighbours in an SQT */

int sqt_sontype(SQT * , int , int );
int sqt_adj(int , int );
int sqt_reflect(int , int );
int reflect_root(int , int );
SQT *root_neighbour(SQT * , int , int , int );

SQT *sqt_neighbour(SQT *sq, int dir, int ilv, int mxl, int iplat, int iroot, int *irf)
{

   SQT *sqt=NULL;
   LEAF *ll=NULL;

   if(ilv == mxl){
      ll = (LEAF *)sq;

      if(sqt_adj(dir, sqt_sontype(sq, ilv, mxl)))

         sqt = sqt_neighbour(ll->parent, dir, ilv - 1, mxl, iplat, iroot, irf);

      else 

         sqt = ll->parent;

   }

   else {

      if(sq->parent){

         if(sqt_adj(dir, sqt_sontype(sq, ilv, mxl)))

            sqt = sqt_neighbour(sq->parent, dir, ilv - 1, mxl, iplat, iroot, irf);

         else sqt = sq->parent;

      } 

      else {*irf = 1; return root_neighbour(sq, dir, iplat, iroot);}

   }

   if((iplat && *irf) && (dir == 'L' || dir == 'R') && (iroot >= 5 && iroot <= 14)) *irf = 0;

   if(ilv == mxl - 1 && sqt){

      switch((*irf) ? reflect_root(dir, sqt_sontype(sq, ilv, mxl)) : 
                      sqt_reflect(dir, sqt_sontype(sq, ilv, mxl))) {
         case 1:
            sqt = (SQT *)((LEAF *)(sqt->tt[0]));
            break;
         case 2:
            sqt = (SQT *)((LEAF *)(sqt->tt[1]));
            break;
         case 3:
            sqt = (SQT *)((LEAF *)(sqt->tt[2]));
            break;
         case 4:
            sqt = (SQT *)((LEAF *)(sqt->tt[3]));
            break;
      }

   }

   else if(sqt){

      switch((*irf) ? reflect_root(dir, sqt_sontype(sq, ilv, mxl)) :
                      sqt_reflect(dir, sqt_sontype(sq, ilv, mxl))) {
         case 1:
            sqt = sqt->tt[0];
            break;
         case 2:
            sqt = sqt->tt[1];
            break;
         case 3:
            sqt = sqt->tt[2];
            break;
         case 4:
            sqt = sqt->tt[3];
            break;
      }

   }

   return sqt;

}
