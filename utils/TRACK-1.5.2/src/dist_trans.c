#include "mem_er.h"

/* function to dilate an object, uses same connectivity as object finding.
   Uses distance transform on the inverse binary image of objects.        */


#define  MIN(A, B) (int)((A < B) ? A : B)
#define  MIN4(A, B, C) (int)(MIN(A, MIN(B, C)))
#define  MIN8(A, B, C, D, E) (int)(MIN(A, MIN(B, MIN(C, MIN(D, E)))))
#define  MINI4(A, B, C) (int)((A > 0) ? MIN(B, C) : 0)
#define  MINI8(A, B, C, D, E) (int)((A > 0) ? MIN(B, MIN(C, MIN(D, E))) : 0)


extern int cc;


void dist_trans(int *oba, int xd, int yd, int mth)

{

   int i,j;
   int *aa, *ab, *ac;

   if(cc == 'e'){

      for(i=yd-2; i>0; i--){

          ac = oba + i * xd;
          aa = ac + 1;
          ab = oba + (i + 1) * xd + 1;


          for(j=1; j<xd-1; j++){

              *aa = MINI4(*aa, *ac + 1, *ab + 1);

              ++aa; ++ac; ++ab;


          }

     }

     for(i=2; i < yd; i++){

         ac = oba + i * xd - 1;
         aa = ac - 1;
         ab = oba + (i-1) * xd - 2;

         for(j=xd-2; j >0; j--){

             *aa = MIN4(*aa, *ac + 1, *ab + 1);
             --aa; --ac; --ab;



         }


     }


   }

   else if(cc == 'v'){

      mth *= 2;

      for(i=yd-2; i>0; i--){

          ac = oba + i * xd;
          aa = ac + 1;
          ab = oba + (i + 1) * xd + 1;


          for(j=1; j<xd-1; j++){

              *aa = MINI8(*aa, *ac + 2, *ab + 2, *(ab-1)+3, *(ab+1)+3);

              ++aa; ++ac; ++ab;


          }

     }

     for(i=2; i < yd; i++){

         ac = oba + i * xd - 1;
         aa = ac - 1;
         ab = oba + (i-1) * xd - 2;

         for(j=xd-2; j >0; j--){

             *aa = MIN8(*aa, *ac + 2, *ab + 2, *(ab+1)+3, *(ab-1)+3);
             --aa; --ac; --ab;



         }


     }



   }


   if(mth){

      for(i=0; i<xd*yd; i++) {
         aa = oba + i; 
         if((*aa) <= mth) *aa = 0;
         else *aa -= mth;
 

      }

   }

   return;

}
