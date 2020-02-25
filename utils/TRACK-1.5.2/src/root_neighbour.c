#include <Stdio.h>
#include "sqt.h"

SQT *root_neighbour(SQT *sq, int dir, int iplat, int iroot)
{

    SQT *sqq=NULL;


    if(!iplat){

/* octahedron */

       if(iroot >= 0 && iroot < 4){

          if(dir == 'L'){
             sqq = (iroot - 1 >= 0) ? sq - 1: sq + 3;
          }
          else if(dir == 'R'){
             sqq = (iroot + 1 <= 3) ? sq + 1: sq - 3;  
          }
          else if(dir == 'V'){
             sqq = sq + 4;
          }

       }
       else {

          if(dir == 'L'){
             sqq = (iroot - 1 >= 4) ? sq - 1: sq + 3;
          }
          else if(dir == 'R'){
             sqq = (iroot + 1 <= 7) ? sq + 1: sq - 3;  
          }
          else if(dir == 'V'){
             sqq = sq - 4;
          }

       }


    }

    else {

/* icosahedron */

       if(iroot >= 0 && iroot < 5){

          if(dir == 'L'){
             sqq = (iroot - 1 >= 0) ? sq - 1: sq + 4;
          }
          else if(dir == 'R'){
             sqq = (iroot + 1 <= 4) ? sq + 1: sq - 4;  
          }
          else if(dir == 'V'){
             sqq = sq + 5;
          }

       }
       else if(iroot >= 5 && iroot < 15){

          if(dir == 'L'){
             sqq = (iroot - 1 >= 5) ? ((sq->ud == 'd') ? sq + 4 : sq - 5) : sq + 9;
          }
          else if(dir == 'R'){
             sqq = (iroot + 1 <= 14) ? ((sq->ud == 'd') ? sq + 5 : sq - 4) : sq - 9;  
          }
          else if(dir == 'V'){
             if(sq->ud == 'd') sqq = sq - 5;
             else sqq = sq + 5 ;
          }


       }

       else {

          if(dir == 'L'){
             sqq = (iroot - 1 >= 15) ? sq - 1: sq + 4;
          }
          else if(dir == 'R'){
             sqq = (iroot + 1 <= 19) ? sq + 1: sq - 4;  
          }
          else if(dir == 'V'){
             sqq = sq - 5;
          }

       }


    }

    return sqq;

}
