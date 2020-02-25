#include <Stdio.h>
#include <stdlib.h>
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"

extern int track_num;

/* function to move a portion of track to another free section of track
   when the feature point seperation criterion is violated.              */

void track_fail(struct frame_objs *fo, struct track_ind *tind, int t1, int k, int fr, int dir)

{

   int i, j;
   int fst=0, lst=0, ptl=0;
   int fs, ls;
   int exfl=0;

   struct track_ind *tr1, *tr2;
   struct track_points *ts, *tf;

   tr1 = tind + t1;

   if(dir == 'f'){

     fst = k + 2;

     if(fst == fr-1) lst = fr;

     else{

        lst = fr;

        for(i=fst+1; i < fr; i++){

           if(((tr1->tp)+i)->feature_id == -1) {lst = i; break;}


        }


     }

     ptl = lst - fst;

   }

   else if(dir == 'b'){

      lst = k - 1;

      if(lst-1 == 0) fst = 0;

      else{

         fst = 0;

         for(i=lst-1; i >= 0; i--){

            if(((tr1->tp)+i)->feature_id == -1) {fst = i+1; break;}

         }

      }

      ptl = lst - fst;

    }


    fs = (fst-1 >= 0) ? fst-1 : fst;
    ls = (lst == fr) ? lst-1 : lst;

    for(i=0; i < track_num; i++){

       tr2 = tind + i;
       exfl = 1;

       if(((tr2->tp)+fs)->feature_id == -1){

            exfl = 0;
            j = fs+1;

            while(j <= ls){

                 if(((tr2->tp)+j)->feature_id != -1) {exfl = 1; break;}

                 j++;

            }


        }

        if(exfl == 0){
                    
            tr1->num_point_track -= ptl;
            tr2->num_point_track += ptl;

            for(j=fst; j < lst; j++){

               tf = (tr1->tp)+j;
               ts = (tr2->tp)+j;
               ts->object_id = tf->object_id;
               ts->feature_id = tf->feature_id;
               tf->object_id = 0;
               tf->feature_id = -1;

            }

            break;

        }

     }

     if(exfl == 1){

        for(i=0; i < track_num; i++){

           tr2 = tind + i;

           if(tr2->num_point_track == 0){

              tr1->num_point_track -= ptl;
              tr2->num_point_track += ptl;

              for(j=fst; j < lst; j++){

                  tf = (tr1->tp)+j;
                  ts = (tr2->tp)+j; 
      
                  ts->object_id = tf->object_id;
                  ts->feature_id = tf->feature_id;
                  tf->object_id = 0;
                  tf->feature_id = -1;

              }

              break;

           }

        }

    }

    return;

}
