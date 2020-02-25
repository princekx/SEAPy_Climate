#include <Stdio.h>
#include <stdlib.h>
#include "st_obj.h"
#include "st_fo.h"


void fpt_intensity(struct frame_objs *ff, float upint1, float upint2)
{
   int i, j, k;

   struct object *ob=NULL;
   struct feature_pts *fpts=NULL, *fptsn=NULL;

   for(i=0; i < ff->obj_num; i++){

      ob = (ff->objs) + i;

      fpts = ob->fet->fpt;

      for(j=0; j < ob->fet->feature_num; j++){ 
 
          fptsn = fpts + j;
          if((fptsn->str - upint1) * (upint2 - fptsn->str) < 0){
             for(k=j+1; k<ob->fet->feature_num; k++){

                 *(fpts + k - 1) = *(fpts + k);

             }

             --(ob->fet->feature_num);
             --(ff->tot_f_f_num);
             --j;

          }

      }

   }

   return;

}
