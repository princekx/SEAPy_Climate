#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "st_obj.h"
#include "st_fo.h"
#include "mem_er.h"

void mfilt_obj(struct frame_objs *fo)
{

   int i, j;

   struct features *fet=NULL;
   struct object *ob=NULL;
   struct feature_pts *fpt=NULL;
   struct feature_pts ftemp;

   fo->tot_f_f_num = 0;

   for(i=0; i < fo->obj_num; i++){

       memset(&ftemp, 0, sizeof(struct feature_pts));

       ob = fo->objs + i;
       fet = ob->fet;

       fpt = fet->fpt;


       if(fet->feature_num > 0){

          for(j=0; j < fet->feature_num; j++){

              if((fpt+j)->str > ftemp.str) memcpy(&ftemp, fpt+j, sizeof(struct feature_pts));

          }

          fet->fpt = (struct feature_pts *)realloc_n(fet->fpt, sizeof(struct feature_pts));
          mem_er((fet->fpt == NULL) ? 0 : 1, sizeof(struct feature_pts));

          memcpy(fet->fpt, &ftemp, sizeof(struct feature_pts));

          fet->feature_num = 1;
          ++(fo->tot_f_f_num);


       }

   }

   return;

}
