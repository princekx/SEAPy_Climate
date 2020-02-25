#include <Stdio.h>
#include "st_obj.h"
#include "st_fo.h"
#include "grid.h"

#define  INV   -1.0

/* function to return object field values to correct values if field invert
   offset and partial invert have been used.                                */

extern float offs, ilat;
extern int y1u;

extern GRID *gr;

void object_realf(struct frame_objs *fo, float sign, int of, int hemi, int tf, int *trans)

{

    int i, j;

    float yy;

    struct object *ob=NULL;
    struct point *ptt=NULL;
    struct feature_pts *fp=NULL;

    for(i=0; i < fo->obj_num; i++){

        ob = (fo->objs) + i;

        for(j=0; j < ob->point_num; j++){

           ptt = (ob->pt) + j;

           yy = *(gr->ygrid + ptt->y + y1u - 2);

           if(trans[2]){

              if(hemi == 'n' && yy > ilat) ptt->val *= INV;
              else if(hemi == 's' && yy < ilat) ptt->val *= INV;

           }

           if(trans[0]) ptt->val *= sign;

           if(trans[1] && of == 'y') ptt->val += offs;

        }


        for(j=0; j < ob->fet->feature_num; j++){

            fp = (ob->fet->fpt)+j;

            if(fp->str < CHECK_PT) continue;

            if(tf == 3) yy = *(gr->ygrid + ((fp->y).ixy) + y1u - 2);

            else yy = (fp->y).xy;

            if(trans[2]){

               if(hemi == 'n' && yy > ilat) fp->str *= INV;
               else if(hemi == 's' && yy < ilat) fp->str *= INV;

            }

            if(trans[0]) fp->str *= sign;

            if(trans[1] && of == 'y') fp->str += offs;

        }

    }



    return;

} 

