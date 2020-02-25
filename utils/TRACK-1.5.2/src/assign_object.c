#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_im.h"
#include "st_fo.h"
#include "boundary.h"
#include "st_obj.h"

void assign_points_to_objects(struct object * , int, int , int , int );
int powi(int , int );


/* function to assign storage for objects and the consituitive points */


void assign_object(struct image *iq, struct frame_objs *ff, int ll, int xp, int yp)


{

    int num=powi(powi(2, ll), 2);
    int i, plus, cpt, id=0;

    struct object *ob=NULL;


    if(ff->objs != NULL){

       for(i=0; i < ff->obj_num; i++){

           ob = (ff->objs)+i;

           if(ob->lab == iq->label->name){

             id=1;
             break;
           }

       }
    }


     if(id == 0){

        ++(ff->obj_num);
        plus=ff->obj_num;

        if(plus == 1) {

           ff->objs=(struct object *)malloc_initl(sizeof(struct object));
           mem_er((ff->objs == NULL) ? 0 : 1, sizeof(struct object));

        }

        else {

             ff->objs=(struct object *)realloc_n(ff->objs, plus*sizeof(struct object));
             mem_er((ff->objs == NULL) ? 0 : 1, plus*sizeof(struct object));

        }

        ob=(ff->objs)+plus-1;
        ob->ext = NULL;
        ob->fet = NULL;
        ob->bound = NULL;

        ob->mem_pt_size = 0;
        ob->point_num = 0;
        ob->bound_num = 0;
        ob->b_or_i = '\0';

        ob->lab = iq->label->name;
        ob->pt = (struct point *)calloc(num, sizeof(struct point));
        mem_er((ob->pt == NULL) ? 0 : 1, num * sizeof(struct point));

        cpt = ob->point_num;
        ob->point_num = num;


        assign_points_to_objects(ob, cpt, ll, xp, yp);

     }

     else if(id == 1){

         cpt = ob->point_num;
         (ob->point_num)+=num;

         ob->pt = (struct point *)realloc_n(ob->pt, (ob->point_num)*sizeof(struct point));
         mem_er((ob->pt == NULL) ? 0 : 1, (ob->point_num)*sizeof(struct point));

        assign_points_to_objects(ob, cpt, ll, xp, yp);

     }

     return;

}
  
