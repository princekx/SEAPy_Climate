#include <Stdio.h>
#include <stdlib.h>
#include "st_obj.h"
#include "st_fo.h"
#include "mem_er.h"
#include "statistic.h"
#include "grid.h"
#include "boundary.h"
#include "proj.h"


/* function to check feature points with objects for re-read object data */

int boundary_find(struct object * , struct boundary_cntl * , int * , int * , int );
int inside(int , struct boundary_pt * , float , float );

extern int tf;
extern int x1u, y1u;

extern GRID *gr;

extern float period;

extern PROJ *pp;

void feature_to_object(struct frame_objs *ff, struct frame_objs *fr, int frd)

{

    int i, j, k;
    int in_obj=0;
    int ic=0;
    int ift=0;

    struct object *ob1=NULL, *ob2=NULL;
    struct feature_pts *fpts1=NULL, *fpts2=NULL;
    struct boundary_pt *btt=NULL;

    float xx, yy;


    for(k=0; k < ff->obj_num; k++){

        ob2 = ff->objs + k;

        boundary_find(ob2, NULL, NULL, NULL, 1);
/* convert boundary to float format */


        for(j=0; j < ob2->bound_num; j++){

           btt = ob2->bound + j;

           (btt->x).xy = obj_xreal((btt->x).ixy + x1u - 2);

           (btt->y).xy = *(gr->ygrid + (btt->y).ixy + y1u - 2);                     

        }

        for(i=0; i < fr->obj_num; i++){

            ob1 = fr->objs + i;

            if(ob1->fet->feature_num){

               for(j=0; j < ob1->fet->feature_num; j++){

                   fpts1 = (ob1->fet->fpt) + j;

                   in_obj = 0;

/*printf("%d %d %d\n", i, j, in_obj); */

                   if(tf == 3){

                      xx = obj_xreal((fpts1->x).ixy + x1u - 2);
                      yy = *(gr->ygrid + (fpts1->y).ixy + y1u - 2);

                   }

                   else {

                      xx = (fpts1->x).xy;
                      yy = (fpts1->y).xy;

                   } 

                   if(frd == 2 && gr->prty){

                       switch(gr->prgr){
                            case 0:
                              (pp->prj2)(&xx, 'x', 0);
                              (pp->prj2)(&yy, 'y', 0);
                              break;
                            case 1:
                               ic=azimuthal(&yy, &xx, gr->alat, gr->alng, 0, gr->prty);
                               break;
			    case 2:
			       ic=conic(&yy, &xx, gr->alat, gr->alng, gr->sp1, gr->sp2, 0, gr->prty, &ift);
			       break;
			    case 3:
			       ic=rotated_cyln(&yy, &xx, gr->alat, gr->alng, 0, &ift);
			       break;
                        }

                   }  

                   if(inside(ob2->bound_num, ob2->bound, xx, yy)) in_obj = 1;
                   else if(inside(ob2->bound_num, ob2->bound, xx + period, yy)){in_obj = 1; xx += period;}
                   else if(inside(ob2->bound_num, ob2->bound, xx - period, yy)){in_obj = 1; xx -= period;}  

/* printf("%d %d %d\n", i, j, in_obj); */            

                   if(in_obj){

                      if(!ob2->fet){
                         ob2->fet = (struct features *)malloc_initl(sizeof(struct features));
                         mem_er((ob2->fet == NULL) ? 0 : 1, sizeof(struct features));
                         ob2->fet->fpt = (struct feature_pts *)malloc_initl(sizeof(struct feature_pts));
                         mem_er((ob2->fet->fpt == NULL) ? 0 : 1, sizeof(struct feature_pts));

                         ob2->fet->feature_num = 1;
                         ff->tot_f_f_num = 1;
                      }

                      else {
                         ++(ob2->fet->feature_num);
                         ++(ff->tot_f_f_num);
                         ob2->fet->fpt = (struct feature_pts *)realloc_n(ob2->fet->fpt, (ob2->fet->feature_num)*sizeof(struct feature_pts));
                         mem_er((ob2->fet->fpt == NULL) ? 0 : 1, (ob2->fet->feature_num)*sizeof(struct feature_pts));
                      }

                      fpts2 = (ob2->fet->fpt) + ob2->fet->feature_num - 1;
                      (fpts2->x).xy = xx;
                      (fpts2->y).xy = yy;
                      fpts2->str = fpts1->str;

                      fpts2->iob = i;
                      fpts2->ifet = j;

/* printf("%d %d\n", fpts2->iob, fpts2->ifet); */

                   }

/* printf("-----%d %d %d\n", i, j, in_obj); */


               }         

           }

        }

        free(ob2->bound);
        ob2->bound_num = 0;

        if(!(ob2->fet)){
           free(ob2->pt);
           ob2->point_num = 0;
        }

    }

    return;

}


void putback_area(struct frame_objs *ff, struct frame_objs *fr)
{

    int i,k;

    struct feature_pts *fpts1=NULL, *fpts2=NULL;
    struct object *ob1=NULL, *ob2=NULL;

    for(k=0; k < ff->obj_num; k++){

        ob2 = ff->objs + k;

        for(i=0; i<ob2->fet->feature_num; i++){

            fpts2 = (ob2->fet->fpt) + i;

            ob1 = fr->objs + fpts2->iob;
            fpts1 = ob1->fet->fpt + fpts2->ifet;

            if(tf == 3){
               (fpts1->ox).ixy = (fpts2->ox).ixy;
               (fpts1->oy).ixy = (fpts2->oy).ixy;
            }
            else {
               (fpts1->ox).xy = (fpts2->ox).xy;
               (fpts1->oy).xy = (fpts2->oy).xy;
               fpts1->ostr = fpts2->ostr;
            }

            fpts1->area = fpts2->area;
            fpts1->r_sh = fpts2->r_sh;
            fpts1->ornt[0] = fpts2->ornt[0];
            fpts1->ornt[1] = fpts2->ornt[1];

        }

    }


    return;

}
