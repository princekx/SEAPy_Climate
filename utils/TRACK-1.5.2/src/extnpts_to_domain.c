#include <Stdio.h>
#include "st_obj.h"
#include "st_fo.h"
#include "boundary.h"
#include "grid.h"

extern int x1u, y1u;
extern int tf;
extern int pb;

extern float period;

extern GRID *gr;

/* function to push external object and feature points back into the
   domain following boundary object external point calculations for
   periodic boundary conditions in X-direction.                       */

void extnpts_to_domain(struct frame_objs *fo, int *llb, int *rrb)

{

  int i, j;
  int ix1=gr->ix - 1;
  int ixx, iyy;
  int ipt;

  float xx1=*(gr->xgrid), xx2=*(gr->xgrid+ix1);

  struct object *ob;
  struct point *ptt, *pt1;
  struct feature_pts *fpt;
  struct boundary_pt *btt=NULL;

  if(pb != 'y') {

    printf("***WARNING***, periodic boundary condition not set for file %s\n", __FILE__);
    return;

  }


  for(i=0; i<fo->obj_num; i++){

      ob = (fo->objs) + i;
 
      if(ob->b_or_i == 'l' || ob->b_or_i == 'r' || ob->b_or_i =='b' ){

         ipt = ob->point_num;

         for(j=0; j<ipt; j++){
 
             ptt = (ob->pt) + j;

             ixx = ptt->x + x1u - 1;
             iyy = ptt->y + y1u - 2;

             if(ixx <= 1) {

                if(ixx < 1) ptt->x += ix1;
                else if(ob->point_num < ob->mem_pt_size &&
                        *(llb+iyy) != *(rrb+iyy)           ){ 
                   pt1 = (ob->pt)+ob->point_num;
                   *pt1 = *ptt;
                   pt1->x += ix1;
                   ++(ob->point_num);

                }

             }

             else if(ixx >= gr->ix){

                if(ixx > gr->ix) ptt->x -= ix1;
                else if(ob->point_num < ob->mem_pt_size && 
                        *(llb+iyy) != *(rrb+iyy)           ){
                   pt1 = (ob->pt)+ob->point_num;
                   *pt1 = *ptt;
                   pt1->x -= ix1;
                   ++(ob->point_num);

                }

             }


         }

         if(ob->bound){


            if(!(fo->b_state)){
               for(j=0; j<ob->bound_num; j++){

                   btt = ob->bound + j;
                   ixx = (btt->x).ixy + x1u - 1;

                   if(ixx < 1) (btt->x).ixy += ix1;
                   else if(ixx > gr->ix) (btt->x).ixy -= ix1;

               }

            }

            else {

               for(j=0; j<ob->bound_num; j++){
                   btt = ob->bound + j;
                   if((btt->x).xy < xx1) (btt->x).xy += period;
                   else if((btt->x).xy > xx2) (btt->x).xy -= period;

               }

            }

         }

         if(tf){

            if(ob->fet){
 
               for(j=0; j < ob->fet->feature_num; j++){

                   fpt = (ob->fet->fpt) + j;

                   if(tf == 3){

                     if((fpt->x).ixy < 1) (fpt->x).ixy += (ix1 - x1u + 1);
                     else if((fpt->x).ixy > gr->ix) (fpt->x).ixy  -= (ix1 + x1u -1);

                   }
  
                   else {

                     if((fpt->x).xy < xx1) (fpt->x).xy += period;
                     else if((fpt->x).xy > xx2) (fpt->x).xy -= period;

                   }

               }

            }

         }


      }

   }

   for(i=0; i < gr->iy; i++) *(llb + i) = *(rrb + i) = 0;

   return;

}
