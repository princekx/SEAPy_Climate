#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_fo.h"
#include "boundary.h"
#include "st_obj.h"
#include "proj.h"
#include "grid.h"


/* function to put filtered objects into arrays for plotting using UNIRAS */

#ifdef  NOUNDERSCORE

void filterpl(float *, float *, float * , int *, int *, int *, int *, 
               float *, int *, int *, int *, float *, float *, int *, int *, 
               int *, int *, int *, int *, float *, float *, float *, int *,
               float *, float *, float *, int * , int * , int * , float * ,
               float * , float * , float * , float * , float * , int * , int * , int * );

#else

void filterpl_(float *, float *, float * , int *, int *, int *, int *, 
               float *, int *, int *, int *, float *, float *, int *, int *, 
               int *, int *, int *, int *, float *, float *, float *, int *,
               float *, float *, float *, int * , int * , int * , float * ,
               float * , float * , float * , float * , float * , int * , int * , int * );

#endif


extern int dfil;
extern int x1u, x2u, y1u, y2u;

extern GRID *gr, *gr1, *gr2;
extern CNTRY *cm;

extern int iian;
extern int *xan, *yan;
extern int an_plt;
extern PROJ *pp;

void filtd(struct frame_objs *ff, int tf, int fruser, int frcount, float thresh)

{

     static int ibi=0;
     static int pbb=0;

     int i, j, k, istart=0, ic=0;
     int totb=0;
     int ibc=0;
     int ift=0;

     int gd=gr->ix*gr->iy;

     int vp=1;

     int *xpos=NULL, *ypos=NULL;
     int *ibb=NULL;

     float  *valf=NULL;                   /* filtered object data */

     float  *z=NULL;
     float  *x=NULL, *y=NULL, *zz=NULL;

     float *xfeat=NULL, *yfeat=NULL;

     float *xor=NULL, *yor=NULL, *zor=NULL;
     float xx, yy;

     float *xbb=NULL, *ybb=NULL, *vbb=NULL;

     struct object *ob=NULL;
     struct point *ptt=NULL;
     struct features *fl=NULL;
     struct boundary_pt *btt=NULL;

     PROJ proj={NULL, NULL};
     
     GRDP g1={0, 0, 0.0, 0.0, 0.0, 0.0};
     GRDP g2={0, 0, 0.0, 0.0, 0.0, 0.0};
     GRDP g22={0, 0, 0.0, 0.0, 0.0, 0.0};

     if(!ibi && (ff->b_state == 0 || ff->b_state == 1)){

         printf("Objects have boundary attributes, do you want to plot these as:\r\n"
                "    '0' -- no boundary plotting,                               \r\n"
                "    '1' -- contoured object with single colour boundary line,  \r\n"
                "    '2' -- boundary outline plot only.                       \n\n");
         scanf("%d", &pbb);

         if(pbb < 0 || pbb > 2){

            printf("****ERROR****, %d not a valid input, defaulting to no boundary plotting\n\n", pbb);
            pbb = 0;

         }

         ibi = 1;

     }


/* assign memory */

     xpos = (int *)calloc(dfil, sizeof(int));
     mem_er((xpos == NULL) ? 0 : 1, dfil * sizeof(int));

     ypos = (int *)calloc(dfil, sizeof(int));
     mem_er((ypos == NULL) ? 0 : 1, dfil * sizeof(int));

     valf = (float *)calloc(dfil, sizeof(float));
     mem_er((valf == NULL) ? 0 : 1, dfil * sizeof(float));

     dfil = 0; 

     z = (float *)calloc(gd, sizeof(float));
     mem_er((z == NULL) ? 0 : 1, gd * sizeof(float));

     x = (float *)calloc(gd, sizeof(float));
     mem_er((x == NULL) ? 0 : 1, gd * sizeof(float));

     y = (float *)calloc(gd, sizeof(float));
     mem_er((y == NULL) ? 0 : 1, gd * sizeof(float)); 

     zz = (float *)calloc(gd, sizeof(float));
     mem_er((zz == NULL) ? 0 : 1, gd * sizeof(float));      

     if(ff->tot_f_f_num){

        xfeat = (float *)calloc(ff->tot_f_f_num, sizeof(float));
        mem_er((xfeat == NULL) ? 0 : 1, (ff->tot_f_f_num) * sizeof(float));

        yfeat = (float *)calloc(ff->tot_f_f_num, sizeof(float));
        mem_er((yfeat == NULL) ? 0 : 1, (ff->tot_f_f_num) * sizeof(float));

        if(an_plt == 'y'){

           xor = (float *)calloc(ff->tot_f_f_num, sizeof(float));
           mem_er((xor == NULL) ? 0 : 1, (ff->tot_f_f_num) * sizeof(float));

           yor = (float *)calloc(ff->tot_f_f_num, sizeof(float));
           mem_er((yor == NULL) ? 0 : 1, (ff->tot_f_f_num) * sizeof(float));

           zor = (float *)calloc(2 * ff->tot_f_f_num, sizeof(float));
           mem_er((zor == NULL) ? 0 : 1, 2 * (ff->tot_f_f_num) * sizeof(float));

        }


     }

     if(pbb){

        for(i=0; i< ff->obj_num; i++){

            ob = (ff->objs)+i;

            totb += ob->bound_num;

        }

        xbb = (float *)calloc(totb, sizeof(float));
        mem_er((xbb == NULL) ? 0 : 1, totb * sizeof(float));

        ybb = (float *)calloc(totb, sizeof(float));
        mem_er((ybb == NULL) ? 0 : 1, totb * sizeof(float)); 

        vbb = (float *)calloc(totb, sizeof(float));
        mem_er((vbb == NULL) ? 0 : 1, totb * sizeof(float));

        ibb = (int *)calloc(totb, sizeof(int));
        mem_er((ibb == NULL) ? 0 : 1, totb * sizeof(int));

        ibc = 0;

        if(ff->b_state == 0){

           for(i=0; i<ff->obj_num; i++){

               ob = (ff->objs)+i;

               for(j=0; j<ob->bound_num; j++){

                   btt = (ob->bound)+j;

                   *(xbb + ibc) = *(gr->xgrid + (btt->x).ixy + x1u - 2);
                   *(ybb + ibc) = *(gr->ygrid + (btt->y).ixy + y1u - 2); 
                   *(vbb + ibc) = btt->val;
                   *(ibb + ibc) = i+1;
                   ++ibc;

               }

           }

        }

        else if(ff->b_state == 1){

           for(i=0; i<ff->obj_num; i++){

               ob = (ff->objs)+i;

               for(j=0; j<ob->bound_num; j++){

                   btt = (ob->bound)+j;

                   *(xbb + ibc) = (btt->x).xy;
                   *(ybb + ibc) = (btt->y).xy; 
                   *(vbb + ibc) = btt->val;
                   *(ibb + ibc) = i+1;
                   ++ibc;

               }

           }


        }


     }

     for(j=0; j < ff->obj_num; j++){

        ob = (ff->objs)+j;

        for(i=0; i < ob->point_num; i++){

            ptt = (ob->pt)+i;

            k=istart+i;

            xpos[k] = ptt->x;
            ypos[k] = ptt->y;
            valf[k] = ptt->val;

        }

        istart += ob->point_num;

     }

     dfil = istart;

     if(tf && xfeat){

        if(tf == 3){

           for(i=0; i < ff->obj_num; i++){

              fl = ((ff->objs)+i)->fet;

              for(j=0; j < fl->feature_num; j++){
 
                 xfeat[ic] = *(gr->xgrid + ((fl->fpt+j)->x).ixy + x1u-2);
                 yfeat[ic] = *(gr->ygrid + ((fl->fpt+j)->y).ixy + y1u-2);

                 ++ic;

              }

            }

        }

        else{

            for(i=0; i < ff->obj_num; i++){

              fl = ((ff->objs)+i)->fet;

              for(j=0; j < fl->feature_num; j++){

                 xfeat[ic] = (((fl->fpt)+j)->x).xy;
                 yfeat[ic] = (((fl->fpt)+j)->y).xy;

                 ++ic;

              }

            }



        }

        if(an_plt == 'y'){

            ic = 0;

            for(i=0; i < ff->obj_num; i++){

              fl = ((ff->objs)+i)->fet;

              for(j=0; j < fl->feature_num; j++){

/* weight by anisotropy */

                 *(xor + ic) = ((fl->fpt)+j)->ornt[0] * ((fl->fpt)+j)->r_sh;
                 *(yor + ic) = ((fl->fpt)+j)->ornt[1] * ((fl->fpt)+j)->r_sh;

                 ++ic;

              }

            }          

            if(gr->prty){
	    
	       g1.prgr = 0;
	       g1.prty = 0;

               if(gr == gr1) {
	          vp = 0;
		  g2.prgr = gr1->prgr;
		  g2.prty = gr1->prty;
	          g2.alat = gr1->alat; g2.alng = gr1->alng;
	          g2.sp1 = gr1->sp1; g2.sp2 = gr1->sp2;	       
	       }
               else if(gr == gr2){
	          vp = 1;
		  g22.prgr = gr2->prgr;
		  g22.prty = gr2->prty;
		  g22.alat = gr2->alat; g22.alng = gr2->alng;
	          g22.sp1 = gr2->sp1; g22.sp2 = gr2->sp2;
	       }
	       else {
	          vp = 0;
	          g2.prgr = 0;
		  g2.prty = gr->prty;
	       } 

               for(i=0; i < ff->tot_f_f_num; i++){

                   xx = xfeat[i];
                   yy = yfeat[i];


                   switch(gr->prgr){
                     case 0:
                        (pp->prj2)(&xx, 'x', 0);
                        (pp->prj2)(&yy, 'y', 0);
                        break;
                     case 1:
                       azimuthal(&yy, &xx, gr->alat, gr->alng, 0, gr->prty);
                       break;
		     case 2:
		       conic(&yy, &xx, gr->alat, gr->alng, gr->sp1, gr->sp2, 0, gr->prty, &ift);
		       break;
                   }


                   if(!vp)vect_proj(xx, yy, xor+i, yor+i, proj, &g1, &g2);
		   else vect_proj(xx, yy, xor+i, yor+i, proj, &g1, &g22);


               }


            } 

        }

     }


#ifdef  NOUNDERSCORE

     filterpl(gr->xgrid, gr->ygrid, z, &gr->ix, &gr->iy, xpos, ypos, valf,
               &dfil, &fruser, cm->cmi, cm->cmxg, cm->cmyg, &cm->dcm, &x1u, 
               &x2u, &y1u, &y2u, &frcount, &thresh, xfeat, yfeat, &ff->tot_f_f_num,
               x, y, zz, &iian, xan, yan, xor, yor, zor, xbb, ybb, vbb, ibb, &totb, &pbb);

#else

     filterpl_(gr->xgrid, gr->ygrid, z, &gr->ix, &gr->iy, xpos, ypos, valf,
               &dfil, &fruser, cm->cmi, cm->cmxg, cm->cmyg, &cm->dcm, &x1u, 
               &x2u, &y1u, &y2u, &frcount, &thresh, xfeat, yfeat, &ff->tot_f_f_num,
               x, y, zz, &iian, xan, yan, xor, yor, zor, xbb, ybb, vbb, ibb, &totb, &pbb);

#endif

     if(ff->tot_f_f_num){

        free(xfeat);
        free(yfeat);


     }

     free(z);
     free(x);
     free(y);
     free(zz);

     free(xpos);
     free(ypos);

     free(valf);

     if(an_plt == 'y'){

       free(xan);
       free(yan);

       free(xor);
       free(yor);

       iian = 0;

     }

     free(xbb);
     free(ybb);
     free(vbb);
     free(ibb);

     return;

}
