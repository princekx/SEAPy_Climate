#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"
#include "st_obj.h"

/* this function constructs the hierarchical data structure 
   for each frame and uses the structure to perform the segmentation, 
   returning a pointer to an array of structures of objects.             */



void alloc_mem_hierarc(struct image ** , int , int , int , int );
int assign_node(int , int , int , int );
int assign_label(struct image *, int , int );
void form_objects(struct image * , struct frame_objs * , int , int , int );
void hprint(struct image ** , int , int , int );
void update(struct image *);


void hierarc_segment(struct image **hierarchy, struct image *bin, int level, struct frame_objs *ff, int xxdim, int yydim, int ch)

{

    int i, j, k, n, xtemp, ytemp, xdim, ydim;
    int perf, plus, nclass=0, flag;

    struct image *p1, *p2, *p3;

    xdim=xxdim;
    ydim=yydim;

    hierarchy[0]=bin;


    alloc_mem_hierarc(hierarchy, perf=1, level, xxdim, yydim);     


/*    for(j=0; j<ydim; j++){
         for(i=0; i<xdim; i++){
             plus=j*xdim+i;
             printf("%d ", (hierarchy[0]+plus)->pval);
         }
         printf("\n");
     }
     printf("\n"); */


     for(i=1; i< level; i++){

         xtemp = xdim;
         ytemp = ydim;
         xdim=xdim/2;
         ydim=ydim/2;

         n=xdim*ydim;
         flag=0;
         if(xdim == 1 && n != 1) {
            flag=1;
            xdim=2;
         }
         if(ydim == 1 && n != 1) {
            flag=2;
            ydim=2;
         }
         p1=hierarchy[i-1];
         p2=p1+xtemp;

         for(j=0; j< ((flag==2) ? 1 : ydim); j++){

             for(k=0; k< ((flag==1) ? 1 : xdim); k++){

                 plus=j*xdim+k;
                 p3=hierarchy[i]+plus;

                 p3->pval= assign_node(p1->pval, (p1+1)->pval, p2->pval, (p2+1)->pval);

                 p1->mother = (p1+1)->mother=p2->mother=(p2+1)->mother=p3;

                 if(p3->pval == -1){
                    p3->nw = p1;
                    p3->ne = (p1+1);
                    p3->sw = p2;
                    p3->se = (p2+1);
                 }

                 else
                    p3->nw=p3->ne=p3->sw=p3->se=NULL;

/*                 printf("%2d ", p3->pval); */
                  
                 p1+=2;
                 p2+=2;


             }

             p1+=xtemp;
             p2+=xtemp;
/*             printf("\n"); */
         }
/*         printf("\n\n"); */


     }

     assign_label(hierarchy[level-1], nclass, ch);

/*     hprint(hierarchy, level, xxdim, yydim); */

     update(hierarchy[level-1]);

/*     hprint(hierarchy, level, xxdim, yydim); */

     form_objects(hierarchy[level-1], ff, level-1, 0, 0);

     alloc_mem_hierarc(hierarchy, perf=0, level, xxdim, yydim);

     return;

}
