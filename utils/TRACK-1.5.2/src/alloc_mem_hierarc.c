#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include <errno.h>
#include "st_im.h"

/* fuction to allocate/de-allocate storage for the data hierarchys */

void alloc_mem_hierarc(struct image **hierarchy, int perf, int level, int xdim, int ydim)

{

   int i, j, k, n, xtemp, ytemp, plus, mem_flag;
   struct image *pp;

   xtemp=xdim;
   ytemp=ydim;
   n=xtemp*ytemp;
   mem_flag = 0;

/* allocate storage for hierarchy pointers */

   if(perf == 1){

      for(i=0; i<level; i++){

          if(i > 0 ) {

               hierarchy[i]=(struct image *)calloc(n, sizeof(struct image));
               mem_er((hierarchy[i] == NULL) ? 0 : 1, n * sizeof(struct image));

           }

/* initialize pointers for dummy nodes introduced for rectanguler regions */

           if(mem_flag == 1) {

              for(j=0; j < ytemp; j++){

                  pp=hierarchy[i]+2*j+1;

                  pp->pval=0;
                  pp->nw=pp->ne=pp->sw=pp->se=NULL;

               }
           }
           if(mem_flag == 2){

               for(j=0; j < xtemp; j++){

                   pp=hierarchy[i]+xtemp+j;

                   pp->pval=0;
                   pp->nw=pp->ne=pp->sw=pp->se=NULL;

                }
            }

/* allocate storage for the equivelence classes  */

          for(j=0; j< ytemp; j++) {

              for(k=0; k< xtemp; k++) {

                  plus=j*xtemp+k;
                  pp=hierarchy[i]+plus;

                  pp->label = pp->stor = (struct eq_class *)malloc_initl(sizeof(struct eq_class));
                  mem_er((pp->label == NULL) ? 0 : 1, sizeof(struct eq_class));

              }
          } 

          xtemp=xtemp/2;
          ytemp=ytemp/2;
          n=xtemp*ytemp;
          mem_flag=0;
          if(xtemp == 1 && n != 1) {
            xtemp=2;
            mem_flag = 1;
            n=2*ytemp;
          }
          if(ytemp == 1 && n != 1) {
                  ytemp=2;
                  mem_flag = 2;
                  n=2*xtemp;
          }

      }

    }

    else if(perf == 0) {

      for(i=0; i<level; i++){

/* free storage for the equivelence classes  */

          for(j=0; j<ytemp; j++) {

              for(k=0; k<xtemp; k++) {

                  plus=j*xtemp+k;
                  pp=hierarchy[i]+plus;

/*                  if(pp->label && pp->label->farther) 
                      pp->label->farther = NULL; */


                  pp->label = NULL;
                  pp->mother = NULL;
                  pp->nw = pp->ne = pp->sw = pp->se = NULL;

                  free(pp->stor);

              }

          }

          xtemp=xtemp/2;
          ytemp=ytemp/2;
          n=xtemp*ytemp;
          if(xtemp == 1 && n != 1) xtemp=2;
          if(ytemp == 1 && n != 1) ytemp=2;


          if(i > 0 ) {

            free(hierarchy[i]);
            hierarchy[i]=NULL;

          }

       }


    }

    return;

}
