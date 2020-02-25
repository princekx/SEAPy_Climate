#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem_er.h"

#include "grid.h"

#define  GL      360.00

/* function to read a utf 1.3 header i.e. 

   xdim  ydim  frames

   xgrid ......->
    .
    .
    .
    |
   \ /

   ygrid ...... <-
    .
    .
    .
   / \
    |

and perform transformations.                           */


extern int init;
extern float period;
extern float *abuf;

extern FILE *finit;

extern int eqsw;
extern int whemi;

void utf3_read_header(float ival, GRID *gr, int *frnum, FILE *fdatin, int *tl, int *gof)

{

    int i, k;
    int in[5];
    int gf=0;
    int dim=0;

    float db, de;

    char ch='0';


/* read data dimensions */

    i = 0;

    fscanf(fdatin, "%d %d %d %d", &in[0], &in[1], &in[2], &in[3]);
    fscanf(fdatin, "%*d %d", &eqsw);

 
    for(i=0; i < 9; i++) fscanf(fdatin, "%*d");

    gr->ix = in[0] + 1;
    gr->iy = in[1] * in[3];

    printf("the grid dimensions are %d * %d \n", gr->ix, gr->iy);

/* is a translation of the grid required */

    printf("do you want to translate the grid 'y' or 'n'?\n\n");
    printf("***WARNING***, only for images periodic in x direction\n");

    if(init){
       fscanf(finit,"\n");
       *tl = getc(finit);
    }
    else{
       scanf("\n");
       *tl = getchar();
       fprintf(finit, "%c\n", *tl);
    }
    if(*tl == 'y'){

      printf("what is the lattitude offset required No. of x grid points?\n");

      if(init) fscanf(finit, "%d", gof);
      else{
         scanf("%d", gof);
         fprintf(finit, "%d\n", *gof);
      }

    }

/* assign memory for grid */

    gr->xgrid = (float * )calloc(gr->ix, sizeof(float));
    mem_er((gr->xgrid == NULL) ? 0 :1, gr->ix * sizeof(float));

    gr->ygrid = (float * )calloc(gr->iy, sizeof(float));
    mem_er((gr->ygrid == NULL) ? 0 : 1, gr->iy * sizeof(float));

/* assign memory for reading field data */

    dim = gr->ix * ((eqsw) ? gr->iy + 1 : gr->iy);

    abuf = (float *)calloc(dim, sizeof(float));
    mem_er((abuf == NULL) ? 0 : 1, dim * sizeof(float));

/* read X grid data */

    if(*gof == 0){

       for(i= *gof; i< (gr->ix)-1; i++){

          fscanf(fdatin, "%f", ((gr->xgrid)+i));

       } 

       *((gr->xgrid) + gr->ix - 1) = GL;

    }

    else if(*gof > 0){

       for(i= *gof; i< gr->ix; i++){

          fscanf(fdatin, "%f", ((gr->xgrid)+i)); 

       }

       for(i=1; i< *gof; i++){

          fscanf(fdatin, "%f", ((gr->xgrid)+i));
          *((gr->xgrid)+i) -= period;

       }

       *(gr->xgrid) = *((gr->xgrid) + gr->ix - 1) - period;

    }

    else {

       gf = -1 * (*gof);

       for(i=0; i < gf; i++){

           k = gr->ix - gf + i - 1;

           fscanf(fdatin, "%f", ((gr->xgrid)+k));
           *((gr->xgrid)+k) += period;

       }

       for(i=gf; i < (gr->ix)-1; i++){

           fscanf(fdatin, "%f", ((gr->xgrid)+i - gf));

       }

       *((gr->xgrid) + gr->ix - 1) = *(gr->xgrid) + period;

     } 

/* read Y grid data */

    for(i= 0; i< in[1]; i++)

         fscanf(fdatin, "%f", ((gr->ygrid)+ gr->iy - i - 1));


    if(in[3] == 2){
    
       for(i= 0; i< in[1]; i++)

           fscanf(fdatin, "%f", ((gr->ygrid)+ gr->iy - in[1] - i - 1));

    }

    if(in[3] == 1){

      if(*gr->ygrid > 0.) whemi = 1;
      else whemi = 2;

    }

    for(i=0; i < in[2]; i++)

        fscanf(fdatin, "%*f");

    fscanf(fdatin, "%*s %*f %*s %f %*s %f", &db, &de);

    ch = '0';
    while(ch != '\n') fscanf(fdatin, "%c", &ch);

    *frnum = (int)((de - db) / ival) + 1;

    ch = '0';
    while(ch != '\n') fscanf(fdatin, "%c", &ch);
    ch = '0';    
    while(ch != '\n') fscanf(fdatin, "%c", &ch);

    return;

}
