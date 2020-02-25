#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "mem_er.h"


#include "grid.h"

#define  MAXCHAR   100

/* function to read a standard header i.e. 

   xdim  ydim  frames

   xgrid ......->
    .
    .
    .
    |
   \ /

   ygrid ...... ->
    .
    .
    .
    |
   \ /
                                               */


int std_x, std_y;

extern int init;
extern float period;
extern float *abuf;
extern int nhp, shp;

extern FILE *finit;

void grid_trans(GRID * , float * , int * , int * );

void std_read_header(GRID *gr, int *frnum, FILE *fdatin, int *tl, int *gof)

{

    int i;

    char text[MAXCHAR];

    float *xgtmp=NULL;

    fscanf(fdatin, "%d %d %d", &std_x, &std_y, frnum); /* read grid dimensions 
                                                     and no. of frames */

    gr->ix = std_x;
    gr->iy = std_y;

/* assign memory for grid */

    gr->xgrid = (float * )calloc(gr->ix + 1, sizeof(float));
    mem_er((gr->xgrid == NULL) ? 0 :1, (gr->ix + 1) * sizeof(float));

    xgtmp = (float * )calloc(gr->ix, sizeof(float));
    mem_er((xgtmp == NULL) ? 0 :1, gr->ix * sizeof(float));

    gr->ygrid = (float * )calloc(gr->iy, sizeof(float));
    mem_er((gr->ygrid == NULL) ? 0 : 1, gr->iy * sizeof(float));

/* assign memory for reading field data */

    abuf = (float *)calloc(std_x * std_y, sizeof(float));
    mem_er((abuf == NULL) ? 0 : 1, std_x * std_y * sizeof(float));


/* read X grid data */

    for(i=0; i < gr->ix; i++) fscanf(fdatin, "%f", (xgtmp+i));

/* read Y grid data */

    for(i= 0; i< gr->iy; i++)

         fscanf(fdatin, "%f", ((gr->ygrid)+i));

    fgets(text, MAXCHAR, fdatin);

    grid_trans(gr, xgtmp, tl, gof);

    free(xgtmp);

    return;

}

