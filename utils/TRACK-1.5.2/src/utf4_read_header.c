#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "mem_er.h"

#include "grid.h"

#define  MAXCHR  150

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

int nhp='y', shp='y';
int i_utf4x, i_utf4y; 

void grid_trans(GRID * , float * , int * , int * );

void utf4_read_header(float ival, GRID *gr, int *frnum, FILE *fdatin, int *tl, int *gof)

{

    int i=0;
    int in[5];
    int yeq='n';
    int iyh;

    float db, de;
    float *gg=NULL;
    float *xgtmp=NULL;

    char aa[MAXCHR];

    *frnum = MAXFRM;

/* read data dimensions */

    fgets(aa, MAXCHR, fdatin);

    sscanf(aa, "%d%d%d%d%*d%d%*d%*d%*d%*d%d", &in[0], &in[1], &in[2], &in[3], &eqsw, &in[4]);

    if(eqsw){
       printf("***WARNING***, this data contains values at the equator\r\n"
              "               do you wish to retain the equator?      \n\n");
       if(init){
          fscanf(finit,"\n");
          yeq = getc(finit);
       }
       else{
          scanf("\n");
          yeq = getchar();
          fprintf(finit, "%c\n", yeq);   
       }

    }

    i_utf4x = in[0];
    gr->ix = in[0];
    gr->iy = in[1] * in[3];

    if(in[3] == 2) gr->iy -= eqsw;

    i_utf4y = gr->iy;

    if(eqsw){
      if(yeq == 'y')
         printf("the grid dimensions are %d * %d \n", gr->ix, gr->iy);
      else
         printf("the grid dimensions are %d * %d \n", gr->ix, gr->iy - 1);
    }
    else
      printf("the grid dimensions are %d * %d \n", gr->ix, gr->iy);


/* assign memory for reading field data */

    abuf = (float *)calloc(gr->ix * gr->iy, sizeof(float));
    mem_er((abuf == NULL) ? 0 : 1, gr->ix * gr->iy * sizeof(float));

    if(eqsw){

      iyh = (gr->iy - 1) / in[3];
      if(yeq == 'n') gr->iy -= 1;

    }

    else iyh = gr->iy / in[3];


/* assign memory for grid */

    gr->xgrid = (float * )calloc(gr->ix + 1, sizeof(float));
    mem_er((gr->xgrid == NULL) ? 0 :1, (gr->ix + 1) * sizeof(float));

    xgtmp = (float * )calloc(gr->ix, sizeof(float));
    mem_er((xgtmp == NULL) ? 0 :1, gr->ix * sizeof(float));

    gr->ygrid = (float * )calloc(gr->iy, sizeof(float));
    mem_er((gr->ygrid == NULL) ? 0 : 1, gr->iy * sizeof(float));

/* read X grid data */

    for(i=0; i < gr->ix; i++) fscanf(fdatin, "%f", (xgtmp+i));




/* read Y grid data */

    gg = (gr->ygrid)+ gr->iy - 1;

    for(i= 0; i< iyh; i++) fscanf(fdatin, "%f", (gg--));
    
    if(eqsw){
       if(yeq == 'y') {fscanf(fdatin, "%f", (gg--)); eqsw = 0;}
       else fscanf(fdatin, "%*f ");
    }

    if(in[3] == 2){

       for(i= 0; i< iyh; i++) fscanf(fdatin, "%f", (gg--));

    }

    if(in[3] == 1 && eqsw && yeq == 'n'){

      if(*gr->ygrid > 0.) whemi = 1;
      else whemi = 2;

    }

    grid_trans(gr, xgtmp, tl, gof);

    for(i=0; i < in[2]; i++) fscanf(fdatin, "%*f");

    for(i=0; i < in[4]; i++) fscanf(fdatin, "%*f");

    fgets(aa, MAXCHR, fdatin);

    fgets(aa, MAXCHR, fdatin);

    if(strstr(aa, "DAYS") || strstr(aa, "days") || strstr(aa, "Days")){

       if(sscanf(aa, "%*s %*f %*s %f %*s %f", &db, &de) == 2) {

          *frnum = (int)((de - db) / ival) + 1;

          printf("****INFORMATION****, No. of frames in this file is -- %d\n", *frnum);

       }

    }

    else {
       printf("****INFORMATION****, default No. of frames for this file is -- %d\n\n", *frnum);
 
    }

    fgets(aa, MAXCHR, fdatin);

    fgets(aa, MAXCHR, fdatin);

    free(xgtmp);

    return;

}
