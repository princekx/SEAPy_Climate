#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "grid.h"
#include "mem_er.h"


/* function to perform grid transformations, translations and inversions to
   a standard form.                                                         */

int neqt='y', ieq=0;

extern int init;
extern float period;
extern int nhp, shp;
extern int std_x, std_y;

extern FILE *finit;

void grid_trans(GRID *gr, float *xgtmp, int *tl, int *gof)
{

    int i, k, gf=0;
    int itr='n';

    float hg0, hg1;
    float add, ftmp; 
    float foff=0.0;

    printf("the current grid dimensions are %d * %d \n\n", gr->ix, gr->iy);

/* is a translation of the grid required */

    printf("do you want to translate the grid 'y' or 'n'?\n\n");


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

      printf("what is the longitude offset required No. of x grid points?\n");

      if(init) fscanf(finit, "%d", gof);
      else{
        scanf("%d", gof);
        fprintf(finit, "%d\n", *gof);
      }
 
    }

    if(fabs(*xgtmp) < TOLGRID) *xgtmp = 0.0;
    if(fabs(*(xgtmp + gr->ix - 1) - period) < TOLGRID) *(xgtmp + gr->ix - 1) = period;

    gr->gcen = 0;

    if(*xgtmp > 0. && *(xgtmp + gr->ix - 1) < 0.){

       for(i=0; i < gr->ix; i++) if(*(xgtmp + i) < 0) *(xgtmp + i) += period;

    }

    if(*xgtmp < 0. || *(xgtmp + gr->ix - 1) > period){

        printf("****WARNING****, grid is not contained in (0, 360) longitude, this maybe because data     \r\n"
               "                 is not defined on the sphere. Do you want to try and correct, 'y' or 'n'.\n\n");

        if(init){
           fscanf(finit,"\n");
           itr = getc(finit);
        }
        else{
           scanf("\n");
           itr = getchar();
           fprintf(finit, "%c\n", itr);
        }
        if(itr == 'y'){

           if((*xgtmp < 0. && *(xgtmp + gr->ix - 1) > period) ||
               (*xgtmp < -period && *(xgtmp + gr->ix - 1) > 0.)){
              printf("****ERROR****, grid cannot be contained in (0, 360).\n\n");
              exit(1);
           }
  
           printf("For correct grid do you want data transformation '1', or simple translation '2'.\r\n"
                  "Transformation only available for data in (-180, 180) otherwize translation.    \n\n");

           if(init) fscanf(finit, "%d", &gr->gcen);
           else {
              scanf("%d", &gr->gcen);
              fprintf(finit, "%d\n", gr->gcen);
           }

           switch(gr->gcen){
              case 1:
                 printf("****INFORMATION****, using grid and data transformation.\n\n");
                 break;
              case 2:
                 printf("****INFORMATION****, using grid translation only.\n\n");
                 break;
              default:
                 printf("****INFORMATION****, chosen option not valid using transformation as default.\n\n");
                 break;
           }

           if(!(gr->gcen == 1 || gr->gcen == 2)) gr->gcen = 1;
           if(gr->gcen == 1 && *(xgtmp + gr->ix - 1) > period) gr->gcen = 2;

              
           if(gr->gcen == 1) {

              gr->nneg = 0;

              for(i=0; i < gr->ix; i++) if(*(xgtmp + i) < 0) ++(gr->nneg);

              gr->wrk = (float * )calloc(gr->ix, sizeof(float));
              mem_er((gr->wrk == NULL) ? 0 :1, gr->ix * sizeof(float));

              for(i=gr->nneg; i < gr->ix; i++) *(gr->wrk + i - gr->nneg) = *(xgtmp + i);
              for(i=0; i<gr->nneg; i++) *(gr->wrk + gr->ix - gr->nneg + i) = *(xgtmp + i) + period ;

              for(i=0; i <gr->ix; i++)*(xgtmp + i) = *(gr->wrk + i);
           
           }

           else {

              if(*xgtmp < 0.){
                foff = fabs(*xgtmp);
                gr->f_off = foff;
                for(i=0; i < gr->ix; i++) *(xgtmp + i) += foff; 
                printf("****INFORMATION****, grid translated by %f\n\n", foff);
              }
              else if(*(xgtmp + gr->ix - 1) > period){
                foff = *(xgtmp + gr->ix - 1) - period;
                gr->f_off = -foff;
                for(i=0; i < gr->ix; i++) *(xgtmp + i) -= foff;
                printf("****INFORMATION****, grid translated by %f\n\n", -foff);
              }

              if(*xgtmp < 0. || *(xgtmp + gr->ix - 1) > period){
                 printf("****ERROR****, grid not contained in (0, 360) longitude. Un-able to correct.\n\n");
                 exit(1);
              }

           }
        }

        else {
           printf("****INFORMATION****, use the Euclidean distance measure if using a projection \r\n"
	          "                     for the tracking, i.e. not using the spherical tracking. \n\n");
        }

    }



/* test for X uniformity */

   for(i=1; i<gr->ix - 1; i++){

       if(*(xgtmp + i + 1) - 2.0 * *(xgtmp + i) + *(xgtmp + i - 1) > TOLGRID){

          printf("****WARNING****, grid is not uniform in the X direction.\r\n"
                 "                 no grid translation possible.          \n\n");
          *gof = 0;
          *tl = 'n';
          gr->igtyp = 3;
          gr->ixfrm = 1;
          break;

       }

   }

/* test for global/periodic grid */

    if(gr->igtyp != 3 && (fabs(*(xgtmp)) > TOLGRID || fabs(*(xgtmp + gr->ix - 1) - period) > TOLGRID)){

      printf("****WARNING****, grid is not periodic.\n\n");

      gr->igtyp = 1;
      add = period;
 
    }

    else add = period + *(xgtmp+1) - *xgtmp;

    hg0 = *(xgtmp+1) - *xgtmp;
    if(fabs(0.5 * hg0 - *xgtmp) < TOLGRID) hg0 = *xgtmp;

    hg1 = period - *(xgtmp + gr->ix - 1);
    if(hg1 < TOLGRID) hg1 = period - *(xgtmp + gr->ix - 2);

    if(gr->igtyp != 3 && fabs(hg1 - hg0) > TOLGRID){

      printf("****WARNING****, grid is not global.\n\n");
      if(*tl == 'y'){

          printf("****WARNING****, grid translation not possible for non-global data.\n\n");
          *tl = 'n';
          *gof = 0;

       }
       gr->igtyp = 2;

     }


     if(*gof == 0){

       for(i= 0; i< gr->ix; i++) *((gr->xgrid)+i) = *(xgtmp + i);


     }

     else if(*gof > 0 && *tl == 'y'){

       for(i= *gof; i< gr->ix; i++){

            *((gr->xgrid)+i) = *(xgtmp + i - *gof);

       }


       for(i=0; i< *gof; i++){

          *((gr->xgrid)+i) = *(xgtmp + gr->ix - *gof + i);
          *((gr->xgrid)+i) -= add;
       }

     }

     else if(*gof < 0 && *tl == 'y') {

       gf = -1 * (*gof);

       for(i=0; i < gf; i++){

           k = gr->ix - gf + i;

           *((gr->xgrid)+k) = *(xgtmp + i);
           *((gr->xgrid)+k) += add;

       }


       for(i=gf; i < gr->ix; i++){

             *((gr->xgrid) + i - gf) = *(xgtmp + i);

       }

     }


/* Y-grid transforms */

/* test for Y uniformity */

    for(i=1; i<gr->iy - 1; i++){

        if(*(gr->ygrid + i + 1) - 2.0 * *(gr->ygrid + i) + *(gr->ygrid + i - 1) > TOLGRID){

          printf("****WARNING****, grid is not uniform in the Y direction.\n\n");

          gr->iyfrm = 1;
          break;

        }

    }

    if(*(gr->ygrid) > *(gr->ygrid + 1)){

       printf("****WARNING****, grid is inverted to expected ordering\r\n"
              "                 grid will be corrected.              \n\n");

       gr->h_inv = 1;

       for(i=0; i < gr->iy/2; i++){

          ftmp = *(gr->ygrid + i);
          *(gr->ygrid + i) = *(gr->ygrid + gr->iy - i - 1);
          *(gr->ygrid + gr->iy - i - 1) = ftmp;

       }

       if(!(gr->wrk)){
           gr->wrk = (float * )calloc(gr->ix, sizeof(float));
           mem_er((gr->wrk == NULL) ? 0 :1, gr->ix * sizeof(float));
       }

    }
    else gr->h_inv = 0;


/* test values are in the range (-90, 90) */

    if(*(gr->ygrid) + 90.0 < -TOLPOLE || *(gr->ygrid + gr->iy - 1) - 90.0 > TOLPOLE){
       printf("****WARNING****, latitude outside of range (-90, 90), data may not be defined on a sphere.\r\n"
              "                 use Euclidean distance measure.                                          \n\n");
    }

/* test for equator */

    for(i=0; i < gr->iy; i++){
        if(fabs(*(gr->ygrid + i)) < TOLGRID){
           printf("****WARNING****, data contains the equator, do you want to retain it, 'y' or 'n'\n\n");
           if(init){
              fscanf(finit,"\n");
              neqt = getc(finit);
           }
           else{
              scanf("\n");
              neqt = getchar();
              fprintf(finit, "%c\n", neqt);   
           }
           if(neqt == 'n') ieq = i;
           break;
        }
    }

    if(neqt == 'n'){
       for(i=ieq; i < gr->iy - 1; i++) *(gr->ygrid + i) =  *(gr->ygrid + i + 1);
       gr->iy -= 1;
    }

    if(fabs(*(gr->ygrid) + 90.0) <= TOLPOLE){

       printf("***WARNING***, data contains a SH pole, do you want to retain it, 'y' or 'n'\n\n");

       *(gr->ygrid) = -90.0;
       
       if(init){
          fscanf(finit,"\n");
          shp = getc(finit);
       }
       else{
          scanf("\n");
          shp = getchar();
          fprintf(finit, "%c\n", shp);   
       }

       if(shp == 'n'){ ++(gr->ygrid); gr->iy -= 1;}

    }

    if(fabs(*(gr->ygrid + gr->iy - 1) - 90.0) <= TOLPOLE){

       printf("***WARNING***, data contains a NH pole, do you want to retain it, 'y' or 'n'\n\n");

       *(gr->ygrid + gr->iy - 1) = 90.0;       

       if(init){
          fscanf(finit,"\n");
          nhp = getc(finit);
       }
       else{
          scanf("\n");
          nhp = getchar();
          fprintf(finit, "%c\n", nhp);   
       }

       if(nhp == 'n') {gr->iy -= 1;}

    }

    if(shp == 'n' || nhp == 'n' || neqt == 'n') 

       printf("\n the NEW grid dimensions are %d * %d \n\n", gr->ix, gr->iy);


    return;

}
