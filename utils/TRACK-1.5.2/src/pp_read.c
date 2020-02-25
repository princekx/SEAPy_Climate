#include <Stdio.h>
#include <stdlib.h>
#include "pp.h"
#include "grid.h"

int pp_header(GRID * , FILE * , int * , int * , int , int );

extern PP *pph;
extern void *databuf;
extern float *abuf;

int pp_read(FILE *fpp, int pr)
{

  int i, j, it, itt;
  long int ss;

  float *ddat=(float *)databuf;

  if(!pp_header(NULL, fpp, NULL, NULL, 0, pr)) return 0;

  fread(&ss,4,1,fpp);
 
  fread(ddat, pph->lbnpt * pph->lbrow * sizeof(float),1,fpp);

  fread(&ss,4,1,fpp);


  if(pr > 0){

     if(pph->bdy < 0){

       for(i=0; i< pph->lbrow; i++){

           it = (pph->lbrow - 1 - i) * pph->lbnpt;
           itt = i * pph->lbnpt;

           for(j=0; j< pph->lbnpt; j++){

               *(abuf + it + j) = *(ddat + itt + j);

           }

       }

     }

     else{

        for(i=0; i < pph->lbrow; i++){

           it = i * pph->lbnpt;

           for(j=0; j < pph->lbnpt; j++){

             *(abuf + it + j) = *(ddat +it + j);

           }

        } 

     }

  }

  else if (pr < 0){

    for(i=0; i < pph->lbrow * pph->lbnpt; i++) *(abuf + i) = *(ddat + i);


  }
   
  return 1;

}

