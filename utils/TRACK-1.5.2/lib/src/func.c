#include <Stdio.h>
#include <Math.h>
#include "bisp.h"

/* this routine evaluates a function and its first derivatives at a
   given point for use in the function optimization routine.         */

#ifdef NOUNDERSCORE

void bisp(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#else

void bisp_(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#endif

extern struct sp_dat *cokn;

double func(double *vv, double *dv, int dimv, int type, int opt)

{

   int i;

   double val=0.;

/* if opt=1 compute derivatives */


#ifdef NOUNDERSCORE
  
   bisp(&val, dv, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof,
         &vv[0], &vv[1], &opt);

#else

   bisp_(&val, dv, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof,
         &vv[0], &vv[1], &opt); 

#endif

   if(type){

      for(i=0; i < dimv; i++) dv[i] *= -1.0;
      val *= -1.0;

   }

   return val;

}
   
