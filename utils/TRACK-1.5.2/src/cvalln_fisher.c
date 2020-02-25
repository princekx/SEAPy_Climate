#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "mem_er.h"
#include "m_values.h"

/* calculate cross-validation function for a fisher kernal function */

double fisher_density(struct dpt * , struct dpt * , double ,double ,int );

double cvalln_fisher(struct dpt *dt, int dtn, double sm, double *plt, float *wght, int im, int ht, int ty, int ms)

{

   int i, j;
   int nf=0;
   int sw=0;

   double cn, cons=0., cnn;
   double fd, ftd, dif;
   double *d1=NULL, *d2=NULL, *d3=NULL;
   double llh=0.;
   double nn, ee, ee1, lgtol=0., lgcn;
   double tsm=0.;
   double wt=0.0;

   struct dpt *dt1, *dt2;

   if(ht < 0 || ht > 3){

      printf("***error***, non-valid switch type in %s \n\n", __FILE__);
      return 0.;

   }

   if(im > 0) sw = 1;

   if(!sw) {

      if(im == 0) {

         d1 = (double * )malloc_initl(sizeof(double));
         mem_er((d1 == NULL) ? 0 : 1, sizeof(double));

      }

      else if(im == -1){

         d2 = (double * )malloc_initl(sizeof(double));
         mem_er((d2 == NULL) ? 0 : 1, sizeof(double));

         d3 = (double * )malloc_initl(sizeof(double));
         mem_er((d3 == NULL) ? 0 : 1, sizeof(double));

      }


   }

   if(!wght) {
      nn = (ms) ? (float)dtn : (float) (dtn - 1);
      cn = 1.0 / (2.0 * FPI * nn);
   }

   else {
      nn = 0.0;
      cn = 1.0 / (2.0 * FPI);
   }



   lgcn = EXPON + log(cn);

   if(ht < 3){

      tsm = sm;

      if(tsm < SMMIN) cons = 0.5 * cn;

      else {

         if(tsm > DBAND) nf = 1;

         ee = exp(tsm);
         ee1 = 1.0 / ee;

         cn *= tsm;

         cons = (nf) ? cn :  cn / (ee - ee1);

      }

      lgtol = -(EXPON + log(cons)) / tsm;

   }

   else if(ht == 3) cons = cn;

   for(i=0; i < dtn; i++) {

       if(wght) nn =0.0;

       dt1 = dt + i;

       fd = 1.0e-100;

       if(!sw){

          if(d1) *d1 = 0.0;
          else *d2 = *d3 = 0.0;

       }

       if(ht < 3) {

          for(j=0; j < dtn; j++){

              if(!ms && i == j) continue;

              dt2 = dt + j;

              if(wght){wt = *(wght + j); nn += wt;}
              else wt = 1.0;

              fd += (ftd = wt * fisher_density(dt1, dt2, tsm, lgtol, nf));

              if(!sw){
                 if(d1) *d1 += ftd * dt2->sdt;
                 else {*d2 += ftd * dt2->vec[0]; *d3 += ftd * dt2->vec[1];}
              }

          }

       }

       else if (ht == 3){

          for(j=0; j < dtn; j++){

              nf = 0;

              if(!ms && i == j) continue;

              tsm = *(plt + j) * sm;

              if(tsm < SMMIN) cnn = 0.5 * cn;

              else {

                 if(tsm > DBAND) {    /* change evaluation of exponetial */

                   nf = 1;
                   cnn = tsm;

                 }

                 else{

                   ee = exp(tsm);
                   ee1 = 1.0 / ee;
                   cnn = tsm / (ee - ee1);

                 }

              }

              lgtol = -(log(cnn) + lgcn) / tsm;

              dt2 = dt + j;

              if(wght){wt = *(wght + j); nn += wt;}
              else wt = 1.0;

              fd += (ftd = wt * cnn * fisher_density(dt1, dt2, tsm, lgtol, nf));

              if(!sw){
                 if(d1) *d1 += ftd * dt2->sdt;
                 else {*d2 += ftd * dt2->vec[0]; *d3 += ftd * dt2->vec[1];}
              }

          }

       }

       if(!sw){

          if(d1) {
              *d1 = (fd > DTOL) ? *d1 / fd : 0.0;
              dif = dt1->sdt - *d1;
              if(ms) dt1->dsdt = *d1;
              llh -= (dif*dif);
          }
          else{ 
             if(fd > DTOL) {*d2 /= fd; *d3 /= fd;}
             else { *d2 = 0.0; *d3 = 0.0;}
             dif = dt1->vec[0] - *d2;
             llh -= (dif*dif);
             dif = dt1->vec[1] - *d3;
             if(ms) {dt1->dvec[0] = *d2; dt1->dvec[1] = *d3;}
             llh -= (dif*dif);
          }

       }

       else {

           if(!wght) fd *= cons;
           else if(nn > TOLWT) fd *= cons / nn;
           else fd = 1.0e-100; 

           llh += log(fd);

           if(ty) *(plt + i) = fd;

       }

   }

   if(!sw){
      if(d1) free(d1);
      if(d2){ free(d2); free(d3);}
   }

   return llh;

}
