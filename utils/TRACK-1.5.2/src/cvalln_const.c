#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "mem_er.h"
#include "m_values.h"

/* calculate cross-validation function for a constant kernal function*/

double const_density(struct dpt * , struct dpt * , double );

double cvalln_const(struct dpt *dt, int dtn, double sm, double *plt, float *wght, int im, int ht, int ty, int ms)

{

   int i, j;

   int sw=0;

   double cn, cons=0., cnn, cd;
   double fd, ftd, dif;
   double *d1=NULL, *d2=NULL, *d3=NULL;
   double llh=0.;
   double nn;
   double tsm=0., ts;
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

   if(ht < 3){

      if(sm < SMMIN) {tsm = -1.0; cons = 1.0;}

      else{

         tsm = (1.0 + sm) / sm;

         cn *= tsm;

         cd = tsm - 1.0;

         cons = cn / cd;

      }

   }

   else if(ht == 3) cons = cn;

   for(i=0; i < dtn; i++) {

       if(wght) nn =0.0;

       dt1 = dt + i;

       fd = (sw) ? 1.0e-5 : 1.0e-20;

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
             

             fd += (ftd = wt * const_density(dt1, dt2, tsm));

             if(!sw){
                if(d1) *d1 += ftd * dt2->sdt;
                else {*d2 += ftd * dt2->vec[0]; *d3 += ftd * dt2->vec[1];}
             }

          }


       }

       else if (ht == 3){


          for(j=0; j < dtn; j++){
 
             if(!ms && i == j) continue;

             ts = *(plt + j) * sm;

             if(ts < SMMIN) {tsm = -1.0; cnn = 1.0;}

             else {

                tsm = (1.0 + ts) / ts;

                cd = tsm - 1.0;

                cnn = tsm / cd;

             }

             dt2 = dt + j;

             if(wght){wt = *(wght + j); nn += wt;}
             else wt = 1.0;

             fd += (ftd = wt * cnn * const_density(dt1, dt2, tsm));

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
           else fd = 1.0e-20;
           
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
