#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "mem_er.h"
#include "m_values.h"
#include "p_vecs.h"

/* calculate cross-validation function for a linear kernal function*/

double non_iso_density(struct dpt * , struct dpt *, struct cvecs * , double , double ,int , int * );
double power_norm(double , double , int );

double cvalln_non_iso(struct dpt *dt, int dtn, double sm, double *plt, float *wght, int im, int ht, int ty, int ms, struct cvecs *pcom, int nni)

{

   int i, j;
   int ireg=0;

   int sw=0;

   double cn, cons=0., cnn;
   double fd, ftd, dif;
   double *d1=NULL, *d2=NULL, *d3=NULL;
   double llh=0.;
   double nn;
   double tsm=0., ts;
   double beta;
   double ddd=0.0;
   double wt=0.0;

   struct dpt *dt1, *dt2;
   struct cvecs *pp=NULL;


   if(ht != 3){

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


   if(!wght){

      nn = (ms) ? (float)dtn : (float) (dtn - 1);

      cn = 1.0 / (FPI2 * nn);

   }

   else {
      nn = 0.0;
      cn = 1.0 / FPI2;
   }

   cons = cn;



   for(i=0; i < dtn; i++) {

       if(wght) nn =0.0;

       dt1 = dt + i;

       fd = (sw) ? 1.0e-5 : 1.0e-20;

       if(!sw){

          if(d1) *d1 = 0.0;
          else *d2 = *d3 = 0.0;

       }



       for(j=0; j < dtn; j++){

           if(!ms && i == j) continue;

           ts = *(plt + j) * sm;

           tsm = (1.0 + ts) / ts;

           dt2 = dt + j;

           pp = pcom + j;

           beta = pp->beta * tsm;

           ddd = non_iso_density(dt1, dt2, pp, tsm, beta, nni, &ireg);

           cnn = (ireg) ? 1.0 / power_norm(tsm, beta, nni) : 0.0;

           if(wght){wt = *(wght + j); nn += wt;}
           else wt = 1.0;


           fd += (ftd = wt * cnn * ddd);

           if(!sw){
              if(d1) *d1 += ftd * dt2->sdt;
              else {*d2 += ftd * dt2->vec[0]; *d3 += ftd * dt2->vec[1];}
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
           else {
             if(nn > TOLWT) fd *= cons / nn;
             else fd = 1.0e-20;
           }

           llh += log(fd);

       }

   }

   if(!sw){
      if(d1) free(d1);
      if(d2){ free(d2); free(d3);}
   }

   return llh;

}
