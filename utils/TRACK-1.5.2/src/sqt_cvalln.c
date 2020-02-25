#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "mem_er.h"
#include "m_values.h"
#include "sqt.h"
#include "tele.h"

/* calculate cross-validation function for a power kernal function based
   on SQT.                                                                */

extern int ipow;           /* exponent for power kernel */

double power_density(struct dpt * , struct dpt * , double , int );
int kernel_insqt(LEAF * , VEC * , struct dpt * , double );
LEAF **sqt_leaf_sample(LEAF ** , struct dpt * , LEAF * , VEC * , double , int * , int * , int );
double tele_kernel(double , double , double , int );

double sqt_cvalln(struct dpt *dt, int dtn, double sm, double *plt, float *wght, int im, int ht, int ty, int ms, LEAF *lf, VEC *gv, int nlf, int onet, int sqt_c, TELE *tele)

{

   int i, j, k, m;
   int did=0, dod=0;
   int insqt=0;
   int nlff=0;

   static int maxlf=0;

   int sw=0;

   float ipp;

   double cn, cons=0., cnn=0.0, cd;
   double ftd, dif;
   static double *fd=NULL;
   static double *d1=NULL, *d2=NULL, *d3=NULL;
   double llh=0.;
   double nn=0.0;
   double tsm=0., ts;
   double wt=1.0;

   double sang=0.0;

   struct dpt *dt1, *dt2;

   LEAF *lff=NULL, *lft=NULL;
   static LEAF **lnf=NULL;

   ipp = ipow + 1;

/* free memory if no data */

   if(!dt || !dtn){
      free(d1);
      free(d2);
      free(d3);

      free(fd);
      free(lnf);
      maxlf = 0;

   }


   if(ht < 0 || ht > 3){

      printf("***error***, non-valid switch type in %s \n\n", __FILE__);
      return 0.;

   }

   fd = (double * )calloc(dtn, sizeof(double));
   mem_er((fd == NULL) ? 0 : 1, dtn * sizeof(double)); 

   if(im > 0) sw = 1;

   if(!sw) {

      if(im == 0 && !d1) {

         d1 = (double * )calloc(dtn, sizeof(double));
         mem_er((d1 == NULL) ? 0 : 1, dtn * sizeof(double));

      }

      else if(im == -1){

         if(!d2){
            d2 = (double * )calloc(dtn, sizeof(double));
            mem_er((d2 == NULL) ? 0 : 1, dtn * sizeof(double));
         }

         if(!d3){
            d3 = (double * )calloc(dtn, sizeof(double));
            mem_er((d3 == NULL) ? 0 : 1, dtn * sizeof(double));
         }

      }


   }

   if(sw){
      for(i=0; i < dtn; i++) *(fd + i) = 1.0e-5;
   }
   else {
      for(i=0; i < dtn; i++) *(fd + i) = 1.0e-20;
      if(d1){
         for(i=0; i < dtn; i++) *(d1 + i) = 0;
      }
      else {
         for(i=0; i < dtn; i++) *(d2 + i) = *(d3 + i) = 0;
      }
   }


   if(!wght){
      nn = (ms) ? (float)dtn : (float) (dtn - 1);
      cn = ipp / (FPI2 * nn);

   }
   else {
      if(!tele) {
         for(i=0; i < dtn; i++) nn += *(wght + i);
         cn = (nn > TOLWT) ? ipp / FPI2 : 0.0;
      }
      else {
         nn = (ms) ? (float)dtn : (float) (dtn - 1);
         cn = ipp / (FPI2 * nn);
      }
   }


   if(ht < 3){

      if(sm < SMMIN) {tsm = -1.0; cons = 1.0; sang = -1.0;}

      else{

         tsm = (1.0 + sm) / sm;

         cn *= tsm;

         cd = tsm - 1.0;

         cons = cn / pow(cd, ipp);

         sang = sin(acos(1.0 / tsm));

         cnn = 1.0;

      }

   }

   else if(ht == 3) cons = cn;

   if(tele) cons /= tele->h;

   for(i=0; i < nlf; i++) {

       lff = lf + i;

       for(j=0; j < lff->ndata; j++){

           did = *(lff->ldata + j);

           if(wght) {
               wt = (!tele) ? *(wght + did) : tele_kernel(tele->tst, *(wght + did), tele->h, tele->ktyp);
           }
           else wt = 1.0;

           dt1 = dt + did;

           nlff = 0;

           if(ht == 3){

              ts = *(plt + did) * sm;

              if(ts < SMMIN) {tsm = -1.0; cnn = 1.0;}

              else {

                 tsm = (1.0 + ts) / ts;

                 cd = tsm - 1.0;

                 cnn = tsm / pow(cd, ipp);

                 sang = sin(acos(1.0 / tsm));

              }

           }


/* test for kernel support entirely within spherical triangle */

           if(sqt_c) {
              insqt = kernel_insqt(lff, gv, dt1, sang);

              if(insqt)
                lnf = sqt_leaf_sample(lnf, dt1, lff, gv, tsm, &nlff, &maxlf, 1);
              else 
                lnf = sqt_leaf_sample(lnf, dt1, lff, gv, tsm, &nlff, &maxlf, 0);

           }


/* test for neighbouring triangles */

           if(!sqt_c) 
             lnf = sqt_leaf_sample(lnf, dt1, lff, gv, tsm, &nlff, &maxlf, 0);


           for(k=0; k < nlff; k++){

              lft = lnf[k];
              lft->ingh = 0;

              for(m=0; m < lft->ndata; m++){

                  dod = *(lft->ldata + m);

                  if(!ms && did == dod) continue;

                  dt2 = dt + dod;

                  *(fd + dod) += (ftd = wt * cnn * power_density(dt1, dt2, tsm, ipow));

                  if(!sw){
                     if(d1) *(d1 + dod) += ftd * dt2->sdt;
                     else {
                       *(d2 + dod) += ftd * dt2->vec[0]; 
                       *(d3 + dod) += ftd * dt2->vec[1];
                     }

                  }

              }

           }


      } 

   }

   if(!sw){

     if(d1) {
        for(i=0; i < dtn; i++){
           ftd = *(fd + i);
           *(d1 + i) = (ftd > DTOL) ? *(d1 + i) / ftd : 0.0;
           dif = (dt + i)->sdt - *(d1 + i);
           if(ms) (dt + i)->dsdt = *(d1 + i);
           llh -= (dif*dif);
       }
     }
     else{ 
        for(i=0; i < dtn; i++){
           ftd = *(fd + i);
           if(ftd > DTOL) {*(d2 + i) /= ftd; *(d3 + i) /= ftd;}
           else { *(d2 + i) = 0.0; *(d3 + i) = 0.0;} 
           dif = (dt + i)->vec[0] - *(d2 + i);
           llh -= (dif*dif);
           dif = (dt + i)->vec[1] - *(d3 + i);
           llh -= (dif*dif);
           if(ms) {(dt + i)->dvec[0] = *(d2 + i); (dt + i)->dvec[1] = *(d3 + i);}
        }
     }


   }

   else { 

     for(i=0; i < dtn; i++){
         if(!wght) *(fd + i) *= cons;
         else {
            if(!tele) {
              *(fd + i) *= cons / ((ms) ? nn : (nn - *(wght + i)));
            }
            else {
              *(fd + i) *= cons;
            }
         }

         ftd = *(fd + i);

         llh += log(ftd);

         if(ty) *(plt + i) = ftd;

     }

  }

   if(onet){
      free(d1);
      free(d2);
      free(d3);
      free(fd);
      free(lnf);
      maxlf = 0;

   }

   return llh;

}

