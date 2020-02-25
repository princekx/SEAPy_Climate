#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "mem_er.h"
#include "constraint.h"

#define TOLPN   0.000001
#define TOLR    0.000001
#define TINY    1.0e-12

/* function to iteratively correct initial vector if constraints are violated */

void init_vec(double * , double * , double * , double * , int , int , int );
void swap(double * , double * , double * , int * , int , int , int , int );
void update_h(double * , double * ,int ,int , ...);
int update_nn1(double * , double * , int , int , int , int , ...);

extern struct cnst cst;

int correct_vect(double *mv, double *nn1, double *n, double *b, double *hh, int *acon, int *q, int dimv)

{

   int i, j, ic=0, iq;
   int ii, rpt, cf;

   double lambda[cst.numc], nl[*q], r[*q], *t1, *np;
   double *p, pn, mpn;
   double sum=0., sum1;

   rpt = 0;

   for(i=0; i < cst.numc; i++) lambda[i] = 0.;   /* initialization */

/* check for non-active constraints */

   for(i=cst.numeq; i < cst.numc; i++){

        cf = 1;

        for(j=cst.numeq; j < *q; j++) if(i+1 == acon[j]) {cf = 0; break;}

        if(cf){

           sum = 0.;

           t1 = cst.nn+i*dimv;

           for(j=0; j < dimv; j++) sum += *(t1+j) * mv[j];

           sum -= *(cst.bb+i);

           if(fabs(sum) < TINY) sum = TINY;


           if(sum < 0.) {rpt = 1; ic = i; break;}

         }

    }

    if(!(rpt)) return rpt;

    else{

      t1 = cst.nn + ic*dimv;

      for(i=0; i < *q; i++){

          nl[i] = 0.;

          ii = i*dimv;

          for(j=0; j < dimv; j++) nl[i] += *(t1+j) * n[ii+j];

      }

      for(i=0; i < *q; i++){

         r[i] = 0.;

         for(j=0; j < *q; j++) r[i] += nn1[j*cst.numc+i] * nl[j];

      }

      sum1 = 0.;
      for(i=0; i < dimv; i++){

          pn = 0.;
          for(j=0; j < *q; j++) pn += n[j*dimv + i] * r[j];

          mpn = *(t1+i) - pn;

          sum1 += mpn*mpn;
      }

      if(sum1 > TOLPN) {             /* add to constraint basis */

        lambda[*q] = sum;

        ii = *q * dimv;
        b[*q] = *(cst.bb+ic);
        acon[*q] = ic+1;
        for(j=0; j < dimv; j++) n[ii+j] = *(t1+j);
        np = n + ii;

        if(*q == 0)*nn1 = 1.; 
        else update_nn1(nn1, n, *q, dimv, cst.numc, 1);

        update_h(hh, np, dimv, 1);

        ++(*q);

        init_vec(mv, nn1, n, lambda, dimv, cst.numc, *q);

      }

      else {                        /* drop and add a constraint */

         iq = cst.numeq;

         for(i=cst.numeq; i < *q; i++) if(r[i] > TOLR) {iq = i+1; break;}

         if(iq > cst.numeq){

            if(iq < *q) swap(n, nn1, b, acon, iq, *q, dimv, cst.numc);

            p = (double * )calloc(dimv*dimv, sizeof(double));
            mem_er((p == NULL) ? 0 : 1, dimv*dimv * sizeof(double));

            update_nn1(nn1, n, *q, dimv, cst.numc, 0, p);

            update_h(hh, n+(*q-1)*dimv, dimv, 0, p);

            --(*q);

            acon[*q] = 0;

            free(p);

            lambda[*q] = sum;

            ii = *q * dimv;
            b[*q] = *(cst.bb+ic);
            acon[*q] = ic+1;

            for(j=0; j < dimv; j++) n[ii+j] = *(t1+j);
            np = n + ii;

            if(*q == 0)*nn1 = 1.;
            else update_nn1(nn1, n, *q, dimv, cst.numc, 1);


            update_h(hh, np, dimv, 1);

            ++(*q);

            init_vec(mv, nn1, n, lambda, dimv, cst.numc, *q);

         }

         else {              /* no feasable solution */

            printf("***error***, no feasable soultion for optimization by %s\n" 
                   "corrected vector returned\n", __FILE__);

            rpt = -1;

         }

      }

      return rpt;

    }

}

