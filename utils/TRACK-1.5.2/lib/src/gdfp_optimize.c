#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "constraint.h"

/* this routine uses the Goldfarb-Davidon-Fletcher-Powell variable
   metric method to find the local minimum (maximum) of a function
   under linear equality and inequality constraints.        
   See:- "Extension of Davidon's Variable Metric Method to Maximization
          Under Linear Inequality and Equality Constraints", D. Goldfarb,
          SIAM J. Appl. Math. 17, pp 739, 1969.                           

   extensions to this algorithm have been incorporated based on the 
   techniques devised by Rosen to find a initial feasable point, and
   for dealing with linear dependance of inequality constraints.

   See:- "The Gradient Projection Method for Nonlinear Programming.
          Part I. Linear Constraints", L. B. Rosen, J. Soc. Indust. Appl.
          Math., 8, pp 181, 1960.                                        */

#define   TOLDIF    0.00001     /* tolerance on gradient  */
#define   TOLINT    0.00001      /* tolerance on the line interval */
#define   TOLC      0.00000001  /* constraint tolerance */
#define   TOLF      0.000001    /* tolerance on function change */
#define   TOLGAM    0.0001      /* tolerance on gamma upper bound */
#define   KMAG1     0.000001    /* gamma upper bounds KMAG1 < KMAG2 */
#define   KMAG2     0.5
#define   SMALL     1.0e-20
#define   MAXLAM    1.0e+8
#define   MAXITER   20


int correct_vect(double * , double * , double * , double * , double * , int * , int * , int );
double func(double * , double * , int , int , int );
void init_vec(double * , double * , double * , double * , int , int , int );
void swap(double * , double * , double * , int * , int , int , int , int );
void update_h(double * , double * ,int ,int , ...);
int update_nn1(double * , double * , int , int , int , int , ...);

extern struct cnst cst;
extern double *ih, *inn;

double gdfp_optimize(double *mvec, double fm, int dimv, int type, int iopt, int *ier)

{

    int ii, jj, i, j, q=0, iq;
    int kk=dimv*dimv;
    int ncc=cst.numc*cst.numc;
    int dc=0, da=0;
    int acon[cst.numc];
    int cf, ic=0;
    int iter, cg, ld;
    int icor=0, rpt;

    double val0, hm;
    double dvec[dimv], h[dimv], hh[kk], p[kk], mv[dimv];
    double mvc[dimv], dv[dimv], y[dimv], sig[dimv];
    double alpha[cst.numc], al[cst.numc];
    double nn1[ncc], *t1;
    double lam, gam, gam0, gam1;
    double f0, f1, f, fmax;
    double w, d, z, dd;

    double n[dimv*cst.numc], b[cst.numc], *np;
    double sum, sum1, aq=1., aa, sm;

/* initialize error flag */

   *ier = 0;
   val0 = 0.0;

/* initialize constraint basis by checking which constraints are active.
   equality constraints should have been pre-checked to ensure the initial
   feasable point satisfies them.                                         */

   for(i=0; i < cst.numc; i++) acon[i] = 0;

   for(i=0; i < dimv; i++) mv[i] = mvec[i];

/* put equality constraints into constraint basis */

   if(cst.numeq > 0){

      for(i=0; i < cst.numeq; i++){

          ii = i*dimv;

          b[i] = *(cst.bb+i);

          t1 = cst.nn + ii;

          acon[i] = i+1;

          for(j=0; j < dimv; j++) n[ii+j] = *(t1+j);

/* check equality constraints are satisfied by the starting vector */

          sum = 0.;

          for(j=0; j < dimv; j++) sum += n[ii+j] * mvec[j];


          alpha[i] = sum -= *(cst.bb+i);

          if(fabs(sum) > TOLC && !(icor)){

             printf("***error*** starting vector does not satisfy equality\r\n"
                    "            constraints for optimization by routine  \r\n"
                    "            gdfp_optimize, correction will be made.\n");

             icor = 1;

          }

      }

    }

/* initilise the Hq0 , nn1 matricies */

    if(iopt == 0 ){    /* if iopt = 1 use a previously initialised 
                             strict equality constraint basis.             */

       for(i=0; i < dimv; i++){

           ii = i*dimv;

           for(j=0; j < dimv; j++) {jj = ii++; hh[jj] = (i == j) ? 1. : 0.;}

       }

       for(i=0; i < cst.numeq; i++) {

          np = n+i*dimv;

          if(q == 0){nn1[0] = 1.; ld = 0;}
          else ld = update_nn1(nn1, n, q, dimv, cst.numc, 1);

          if(ld){

            printf("***error***, linearly dependant set of equality constraints\r\n"
                   "             choose a smaller linearly independant set.       \n");
            exit(1);

          }

          update_h(hh, np, dimv, 1);

          ++q;

       }

       if(cst.numeq){

          for(i=0; i < kk; i++) *(ih + i) = hh[i];

          for(i=0; i < ncc; i++) *(inn + i) = *(nn1+i);

       }


    }

    else {

       for(i=0; i < kk; i++) hh[i] = *(ih + i);
       for(i=0; i < ncc; i++) *(nn1+i) = *(inn + i);
    }

/* correct initial feasable point if required */

    if(icor){

       init_vec(mvec, nn1, n, alpha, dimv, cst.numc, cst.numeq);

       for(i=0; i < cst.numeq; i++) alpha[i] = 0.;

    }

    icor = 0;

/* check which inequality constraints satisfy equality and put in 
   constraint basis.                                                */

    for(i=cst.numeq; i < cst.numc; i++){

        sum = 0.;

        t1 = cst.nn+i*dimv;

        for(j=0; j < dimv; j++) sum += *(t1+j) * mvec[j];

        sum -= *(cst.bb+i);


        if(sum < TOLC ){

          if(sum < 0. ) icor = 1;

          ii = q*dimv;
          jj = i *dimv;
          b[q] = *(cst.bb+i);
          t1 = cst.nn + jj;
          acon[q] = i+1;
          for(j=0; j < dimv; j++) n[ii+j] = *(t1+j);
          np = n + ii;

          if(q == 0){nn1[0] = 1.; ld = 0;}
          else ld = update_nn1(nn1, n, q, dimv, cst.numc, 1);

          if(ld) acon[q] = 0;

          else{

            update_h(hh, np, dimv, 1);
            alpha[q] = sum;
            ++q;

          }

        }

    }

/* correct feasable point for ineqaulity constarints if required */

    if(icor){

       init_vec(mvec, nn1, n, alpha, dimv, cst.numc, q);

       rpt = 1;

       while(rpt){

          rpt = correct_vect(mvec, nn1, n, b, hh, acon, &q, dimv);

          if(rpt < 0) return func(mvec, dvec, dimv, type, 0);

       }
      

    }


/* compute gradient at the initial point */

    val0 = func(mvec, dvec, dimv, type, 1);

    iter = 0;

loop_start:

    for(i=0; i<dimv; i++) dv[i] = dvec[i];

    dc = 0;

/* compute h = Hg */

    for(i=0; i < dimv; i++){

       h[i] = 0.;

       for(j=0; j < dimv; j++) h[i] += hh[j*dimv + i] * dvec[j];

       if(fabs(h[i]) < TOLDIF) ++dc;

     }

/* compute alpha */

     da = 0;

     for(i=0; i < q ; i++){

       al[i] = 0.;
       ii = i * dimv;

       for(j=0; j < dimv; j++) al[i] += n[ii + j] * dvec[j]; 

     }

     for(i=cst.numeq; i < q; i++){

       alpha[i] = 0.;

       for(j=0; j < q; j++) alpha[i] += al[j] * nn1[j*cst.numc + i];

       if(alpha[i] <= 0.) ++da;

     }

/* check for constrained solution */

     if(dc == dimv && da == q-cst.numeq) return val0;

/* end of step 1 */

/* begin step 2 */

     sum = 0.;

     for(i=0; i < dimv; i++) sum += h[i]*h[i];

     hm = sum;

     sum = sqrt(fabs(sum));

     aq = 0.;
     iq = 0;

     for(i=cst.numeq; i < q; i++){

       if((aa=0.5*alpha[i]/ sqrt(fabs(nn1[i*cst.numc+i]))) >= aq){aq = aa; iq = i+1;}

     }

     if(sum <= aq){

/* drop a constraint from the basis and update operators */

       if(iq < q) swap(n, nn1, b, acon, iq, q, dimv, cst.numc); /* swap constraints */

       update_nn1(nn1, n, q, dimv, cst.numc, 0, p);

       update_h(hh, n+(q-1)*dimv, dimv, 0, p);

       --q;

       acon[q] = 0;
       goto loop_start;

    }

/* begin step 3 */

    lam = MAXLAM;

    for(i=cst.numeq; i < cst.numc; i++){

       cf = 1;

       for(j=cst.numeq; j < q; j++)

           if(i+1 == acon[j]) {cf = 0; break;}

       if(cf){

          sum = 0.;
          t1 = cst.nn + i*dimv;
          for(j=0; j < dimv; j++) sum -= *(t1+j) * mvec[j];

          sum += *(cst.bb + i);

          sum1 = 0.;

          for(j=0; j < dimv; j++) sum1 += *(t1+j) * h[j];

          sum /= sum1;

          if(sum > 0. && sum < lam) {lam = sum; ic = i;}

       }

    }

    if(lam > MAXLAM) {

      printf("***error***, in %s, optimization may be unbounded in direction of search\n", __FILE__);
      *ier = 1;
      return 0.0;


    }

    for(i=0; i < dimv; i++) mvc[i] = mvec[i];

/* calculate the maximum along a line */

    gam0 = 0.;

    f0 = val0;

/* compute (g)Hg and initial point */

    sum = 0.;

    for(i=0; i<dimv; i++) sum += h[i] * dv[i];

    if(sum > 0.) {

       if(f0 > fm)gam1 = 2.0 * (f0 - fm) / sum;
       else gam1 = KMAG2;

       if(fabs(f0-fm) < TOLGAM) gam1 = KMAG1;

       if(gam1 > KMAG2) gam1 = KMAG2;

       if(gam1 > lam) gam1 = lam;

    }

    else{

        printf("***Warning*** in %s line %d, no max. on line within limits\n",__FILE__, __LINE__);
        *ier = 2;
        return val0;

    }

    sum1 = 1.;

    while(sum1 > 0.){

        for(i=0; i<dimv; i++) mv[i] = mvc[i] + gam1 *h[i];

        f1 = func(mv, dvec, dimv, type, 1);

        sum1 = 0.;
 
        for(i=0; i < dimv; i++) sum1 += dvec[i] * h[i];

        if((lam - gam1) < TOLINT)break;

        if(sum1 < 0. || f1 < f0) break;
        else gam1 *= 2.;

        if(gam1 > lam) {gam1 = lam;}

    }

    if(sum1 > 0. && f1 > f0)

       for(i=0; i < dimv; i++) {sig[i] = gam1 * h[i]; mvec[i] += sig[i];}

    else{

      gam = gam1;

      while(gam > TOLINT){

        z = ((f0 - f1) * 3.0 / ((gam1 > SMALL) ? gam1 : SMALL)) + sum + sum1;
        dd = z*z - sum * sum1;
        if(dd < SMALL) dd = SMALL;

        w = sqrt(dd);

        dd = sum - sum1 + 2.0 * w;
        dd = (dd < SMALL) ? SMALL : dd;
        d = (w - sum1 + z)/dd;

        if(d < 0.) d = 0.;
        gam = gam1 * (1.0 - d);
 
        for(i=0; i<dimv; i++) mv[i] = mvc[i] + gam *h[i];

        if(gam < TOLINT) break;

        f = func(mv, dvec, dimv, type, 1);
 
        sm = 0.;

        for(i=0; i < dimv; i++) sm += dvec[i] * h[i];

        fmax = (f0 < f1) ? f1 : f0;

        cg = 1;

        if(f <= fmax){

           if(f0 < f1) cg = 1;
           else cg = 0;

        }

        if((sm > 0. && cg == 1) || (f <= fmax && cg == 1)){

          f0 = f;
          sum = sm;
          gam1 -= gam;
          for(i=0; i<dimv; i++) mvc[i] = mv[i];

        }

        else {

          f1 = f;
          sum1 = sm;
          gam1 = gam;


        }

      }

      for(i=0; i < dimv; i++){sig[i] = mv[i] - mvec[i]; mvec[i] = mv[i];}

      sum = 0.;
      for(i=0; i<dimv; i++) sum += sig[i] * h[i];

      gam1 = sum / hm;

    }

    fm = f = val0;

    val0 = func(mvec, dvec, dimv, type, 1);

    if(f-val0 > TOLF){

       printf("***Warning***, in file %s, line %d, an illegal step was taken during local optimization\n", __FILE__, __LINE__);
       *ier = 3;

       return val0;

    }
     
    sum = 0.;

    for(i=0; i < dimv; i++) sum += dvec[i] * h[i];

/* update H */

    if((lam - gam1) < TOLINT && sum >= 0.) {   /* add to constraint basis */

          ii = q*dimv;
          jj = ic *dimv;
          b[q] = *(cst.bb+ic);
          t1 = cst.nn + jj;
          acon[q] = ic+1;
          for(j=0; j < dimv; j++) n[ii+j] = *(t1+j);

          np = n + ii;

          if(q == 0){ nn1[0] = 1.; ld = 0;}
          else ld = update_nn1(nn1, n, q, dimv, cst.numc, 1);

          if(ld)acon[q] = 0;

          else{

            update_h(hh, np, dimv, 1);

            ++q;

          }

/* check and correct active constraints, may be needed if to much 
   roundof error builds up in which case un-comment between dotted 
   lines.                                                             */

/* ....................................................................         
          icor = 0;
          for(i=0; i < q; i++){

             ii = i * dimv;
             alpha[i] = 0.;
             for(j=0; j < dimv; j++) alpha[i] += n[ii+j] * mvec[j];
             alpha[i] -= b[i];
             if(alpha[i] < 0.) icor = 1;

          }

          if(icor){
  

            init_vec(mvec, nn1, n, alpha, dimv, cst.numc, q);

            rpt = 1;

            while(rpt){

               rpt = correct_vect(mvec, nn1, n, b, hh, acon, &q, dimv);

               if(rpt < 0) {

                  for(i=0; i < dimv; i++) mvec[i] = mv[i];

                  return func(mvec, dvec, dimv, type, 0);

               }

            }

          } 

........................................................................*/

    }

    else { 


          for(i=0; i < dimv; i++) y[i] = dvec[i] - dv[i];
 
          update_h(hh, y, dimv, 1);

          sum = 0.;
          for(i=0; i < dimv; i++) sum += sig[i] * y[i];

          if(fabs(sum) < SMALL) sum = (sum < 0.) ? -SMALL : SMALL;

          for(i=0; i < dimv; i++){

              for(j=0; j < dimv; j++) *(hh + j*dimv + i) -= sig[i] * sig[j] / sum;

          }


    }

    ++iter;

    if(iter > MAXITER) {

      printf("***WARNING***, maximum number of iterations in %s exceeded\n", __FILE__);
      *ier = 4;

      return val0;

    }

    goto loop_start;                
 
    return val0;

}
