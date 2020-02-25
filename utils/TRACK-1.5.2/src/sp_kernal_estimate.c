#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "mem_er.h"
#include "statistic.h"
#include "m_values.h"
#include "p_vecs.h"
#include "sqt.h"
#include "tele.h"


#define  BETA  0.45

/* function to perform various kernal estimates, densities and means.
   

   For the cross-validation algorithm see:-

        N. I. Fisher, T. Lewis and B. J. J. Embleton, 
        "Statistical Analysis of Spherical Data", Cambridge University Press,
        1987.

   or:-

        P. J. Diggle and N. I. Fisher, "Sphere: A contouring Program for
        Spherical Data", Computers and Geosciences, 11, pp 725-766, 1985.

   For the basis of the adaptive scheme see:-

        B. W. Silverman, Density Estimation for Statistics and Data Analysis",
        Chapman and Hall, 1990.

                                                                             
   Key for variable 'im', controls type of kernal estimation:-

        im > 0 .... single density estimation of smoothing paramater, i.e.,
                    smoothing parameter determined and used for regression
                    as well.

           im = 1,  density estimation.
           im = 2,  scaler regression estimation.
           im = 3,  vector (2D) regression estimation.

        im <= 0 ... individual regression estimation of smoothing
                    parameter.

           im = 0,  scaler regression estimation.
           im = -1,  vector (2D) regression estimation.

                                                                              */

extern int tom;
extern int trd;
int ipow=1;           /* exponent for power kernel */


double non_lin_sol(double );
double cvalln_fisher(struct dpt * , int , double , double * , float * , int , int , int , int );
double eval_fisher_stat(double *, double ** , struct dpt * , struct dpt * , float * , int , int , int , int , int , int );
double cvalln_linear(struct dpt * , int , double , double * , float * ,int , int , int , int );
double eval_linear_stat(double *, double ** , struct dpt * , struct dpt * , float * , int , int , int , int , int , int );
double cvalln_quad(struct dpt * , int , double , double * , float * ,int , int , int , int );
double eval_quad_stat(double *, double ** , struct dpt * , struct dpt * , float * , int , int , int , int , int , int );
double cross_validate(struct dpt * , int , double , double * , float * , int , int , double (*cvalln)(), int , struct cvecs * , int , LEAF * , VEC * , int , int );
double cvalln_const(struct dpt * , int , double , double * , float * ,int , int , int , int );
double eval_const_stat(double *, double ** , struct dpt * , struct dpt * , float * , int , int , int , int , int , int );
double cvalln_non_iso(struct dpt * , int , double , double * , float * ,int , int , int , int , struct cvecs * ,int );
double eval_non_iso_stat(double *, double ** , struct dpt * , struct dpt * , float * , int , int , int , int , int , int , struct cvecs * ,int );

double sqt_cvalln(struct dpt * , int , double , double * , float * , int , int , int , int , LEAF * , VEC * , int , int, int , TELE * );
double sqt_eval_stat(double *, double ** , struct dpt * , struct dpt * , float * , int , int , int , int , int , int , LEAF * , VEC * , int , int , TELE * );

void residual(double (*cvalln)(), struct dpt * , int , double , double * , float * , int , int , struct cvecs * , int , LEAF * , VEC * , int , int );

double non_iso_adapt(struct dpt * , int , double , double * , float * , struct cvecs * , int );

void regfal(double * , double * , double);


float sp_kernal_estimate(double **den, struct dpt *st, struct dpt *dt, float *wght, int dtn, int ptnum, int im, int *kty, int ims, float *nn, LEAF *lf, VEC *gv, int nlf, int ninfo, TELE *tele)

{

    int i;
    int ms=0;
    int swrn=0;
    int nni=0;
    int sqt_c=0;

    static int ht;
    static int ks=0;
    static int icv='n';

    static int ind=1;

    float n1, n2;
    float hbw=0.;

    static float al;

    static double sm, dsm;

    float plmin, plmax;

    double dsm1, dsm2, tsm=0.;
    double r;
    double sx1, sy1, sz1, sx2, sy2, sz2;
    double *plt=NULL, *pl;
    double gg;
    double newk, newb;
    float *cv_wt=NULL;

    struct cvecs *pcom=NULL, *pp=NULL;


    static double (*eval_sample_stat)()=NULL;
    static double (*cvalln)()=NULL;
 
    struct dpt *dtt;

    if(!trd && ind == 0) ind = 1;

    if(tom != 'g'){

       printf("***WARNING***, incorrect distance measure specified\r\n"
              "               changing to geodesic measure\n\n");
       tom = 'g';

    }

    if(dtn <= 1 || !dt) {

       printf("***WARNING***, no data for statistical analysis in file.\n\n");
       return 0.0;

    }


    if(ind){

       if(wght){
         if(tele) {
            printf("****INFORMATION****, additional raw weight data is available, this can only       \r\n"
                   "                     currently be used by the SQT estimator, but cannot currently \r\n"
                   "                     be used in cross-validation or robust estimation.            \n\n");
         }
         else {
            printf("Additional weight data is available, do you want to use this\r\n"
                   "in cross-validation or robust estimation of means, 'y' or 'n'\n\n");
            scanf("\n");
            if(getchar() == 'y') cv_wt = wght;
            else cv_wt = NULL;
         }
       }

       if(!lf) {

          printf("which smoothing kernal do you require,                \r\n"
                 "        input ' 0', for exponential (fisher),         \r\n"
                 "        input ' 1', for constant (spherical cap)      \r\n"
                 "        input ' 2', for linear,                       \r\n"
                 "        input ' 3', for quadratic.                    \r\n"
                 "        input '10', for adaptive non-isotropic power. \n\n");

          scanf("%d", &ks);

          switch(ks){
             case 0:
               eval_sample_stat = eval_fisher_stat;
               cvalln = cvalln_fisher;
               *kty = 1;
               break;
             case 1:
               eval_sample_stat = eval_const_stat;
               cvalln = cvalln_const;
               *kty = 2;
               break;
             case 2:
               eval_sample_stat = eval_linear_stat;
               cvalln = cvalln_linear;
               *kty = 3;
               break;
             case 3:
               eval_sample_stat = eval_quad_stat;
               cvalln = cvalln_quad;
               *kty = 4;
               break;
             case 10:
               eval_sample_stat = eval_non_iso_stat;
               cvalln = cvalln_non_iso;
               *kty = 5;
               ht = 3;
               printf("what index is required for the non-isotropic power kernel, 0-4. Integer only for the moment!\n");
               scanf("%d", &nni);
               if(nni < 0 || nni > 4){
                  printf("****WARNING****, interger value out of range, defaulting to 1\n\n");
                  nni = 1;
               }
               break;
             default:
               printf("***error***, illegal index for function pointer assignment in %s\n\n", __FILE__);
               exit(1);
          }

       }

       else {
          eval_sample_stat = sqt_eval_stat;
          cvalln = sqt_cvalln;
          ks = 20;
          *kty = 6;
          printf("The SQT estimator uses the power kernel, what exponent is required (integer).\n\n");
          scanf("%d", &ipow); 

          printf("Do you want to test for kernels entirely encompased by a spherical   \r\n"
                 "triangle, 'y' or 'n', this may improve the estimation for a coarse   \r\n"
                 "decomposition but will probably slow the estimation down for a fine  \r\n"
                 "decomposition.                                                       \n\n");
          scanf("\n");
          if(getchar() == 'y') sqt_c = 1;

       }


       if(ks != 10){

          printf("do you want the default global smoothing parameter '0', \r\n"
                 "            user chosen global smoothing parameter '1', \r\n"
                 "                a value chosen by cross-validation '2', \r\n"
                 "                             or adaptive smoothing '3'. \n\n");

          scanf("%d", &ht);

       }

    }


    if(ht){

       if(ht == 3){

          if(ind && ninfo) printf("***INFORMATION***, see Silverman for advise on parameter settings\n\n");

       }


       if(ind){

          printf("input an initial value for the global smoothing parameter\r\n"
                 "also used for upper bracket value for cross-validation\n\n");

          scanf("%lf", &sm);


       }

    }

/* compute default value for global smoothing parameter based on fisher density */


    if(!ks){

       sx1 = sy1 = sz1 = 0.;
       sx2 = sy2 = sz2 = 0.;
       n1 = n2 = 0.;

       for(i=0; i < dtn; i++){

          dtt = dt + i;

          if(dtt->ydt >= 0.){
             sx1 += dtt->xdt;
             sy1 += dtt->ydt;
             sz1 += dtt->zdt;
             n1 += 1.;
          }
          else {
             sx2 += dtt->xdt;
             sy2 += dtt->ydt;
             sz2 += dtt->zdt;
             n2 += 1.;
          }


       }

       r = sqrt(sx1 * sx1 + sy1 * sy1 + sz1 * sz1) / n1;

       dsm1 = non_lin_sol(r) * pow(n1, INV3);

       r = sqrt(sx2 * sx2 + sy2 * sy2 + sz2 * sz2) / n2;

       dsm2 = non_lin_sol(r) * pow(n2, INV3);

       dsm = (dsm1 < dsm2) ? dsm2 : dsm1;

       printf("Default global smoothing factor based on fisher density = %f\n\n", (float)dsm);

    }

    else if(ind){


       printf("input a default smoothing value\n\n");
       scanf("%lf", &dsm);

    }

    if(!ims){

       if(im > 1 || im <= 0){

          printf("Do you want a robust estimate, 'y' or 'n'\n");
          scanf("\n");
          if(getchar() == 'y') ms = 1;
          else ms = 0;

       }

    }

    if(sm < dsm){

      printf("***WARNING***, value of smoothing parameter       \r\n"
             "               is smaller than default value,     \r\n"
             "               default value is used.             \n\n");

      sm = dsm;

    }

    if(ht == 0 || ht == 1){

       if(!ht) printf("***WARNING***, using the default value may lead to gross oversmoothing\n\n");

       tsm = (ht) ? sm : dsm;


       if(ms) residual(cvalln, dt, dtn, tsm, NULL, cv_wt, ht, im, NULL, nni, lf, gv, nlf, sqt_c);

       if(lf)

         *nn = (*eval_sample_stat)(&tsm, den, st, dt, wght, dtn, ptnum, im, ht, ind, ms, lf, gv, nlf, sqt_c, tele);

       else
 
         *nn = (*eval_sample_stat)(&tsm, den, st, dt, wght, dtn, ptnum, im, ht, ind, ms);


    }

    else if(ht == 2){

        tsm = cross_validate(dt, dtn, sm, NULL, cv_wt, im, ht, cvalln, ind, NULL, nni, lf, gv, nlf, sqt_c);

        if(tsm < dsm) {

          printf("***WARNING***, value of smoothing parameter from \r\n"
                 "               cross-validation is smaller than default\r\n"
                 "               value, default value is used.\n\n");

          tsm = dsm;

        }

        if(ms) residual(cvalln, dt, dtn, tsm, NULL, cv_wt, ht, im, NULL, nni, lf, gv, nlf, sqt_c);

        if(lf)

           *nn = (*eval_sample_stat)(&tsm, den, st, dt, wght, dtn, ptnum, im, ht, ind, ms, lf, gv, nlf, sqt_c, tele);

        else 

           *nn = (*eval_sample_stat)(&tsm, den, st, dt, wght, dtn, ptnum, im, ht, ind, ms);
 
    }

    else if(ht == 3){

/* compute a pilot estimate using a user chosen smoothing parameter */

        tsm = sm;

/* assign memory for pilot estimate */

        plt = (double * )calloc(dtn, sizeof(double));
        mem_er((plt == NULL) ? 0 : 1, dtn * sizeof(double));


/* assign storage for eigen-vectors and eigenvalues. */

        if(ks == 10){

           pcom = (struct cvecs * )calloc(dtn, sizeof(struct cvecs));
           if((pcom == NULL) ? 0 : 1);
           gg = non_iso_adapt(dt, dtn, tsm, plt, cv_wt, pcom, ind);


/* for non-isotropic kernel estimation assign relative values of beta */


           for(i=0; i< dtn; i++){

               pp = pcom + i;

               if(pp->p3[3] < EIGTOL) pp->beta = 0.0;
               else pp->beta = 0.5 * (1.0 - (pp->p2[3] / pp->p3[3]));

               if(pp->beta > BETA) pp->beta=BETA;

           }

        }

        else {

            if(lf)

               gg = (*cvalln)(dt, dtn, tsm, plt, cv_wt, 1, 1, 1, 0, lf, gv, nlf, 1, sqt_c, tele);

            else

               gg = (*cvalln)(dt, dtn, tsm, plt, cv_wt, 1, 1, 1, 0);

        }

        gg = exp(gg / (float)dtn);

        if(ind){

           printf("what sensitivity parameter is required for pilot estimate 0<=alpha<=1\n\n");

           scanf("%f", &al);

        }

/* calculate local band width factors */

        for(i=0; i < dtn; i++){

           pl = plt+i;
           *pl /= gg;
           *pl = pow(*pl, (double)al);

        }

/* find min and max band width factors */

        plmin = plmax = *plt;

        for(i=0; i < dtn; i++){

           pl = plt+i;
           if(*pl < plmin) plmin = *pl;
           else if(*pl > plmax) plmax = *pl;


        }

        if(ninfo) {
           printf("****INFORMATION****, the range of the adaptive band width\r\n"
                  "                     factors is %f to %f\n\n", plmin, plmax);
        }


        if(ind){

          printf("do you want cross validation, 'y' or 'n'\r\n"
                 "***WARNING***, this can be computationally expensive\n\n");

          scanf("\n");
          icv=getchar();

        }

        if(icv == 'y')

           tsm = cross_validate(dt, dtn, sm, plt, cv_wt, im, ht, cvalln, ind, pcom, nni, lf, gv, nlf, sqt_c);


        if(ms) residual(cvalln, dt, dtn, tsm, plt, cv_wt, ht, im, pcom, nni, lf, gv, nlf, sqt_c);


/* scale smoothing factors with global smoothing parameter and check
   against lower bound value.                                           */


        for(i=0; i < dtn; i++) {

           pl = plt + i;

           *pl *= tsm;

           if(*pl < dsm){swrn = 1; *pl = dsm;}

        } 

        if(ind && swrn){

           printf("***WARNING***, some local bandwidth are less         \r\n"
                  "               than the default value and have been  \r\n"
                  "               adjusted to the default value.        \n\n");


        }


        if(ks == 10){

/* coorect bandwidth */

           for(i=0; i< dtn; i++){

               pl = plt + i;
               pp = pcom + i;

               newk = (1.0 + *pl) / *pl;
               newb = pp->beta * newk;


               regfal(&newk, &newb, BETA);

               pp->beta = newb / newk;
               *pl = 1.0 / (newk - 1.0);
           }

           *nn = (*eval_sample_stat)(plt, den, st, dt, wght, dtn, ptnum, im, ht, ind, ms, pcom, nni); 

           free(pcom);

        }

        else{      

           if(lf)

              *nn = (*eval_sample_stat)(plt, den, st, dt, wght, dtn, ptnum, im, ht, ind, ms, lf, gv, nlf, sqt_c, tele);

           else

              *nn = (*eval_sample_stat)(plt, den, st, dt, wght, dtn, ptnum, im, ht, ind, ms);


        }

        free(plt);

    }

    switch(ks){
       case 0:
          hbw = acos(1.0 - log(pow(2.0, 1.0/tsm))) / FP_PI;
          break;
       case 1: case 2: case 3: case 10: case 20:
          hbw = acos(tsm / (1.0 + tsm)) / FP_PI;
          break;
       default:
          printf("***error***, invalid switch in %s\n", __FILE__);
          break;

    }

    if(ind){

       printf("***Information***, the half band width (HBW) for a global \r\n"
              "                   smoothing parameter %f is equal to %f  \r\n"
              "                   degrees\n\n", (float)tsm, hbw);

    }

    if(trd && ind == 1) ind = 0;

    return (float)hbw;

}

/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

              SUPPORT ROUTINES

   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

#define  MAXIT   10
#define  DXTOL   0.0001 

/* newton-raphson solution for root of fisher density maximum likelihood
   function.

           coth(x) - 1/x = R/n

                                                                           */

double non_lin_sol(double r)

{
    int it=0;

    double rt=r, sch, csh, dt;

    while(it++ < MAXIT){


       sch = sinh(rt);
       csh = cosh(rt);

/* if no solution returns 1 by default */

       if(rt <= 0.) return rt = 1.0;

       dt = ((csh/sch) - (1.0/rt) - r) / ((1.0/(sch*sch)) - (1.0 / (rt * rt)));

       rt += dt;

       if(fabs(dt) < DXTOL) break;

    }


    return rt;

}



