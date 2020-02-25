#include <Stdio.h>
#include "statistic.h"
#include "p_vecs.h"
#include "sqt.h"

/* cross-validation part, using golden search method to optimize
   log-likelihood or least squares. Based on Diggle and Fisher:-

    P.J. Diggle and N.I. Fisher, "Sphere: A Contouring Program
    for Spherical Data", Computers and Geosciences, 11, 1985.    */


#define     RR        0.61803399
#define     CC        0.38196601
#define  SMTOL        0.001

double cross_validate(struct dpt *dt, int dtn, double sm, double *plt, float *wght, int im, int ht, double (*cvalln)(), int ind, struct cvecs *pcom, int nni, LEAF *lf, VEC *gv, int nlf, int sqt_c)

{

/*    static int ibr; */

    float tol;

    static double a, b;

    double fa=0., fb=0., fc, intv, xi1, xi2;
    double f1, f2;

    a = SMTOL;

    b = sm;

    if(ind){

      printf("specify a bracket for the golden search?\n\n");

      printf("input two values such that s1 < s2\n");
      scanf("%lf %lf", &a, &b);


    }

    if(pcom){

      fa = (*cvalln)(dt, dtn, a, plt, wght, im, ht, 0, 0, pcom, nni);
      fb = (*cvalln)(dt, dtn, b, plt, wght, im, ht, 0, 0, pcom, nni);

    }


    else if(lf){

      fa = (*cvalln)(dt, dtn, a, plt, wght, im, ht, 0, 0, lf, gv, nlf, 0, sqt_c, NULL);
      fb = (*cvalln)(dt, dtn, b, plt, wght, im, ht, 0, 0, lf, gv, nlf, 0, sqt_c, NULL);

    }

    else {

      fa = (*cvalln)(dt, dtn, a, plt, wght, im, ht, 0, 0);
      fb = (*cvalln)(dt, dtn, b, plt, wght, im, ht, 0, 0);

    }

printf("%f %f %f %f\n", a, fa, b, fb);

    while(fb > fa){

           b += b-a;
           fc = fb;

           if(pcom)

              fb = (*cvalln)(dt, dtn, b, plt, wght, im, ht, 0, 0, pcom, nni);

           else if(lf)

              fb = (*cvalln)(dt, dtn, b, plt, wght, im, ht, 0, 0, lf, gv, nlf, 0, sqt_c, NULL);

           else

              fb = (*cvalln)(dt, dtn, b, plt, wght, im, ht, 0, 0);


           if(fb < fc) break;

    }

    printf("***Information***, initial bracket for optimiztion \r\n"
           "                   is %f, %f\n\n", (float)a, (float)b);

    intv = b - a;
        
    xi1 = a + CC * intv;
    xi2 = a + RR * intv;

    if(pcom){

       f1 = (*cvalln)(dt, dtn, xi1, plt, wght, im, ht, 0, 0, pcom, nni);
       f2 = (*cvalln)(dt, dtn, xi2, plt, wght, im, ht, 0, 0, pcom, nni);

    }

    else if(lf){

       f1 = (*cvalln)(dt, dtn, xi1, plt, wght, im, ht, 0, 0, lf, gv, nlf, 0, sqt_c, NULL);
       f2 = (*cvalln)(dt, dtn, xi2, plt, wght, im, ht, 0, 0, lf, gv, nlf, 0, sqt_c, NULL);

    }

    else {

       f1 = (*cvalln)(dt, dtn, xi1, plt, wght, im, ht, 0, 0);
       f2 = (*cvalln)(dt, dtn, xi2, plt, wght, im, ht, 0, 0);

    }

/* printf("%f %f %f %f\n", a, fa, b, fb); */

    tol = (a + b) / 40.0;

    while(intv > tol){

          if(f1 > f2){

             b = xi2;
             fb = f2;
             xi2 = xi1;
             f2 = f1;
             intv = b - a;
             if(intv <= tol) break;
             xi1 = a + CC * intv;

             if(pcom)

                f1 = (*cvalln)(dt, dtn, xi1, plt, wght, im, ht, 0, 0, pcom, nni);

             else if(lf)

                f1 = (*cvalln)(dt, dtn, xi1, plt, wght, im, ht, 0, 0, lf, gv, nlf, 0, sqt_c, NULL);

             else

                f1 = (*cvalln)(dt, dtn, xi1, plt, wght, im, ht, 0, 0);

          }

          else{

             a = xi1;
             fa = f1;
             xi1 = xi2;
             f1 = f2;
             intv = b - a;
             if(intv <= tol) break;
             xi2 = a + RR * intv;

             if(pcom)

                f2 = (*cvalln)(dt, dtn, xi2, plt, wght, im, ht, 0, 0, pcom, nni);

             else if(lf)

                f2 = (*cvalln)(dt, dtn, xi2, plt, wght, im, ht, 0, 0, lf, gv, nlf, 0, sqt_c, NULL);

             else

                f2 = (*cvalln)(dt, dtn, xi2, plt, wght, im, ht, 0, 0);

printf("%f %f %f %f\n", a, fa, b, fb);
  
          }

          tol = (a + b) / 40.0;

    }

    if(lf) (*cvalln)(NULL, 0, 0, NULL, NULL, 0, 0, 0, 0, NULL, NULL, 0, 0, 0, NULL);

    printf("Final values of bracket are %f, %f\n\n", (float)a, (float)b);

    return 0.5 * (a + b);

}
