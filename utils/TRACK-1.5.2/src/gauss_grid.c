#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "grid.h"
#include "m_values.h"
#include "mem_er.h"

#define  LTOL 1.0e-10
#define  NIT  20

double legendre_node(int , int , double * );


/* function to define a Gaussian grid for fitting data from a non-Gaussian grid */

GRID *gauss_grid(int ntr)

{

    int nlat, nlng;
    int nlat2;
    int i;

    GRID *newg;

    float dx;

    double dum=0.0;

    printf("****INFORMATION****, quadratic gaussian grid only.\n\n");

    nlng = 3 * ntr + 1;
    if(nlng % 2) nlng += 1;
    if(nlng % 4) nlng += 2;

    nlat = nlng / 2;
    nlat2 = nlat / 2;

    printf("NLNG=%d; NLAT=%d\n\n", nlng, nlat);

    ++nlng;

/* allocate memory */

    newg = (GRID *)malloc_initl(sizeof(GRID));
    mem_er((newg == NULL) ? 0 : 1, sizeof(GRID));

    newg->ix = nlng;
    newg->iy = nlat;

    newg->xgrid = (float *)calloc(nlng, sizeof(float));
    mem_er((newg->xgrid == NULL) ? 0 : 1, nlng * sizeof(float));

    newg->ygrid = (float *)calloc(nlat, sizeof(float));
    mem_er((newg->ygrid == NULL) ? 0 : 1, nlat * sizeof(float));


    dx = 360.0 / (nlng - 1);

    for(i=0; i < nlng; i++){

       *(newg->xgrid + i) = (float) i * dx;

    }


    for(i=0; i < nlat2; i++){

        *(newg->ygrid + nlat - i - 1) = (FP_PI2 - acos(legendre_node(i+1, nlat, &dum))) / FP_PI;
        *(newg->ygrid + i) =  - *(newg->ygrid + nlat - i - 1);

    }


    return newg;


}


/* function to compute the Legendre nodes */

double legendre_node(int ilat, int nlat, double *wght)

{

    int i, j, nit;
    

    double sa, sb, sc;
    double bn, th, dth, xx, gg;
    double amm, amn, ann, em;
    double re, aa=0., dd, dlt;
    
    sa = sqrt(0.5);
    sb = sqrt(1.5);
    sc = sqrt(1.0 / 3.0);

    bn = (double) nlat;

    gg = (2.0 * bn + 1.0);

    th = FPI * (2.0 * (double) ilat - 0.5) / gg;
    dth = th + (cos(th) / (8.0 * bn * bn * sin(th)));
    xx = cos(dth);

    nit = 0;

    for(i=0; i<NIT; i++){

        ++nit;
/* iteration for node */

        amm = sa;
        amn = xx * sb;
        em = sc;
        ann = 1.0;

        for(j=1; j < nlat; j++){

            ann += 1.0;
            re = sqrt(4.0 * ann * ann - 1.0) / ann;
            aa = re * (xx * amn - em * amm);
            amm = amn;
            amn = aa;
            em = 1.0 / re;
        }

        dd = gg * em * amm - xx * ann * aa;
        dlt = (1.0 - xx * xx) * aa / dd;
        if(fabs(dlt) < LTOL) {
           *wght = gg * (1.0 - xx * xx) / (dd * dd);
           break;
        }

        xx = xx - dlt;

    }

    if(nit >= NIT){

       printf("****ERROR****, no convergence for Lengendre nodes.\n\n");
       exit(1);

    }

    return xx;

}

/* function to compute the normalized associated Legendre functions for triangular truncation */

void pnm(double *aln, double *dln, double xx, int m, int ider)

{

    int i, j;
    int ncof=0;
    int m1 = m + 1;

    int nn2, nn21, n4n2, n4n21;

    long int m2, n2, n21;

    double difx;
    double sqfpi;

    double qmm, qmm1, qqmm, qqmm1;
    double fct1, fct2, fct, fcc;
    double tm3, c1, c2, qm2;

    double *pnml=NULL;
    double *a1=NULL, *a2=NULL;
    double *tmp=NULL, *tmp2=NULL;

    tmp = aln;
    tmp2 = dln;

    sqfpi = 1.0 / sqrt(4.0 * FPI);

    difx = -sqrt(1.0 - xx * xx) / 2.0;
    fct1 = 1.0;
    fct2 = 0.0;
    fct = 1.0;
    fcc = 1.0;

    qmm = 1.0;

    *(aln++) = sqfpi * qmm;
    if(ider) {
      pnml = (double *)calloc(m1, sizeof(double));
      mem_er((pnml == NULL) ? 0 : 1, m1*sizeof(double));
    }

    ++ncof;

    for(i=0; i < m1; i++){

        tm3 = sqrt((double)(2 * i + 3));

        if(i > 0){

           fct1 += 2.0;
           fct2 += 2.0;
           fct = fct1 * fct2;

           qmm *= sqrt(fct) * difx / (double) i;

           *(aln++) = sqfpi * qmm;
           ++ncof;

        }


        if(i == m) {
          if(ider) *(pnml + m) = sqfpi * xx * tm3 * qmm;
          continue;
        }

        qmm1 = xx * tm3 * qmm;

        *(aln++) = sqfpi * qmm1;
        ++ncof;

        m2 = i * i;

        qqmm = qmm;
        qqmm1 = qmm1;

        for(j=i+2; j < m1; j++){

            n2 = j * j;
            n21 = (j-1) * (j - 1);
            c1 = sqrt((double)(4 * n2 - 1) / (double)(n2 - m2));
            c2 = sqrt((double)(n21 - m2) / (double)(4 * n21 - 1));


            qm2 = c1 * (xx * qqmm1 - c2 * qqmm);
            qqmm = qqmm1;
            qqmm1 = qm2;
            *(aln++) = sqfpi * qm2;
            ++ncof;

        }

        if(ider){
            n2 = m1 * m1;
            n21 = (m1-1) * (m1 - 1);
            c1 = sqrt((double)(4 * n2 - 1) / (double)(n2 - m2));
            c2 = sqrt((double)(n21 - m2) / (double)(4 * n21 - 1));

           *(pnml + i) = sqfpi * c1 * (xx * qqmm1 - c2 * qqmm);
        }
        
      
    }

/* compute derivatives wrt to latitude*/

    if(ider){

       aln = tmp;

       a1 = aln;
       a2 = aln + 2;

       nn2 = m * m;
       nn21 = m1 * m1;
       n4n2 = (double)(4 * nn2 - 1);
       n4n21 = (double)(4 * nn21 - 1);

       *(dln++) = 0.0;

       for(i=0; i < m1; i++){

          m2 = i * i;

          if(i == m)continue;

          if(i > 0){

             n2 = (i + 1) * (i + 1);
             c1 = sqrt((double)(n2 - m2) / (double)(4 * n2 - 1));
             *(dln++) = - ((double) i) * c1 * *(a2++);
             ++a1;

          }

          for(j = i+1; j < m; j++){
              n2 = j * j;
              n21 = (j+1) * (j + 1);

              c1 = sqrt((double)(n2 - m2) / (double)(4 * n2 - 1));
              c2 = sqrt((double)(n21 - m2) / (double)(4 * n21 - 1));
             *(dln++) = ((double) ( j + 1)) * c1 * *(a1++) - ((double) j) * c2 * *(a2++);
          
          }

          c1 = sqrt((double)(nn2 - m2) / n4n2);
          c2 = sqrt((double)(nn21 - m2) / n4n21);
         *(dln++) = ((double) m1) * c1 * *(a1++) - ((double) m) * c2 * *(pnml + i);
         ++a2;
          
       }

       c2 = sqrt((double)(nn21 - nn2) / n4n21);
       *dln =  - ((double) m) * c2 * *(pnml + m);
       
    }

    return;

}
