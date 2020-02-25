#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "m_values.h"
#include "mem_er.h"
#include "sqt.h"
#include "tele.h"

/* evaluate the sample statistics for a linear kernal function */

extern int ipow;

double power_density(struct dpt * , struct dpt * , double , int );
int kernel_insqt(LEAF * , VEC * , struct dpt * , double );
LEAF **sqt_leaf_sample(LEAF ** , struct dpt * , LEAF * , VEC * , double , int * , int * , int );
double tele_kernel(double , double , double , int );

double sqt_eval_stat(double *plt, double **den, struct dpt *st, struct dpt *dt, float *wght, int dtn, int ptnum, int im, int ht, int ind, int ms, LEAF *lf, VEC *gv, int nlf, int sqt_c, TELE *tele)

{

    int i, j, k, m;
    int did, dod;
    int nlff=0;
    int maxlf=0;
    int insqt=0;

    float cf0, cf1[2];

    static int iv=0;

    float ipp;

    double sang=0.0;

    double *d[5]={NULL,NULL,NULL,NULL,NULL};
    double nn=0.0, tsm=1., ts;
    double wt=0.0;
    double *x2=NULL, *mean=NULL, *dent=NULL;
    double dtmp, mss;

    double cn=0.0, con=0., cnn=1.0, dd, cd;

    struct dpt *dt1, *dt2;

    LEAF *lff=NULL, *lft=NULL;
    LEAF **lnf=NULL;

    ipp = ipow + 1;

    if(ht < 0 || ht > 3){

      printf("***error***, non-valid switch for density estimation type in %s \n\n", __FILE__);
      return nn;

    }

    if(ms){

       switch(im){
          case 0: case 2: case 3: 
             printf("What cutoff for the residuals is required for the single regression robust smoother\n\n");
             scanf("%f", &cf0);
             if(im == 0 || im == 2) break;
          case -1:
             printf("What cutoffs for the residuals are required for the vector regression robust smoother in the X and Y directions\n\n");
             scanf("%f %f", &cf1[0], &cf1[1]);
             break;

       }

    }

    if(im > 1 || im == 0){

       if(ind){

          printf("do you wish to estimate the standard deviation as well as the mean??\r\n"
                 "Input '1' for yes or '0' for no\n\n");

          scanf("%d", &iv);

       }

    }
   


    if(!wght){

       nn = (float) dtn;
       cn = ipp / (FPI2 * nn);

    }
    else{
       if(!tele) {
          nn = 0.0;
          for(j=0; j < dtn; j++) nn += *(wght + j);
          cn = (nn > TOLWT) ? ipp / (FPI2 * nn) : 0.0;
       }
       else {
          nn = (float) dtn;
          cn = ipp / (FPI2 * nn);
       }

    }

    if(ht < 3){

      if(*plt < SMMIN) { tsm = -1.0; con = 1.0; sang = -1.0;}

      else {

         tsm = (1.0 + *plt)/ *plt;

         cn *= tsm;

         cd = 1.0 - tsm;

         con = cn / pow(cd, ipp);
         sang = sin(acos(1.0 / tsm));

         cnn = 1.0;

      }


    }

    else if(ht == 3) con = cn;

    if(tele) con /= tele->h;


/* initialize arrays */

    if(im > 1){

      for(i=0; i < ptnum; i++) *(den[1] + i) = 0.0;
      for(i=0; i < ptnum; i++) *(den[2] + i) = 0.0;
      d[1] = den[1];
      d[2] = den[2];
      if(im == 3){ 
        for(i=0; i < ptnum; i++) *(den[3] + i) = 0.0;
        for(i=0; i < ptnum; i++) *(den[4] + i) = 0.0;
        d[3] = den[3]; d[4] = den[4];
      }

    }
    else if(im <= 0){

      if(im == 0) {
         for(i=0; i < ptnum; i++) *(den[1] + i) = 0.0;
         for(i=0; i < ptnum; i++) *(den[2] + i) = 0.0;
         d[1] = den[1]; d[2] = den[2];
      }
      else {
         for(i=0; i < ptnum; i++) *(den[3] + i) = 0.0;
         for(i=0; i < ptnum; i++) *(den[4] + i) = 0.0;
         d[3] = den[3]; d[4] = den[4];
      }

    }

    d[0] = den[0];
    for(i=0; i < ptnum; i++) *(den[0] + i) = 0.0;


    if(iv){
       x2 = (double *)calloc(ptnum, sizeof(double));
       mem_er((x2 == NULL) ? 0 : 1, ptnum * sizeof(double));
    }

    for(i=0; i < nlf; i++) {

       lff = lf + i;

       for(j=0; j < lff->ndata; j++){

           did = *(lff->ldata + j);

           if(wght) {
              wt = (!tele) ? *(wght + did) : tele_kernel(tele->tst, *(wght + did), tele->h, tele->ktyp);

           }
           else wt = 1.0;

           dt1 = dt + did;
           if(d[1] && ms && fabs(dt1->sdt - dt1->dsdt) / cf0 > 1.) dt1->sdt = dt1->dsdt;
           if(d[3] && ms){
              if(fabs(dt1->vec[0] - dt1->dvec[0]) / cf1[0] > 1.) dt1->vec[0] = dt1->dvec[0];
              if(fabs(dt1->vec[1] - dt1->dvec[1]) / cf1[1] > 1.) dt1->vec[1] = dt1->dvec[1];
           }

           nlff = 0;

           if(ht == 3){

              ts = *(plt + did);

              if(ts < SMMIN) {tsm = -1.0; cnn = 1.0; sang = -1.0;}

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

              for(m=0; m < lft->ng; m++){

                  dod = *(lft->lgrid + m);

                  dt2 = st + dod;

                  *(d[0] + dod) += (dd = wt * (dtmp = cnn * power_density(dt1, dt2, tsm, ipow)));

                  if(d[1]) {
                    *(d[1] + dod) += dt1->sdt * dd;
                    if(iv) {
                      *(x2 + dod) += dt1->sdt * dt1->sdt * dd;
                    }

                  }

                  if(d[3]) {
                    *(d[3] + dod) += dt1->vec[0] * dd; 
                    *(d[4] + dod) += dt1->vec[1] * dd;
                  }

              }

           }

       }

    }

    for(i=0; i < ptnum; i++){

       *(d[0] + i) *= con;

    }

    if(d[1]){

       for(i=0; i < ptnum; i++){
           dtmp = *(d[0] + i);
           *(d[1] + i) = (dtmp > DTOL) ? *(d[1] + i) * con / dtmp : 0.0;
       }

       if(iv){

          for(i=0; i < ptnum; i++){
             dtmp = *(d[0] + i);
             mss = *(d[1] + i);
             *(d[2] + i) = (dtmp > DTOL) ? (*(x2 + i) * con / dtmp) - mss * mss : 0.0;
             *(d[2] + i) = (*(d[2] + i) > 0.0) ? sqrt(*(d[2] + i)) : 0.0;
          }

       }

    }

    if(d[3]){
       for(i=0; i < ptnum; i++){
           if(*(d[0] + i) > DTOL){
              *(d[3] + i) *= con / *(d[0] + i);
              *(d[4] + i) *= con / *(d[0] + i);
           }
           else {*(d[3] + i) = *(d[4] + i) = 0.0;}
       }

    }

    free(x2);
    free(mean);
    free(dent);
    free(lnf);

    for(i=0; i < nlf; i++) {
        lff = lf + i;
        free(lff->ldata);
        lff->ndata = 0;
    }

    return nn;

}

