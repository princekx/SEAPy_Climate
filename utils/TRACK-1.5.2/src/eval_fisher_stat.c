#include <Stdio.h>
#include <Math.h>
#include "statistic.h"
#include "m_values.h"

/* evaluate the sample statistics for a fisher kernal function */

#define   ETOL   -1.0e+100

double variance_fisher(double , struct dpt * , struct dpt * , double * , float * , int , int , double );
double fisher_density(struct dpt * , struct dpt * , double ,double , int );

double eval_fisher_stat(double *plt, double **den, struct dpt *st, struct dpt *dt, float *wght, int dtn, int ptnum, int im, int ht, int ind, int ms)

{

    int i, j;
    int nf=0;

    float cf0, cf1[2];

    static int iv=0;

    double nn=0.0, ex, ex1, tsm=0.;

    double cn=0.0, con=0., cnn, dd, lgtol;
    double *d[5]={NULL,NULL,NULL,NULL, NULL};
    double wt=0.0;
    double mean=0.0, dtmp=0.0, dent=0.0;

    struct dpt *dtt, *spt;

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


    lgtol = ETOL;

    if(!wght){
       nn = (float) dtn;
       cn = 1.0 / (2.0 * FPI * nn);

    }

    else{
       nn = 0.0;
       for(j=0; j < dtn; j++) nn += *(wght + j);
       cn = (nn > TOLWT) ? 1.0 / (2.0 * FPI * nn) : 0.0;
    }

    if(ht < 3){

      tsm = *plt;

      if(tsm < SMMIN) con = 0.5 * cn;

      else {

         if(tsm > DBAND) nf = 1;  /* change evaluation of exponetial */

/* compute c/4*pi*n*sinh(c) for use with fisher density */

         ex = exp(tsm);
         ex1 = 1.0 / ex;

         cn *= tsm;

         con = (nf) ? cn : cn / (ex - ex1);

      }

    }

    else if(ht == 3){

      con = cn;

    }
   

    if(im > 1 || im == 0){

      if(ind){


          printf("do you wish to estimate the standard deviation as well as the mean??\r\n"
                 "Input '1' for yes or '0' for no\n\n");
 
          scanf("%d", &iv);


      }

    }

    for(i=0; i < ptnum; i++){

       mean = 0.0;
       dent = 0.0;

       d[0] = den[0] + i;
       *d[0] = 0.0;
       spt = st + i;

       if(im > 1){

         d[1] = den[1] + i;
         d[2] = den[2] + i;
         if(im == 3){ d[3] = den[3] + i; d[4] = den[4] + i;}

       }
       else if(im <= 0){

         if(im == 0) {d[1] = den[1] + i; d[2] = den[2] + i;}
         else { d[3] = den[3] + i; d[4] = den[4] + i;}

       }

       if(d[1]) *d[1] = *d[2] = 0.0;
       if(d[3]) *d[3] = *d[4] = 0.0;


       if(ht < 3){

         for(j=0; j < dtn; j++){

            dtt = dt + j;

            if(wght)wt = *(wght + j);
            else wt = 1.0;


            *d[0] += (dd = wt * (dtmp = fisher_density(spt, dtt, tsm, lgtol, nf)));
            dent += dtmp;

            if(d[1]){

              if(ms && fabs(dtt->sdt - dtt->dsdt) / cf0 > 1.) dtt->sdt = dtt->dsdt;

              *d[1] += dtt->sdt * dd;
              mean += dtt->sdt * dtmp;

            }

            if(d[3]) {

              if(ms && fabs(dtt->vec[0] - dtt->dvec[0]) / cf1[0] > 1.) dtt->vec[0] = dtt->dvec[0];
              if(ms && fabs(dtt->vec[1] - dtt->dvec[1]) / cf1[1] > 1.) dtt->vec[1] = dtt->dvec[1];

              *d[3] += dtt->vec[0] * dd; 
              *d[4] += dtt->vec[1] * dd;

            }
     
         }

       }

       else if(ht == 3) {

         for(j=0; j < dtn; j++){

            nf = 0;

            tsm = *(plt + j);

            if(tsm < SMMIN) cnn = 0.5 * cn;

            else {

               if(tsm > DBAND) {    /* change evaluation of exponetial */

                  nf = 1;
                  cnn = tsm;

               }

               else{

                  ex = exp(tsm);
                  ex1 = 1.0 / ex;
                  cnn = tsm / (ex - ex1);

               }

            }


            dtt = dt + j;

            if(wght)wt = *(wght + j);
            else wt = 1.0;

            *d[0] += (dd = wt * (dtmp = cnn * fisher_density(spt, dtt, tsm, lgtol, nf)));
            dent += dtmp;

            if(d[1]) {

              if(ms && fabs(dtt->sdt - dtt->dsdt) / cf0 > 1.) dtt->sdt = dtt->dsdt;

              *d[1] += dtt->sdt * dd;
              mean += dtt->sdt * dtmp;

            }

            if(d[3]) {

              if(ms && fabs(dtt->vec[0] - dtt->dvec[0]) / cf1[0] > 1.) dtt->vec[0] = dtt->dvec[0];
              if(ms && fabs(dtt->vec[1] - dtt->dvec[1]) / cf1[1] > 1.) dtt->vec[1] = dtt->dvec[1];

               *d[3] += dtt->vec[0] * dd; 
               *d[4] += dtt->vec[1] * dd;

            }    

         }

       }



       *d[0] *= con;
       dent *= con;

       if(d[1]) {
          *d[1] = (*d[0] > DTOL) ? *d[1] * con / *d[0] : 0.0;
          mean = (dent > DTOL) ? mean * con / dent : 0.0;
          if(iv) *d[2] = sqrt(variance_fisher(mean, spt, dt, plt, wght, dtn, ht, nn));
       }
       if(d[3]){
          if(*d[0] > DTOL){*d[3] *= con / *d[0]; *d[4] *= con / *d[0];}
          else { *d[3] = 0.0; *d[4] = 0.0;}
       }


    }

    return nn;

}

/* function to estimate the variance at a point */

double fisher_density(struct dpt * , struct dpt * , double ,double , int );

double variance_fisher(double mean, struct dpt *spt, struct dpt *dt, double *plt, float *wght, int dtn, int ht, double nn)

{

    int j;
    int nf=0;

    double var=0.;
    double ex, ex1, tsm=0.;
    double wt=0.0;

    double cn, con=0., cnn, dd, lgtol;

    double den=0., dif;

    struct dpt *dtt;



    lgtol = ETOL;

    if(nn > TOLWT) cn = 1.0 / (2.0 * FPI * nn);
    else cn = 0.0;


    if(ht < 3){

      tsm = *plt;

      if(tsm < SMMIN) con = 0.5 * cn;

      else {

         if(tsm > DBAND) nf = 1;  /* change evaluation of exponetial */

/* compute c/4*pi*n*sinh(c) for use with fisher density */

         ex = exp(tsm);
         ex1 = 1.0 / ex;

         cn *= tsm;

         con = (nf) ? cn : cn / (ex - ex1);

      }

      for(j=0; j < dtn; j++){

         dtt = dt + j;

         if(wght)wt = *(wght + j);
         else wt = 1.0;

         den += (dd = wt * fisher_density(spt, dtt, tsm, lgtol, nf));

         dif = dtt->sdt - mean;

         var += dif * dif * dd;
     
      }


    }

    else if(ht == 3){

      con = cn;


      for(j=0; j < dtn; j++){

         nf = 0;

         tsm = *(plt + j);

         if(tsm < SMMIN) cnn = 0.5 * cn;

         else {

            if(tsm > DBAND) {    /* change evaluation of exponetial */

               nf = 1;
               cnn = tsm;

            }

            else{

               ex = exp(tsm);
               ex1 = 1.0 / ex;
               cnn = tsm / (ex - ex1);

            }

         }


         dtt = dt + j;

         if(wght)wt = *(wght + j);
         else wt = 1.0;

         den += (dd = wt * cnn * fisher_density(spt, dtt, tsm, lgtol, nf));

         dif = dtt->sdt - mean;

         var += dif * dif * dd;    


      }


    }



    den *= con;

    var = (den > DTOL) ? var * con / den : 0.0;


    return var;

}
