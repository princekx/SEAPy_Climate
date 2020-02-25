#include <Stdio.h>
#include <Math.h>
#include "statistic.h"
#include "m_values.h"
#include "p_vecs.h"

/* evaluate the sample statistics for a linear kernal function */

double variance_non_iso(double , struct dpt * , struct dpt * , double * , float * , int , int , struct cvecs * , int , double );
double non_iso_density(struct dpt * , struct dpt * , struct cvecs * ,double , double , int , int * );

double power_norm(double , double , int );

double eval_non_iso_stat(double *plt, double **den, struct dpt *st, struct dpt *dt, float *wght, int dtn, int ptnum, int im, int ht, int ind, int ms, struct cvecs *pcom, int nni)

{

    int i, j;
    int ireg=0;

    float cf0, cf1[2];

    static int iv=0;

    double nn=0.0, tsm=1., ts;
    double mean=0.0, dtmp=0.0, dent=0.0;

    double cn=0.0, con=0., cnn, dd;
    double *d[5]={NULL,NULL,NULL,NULL,NULL};

    double beta;

    double ddd=0.0;
    double wt=0.0;

    struct dpt *dtt, *spt;
    struct cvecs *pp=NULL;

    if(ht != 3){

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


    if(!wght){

       nn = (float) dtn;
       cn = 1.0 / (FPI2 * nn);
    }
    else{
       nn = 0.0;
       for(j=0; j < dtn; j++) nn += *(wght + j);
       cn = (nn > TOLWT) ? 1.0 / (FPI2 * nn) : 0.0;

    }

    con = cn;

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

       for(j=0; j < dtn; j++){ 

          ts = *(plt + j);

          tsm = (1.0 + ts) / ts;

          pp = pcom + j;

          beta = pp->beta * tsm;

          dtt = dt + j;

          if(wght)wt = *(wght + j);
          else wt = 1.0;

          ddd = non_iso_density(spt, dtt, pp, tsm, beta, nni, &ireg);

          cnn = (ireg) ? 1.0 / power_norm(tsm, beta, nni) : 0.0; 


          *d[0] += (dd = wt * (dtmp = cnn * ddd));
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


       *d[0] *= con;
       dent *= con; 

       if(d[1]) {

          *d[1] = (*d[0] > DTOL) ? *d[1] * con / *d[0] : 0.0;
           mean = (dent > DTOL) ? mean * con / dent : 0.0;
          if(iv) *d[2] = sqrt(variance_non_iso(mean, spt, dt, plt, wght, dtn, ht, pcom, nni, nn));
       }
       if(d[3]){

          if(*d[0] > DTOL){*d[3] *= con / *d[0]; *d[4] *= con / *d[0];}
          else { *d[3] = 0.0; *d[4] = 0.0;}

       }

    }

    return nn;

}



/* function to estimate the variance at a point */

double non_iso_density(struct dpt * , struct dpt * , struct cvecs * ,double , double , int , int * );

double power_norm(double , double , int );

double variance_non_iso(double mean, struct dpt *spt, struct dpt *dt, double *plt, float *wght, int dtn, int ht, struct cvecs *pcom, int nni, double nn)

{

    int j;
    int ireg=0;

    double var=0.;
    double den=0., dif;
    double tsm=0., ts;
    double beta;
    double ddd=0.0;
    double wt=0.0;

    double cn, con=0., cnn, dd, cd;
    struct dpt *dtt;
    struct cvecs *pp=NULL;

    if(ht != 3){

      printf("***error***, non-valid switch for density estimation type in %s \n\n", __FILE__);
      return var;

    }


    if(nn > TOLWT) cn = 1.0 / (FPI2 * nn);
    else cn = 0.0;

    con = cn;

    for(j=0; j < dtn; j++){

       ts = *(plt + j);


       tsm = (1.0 + ts) / ts;

       cd = 1.0 - tsm;

       pp = pcom + j;

       beta = pp->beta * tsm;

       dtt = dt + j;

       if(wght){wt = *(wght + j); nn += wt;}
       else wt = 1.0;

       ddd = non_iso_density(spt, dtt, pp, tsm, beta, nni, &ireg);

       cnn = (ireg) ? 1.0 / power_norm(tsm, beta, nni) : 0.0;

       den += (dd = wt * cnn * ddd);

       dif = dtt->sdt - mean;

       var += dif * dif * dd;

    }


    den *= con;

    var = (den > DTOL) ? var * con / den : 0.0;

    return var;

}
