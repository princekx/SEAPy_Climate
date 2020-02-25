#include <Stdio.h>
#include <Math.h>
#include "p_vecs.h"
#include "statistic.h"

#define DENSS  0.0

void householder(double **, int , double * , double * );
void ql(double * , double * , int , double ** );
void eigsort(double *, double **, int );


/* function to estimate principle component eigen-vectors and eigen-values
   of a scatter matrix for use in non-isotropic kernel functions.          */


void estimate_vectors(double (*density)(), double dnorm, struct cvecs *pp, double *plt, float *wght, double tsm, int np)

{


    int i, j,k,l;

    double mm[3][3], yy[3], sc[3][3], *a[3];
    double d[3], e[3];
    double *cp=NULL, *cpp=NULL;
    double con, ww;
    double dens;
    double sum=0.;
    double wt;

    struct dpt pp1, pp2;

    for(i=0; i<3; i++) a[i] = &sc[i][0];

    con = dnorm;


    for(i=0; i < np; i++){

        for(j=0; j < 3; j++) {

           for(k=0; k < 3; k++)sc[j][k] = 0.0;
           d[j] = e[j] = 0.0;

        }

        cp = pp[i].p1;
        *(cp+3) = 0.0;

        pp1.xdt = *cp;
        pp1.ydt = *(cp + 1);
        pp1.zdt = *(cp + 2);

/* Compute projection operator */

        for(j=0; j < 3; j++){

            for(k=0; k < 3; k++){

                mm[j][k] = ((j == k) ? 1.0 : 0.0 ) - *(cp + j) * *(cp + k);
 
            }

        }



        dens = DENSS;

/* Estimate Scatter matrix */

        for(j=0; j < np; j++){

            cpp = pp[j].p1;

            pp2.xdt = *cpp;
            pp2.ydt = *(cpp + 1);
            pp2.zdt = *(cpp + 2);

            for(k=0; k < 3; k++){

               yy[k] = 0.0;

               for(l=0; l < 3; l++) yy[k] += mm[k][l] * *(cpp + l);

            }

            ww = (*density)(&pp1, &pp2, tsm);

            wt = (wght) ? *(wght + j) : 1.0;

            dens += wt * ww;

/* Compute contribution to scatter matrix */

            for(k=0; k < 3; k++){

                for(l=0; l < 3; l++){

                    sc[k][l] += yy[k] * yy[l] * ww *wt;

                }


            }

        }

        dens *= con;

        if(plt) *(plt + i) = dens;

        for(k=0; k < 3; k++){

            for(l=0; l < 3; l++){

                sc[k][l] *= con;

            }


        }


/* Scale scatter matrix with density */


        if(dens > DTOL){

           for(j=0; j<3; j++) {
               for(k=0; k<3; k++) sc[j][k] /= dens; 
           } 

        }

        else {pp[i].p1[3] = pp[i].p2[3] = pp[i].p3[3] =0.0; continue;}



/* Compute eigenvalues and eigenvectors */


        householder(a, 3, d, e);
        ql(d, e, 3, a);
        eigsort(d, a, 3);

/* sort eigenvalues and eigenvectors for a pair of zero degenerate eigenvalues */

        if(fabs(d[1]) < EIGTOL && fabs(d[2]) < EIGTOL){

           cp = pp[i].p1;

           sum = 0.;

           for(j=0; j<3; j++) sum += *(cp + j) * sc[j][1]; 


           if(fabs(sqrt(sum) - 1.0) < EIGTOL) {

              for(j=0; j<3; j++){

                 sum = sc[j][1];
                 sc[j][1] = sc[j][2];
                 sc[j][2] = sum;

              }

           }

        }
        
        cp = pp[i].p3;

        for(j=0; j<3; j++) *(cp+j) = sc[j][0];
        *(cp+3) = d[0];

        cp = pp[i].p2;
        for(j=0; j<3; j++) *(cp+j) = sc[j][1];
        *(cp+3) = d[1]; 



    }


    return;

}
