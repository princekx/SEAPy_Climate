#include <Stdio.h>
#include <Math.h>
#include "m_values.h"

/* evaluation of the normalization of the non-isotropic power kernal function. */

#define  TOLINT   1.0e-12
#define  LARGEN   1.0

double power_int_mn(double , int , int );
double legendre_p(int , double , int );

double power_norm(double kk, double beta, int m)

{

    int n=0;

    double alpha, alpha2;
    double aa=1.0, m1=1.0;
    double sum=0.0;
    double dinc=LARGEN;

    alpha = 2.0 * beta / kk;
    alpha2 = alpha * alpha;

    while(fabs(dinc) > TOLINT){

        dinc = aa * m1 * legendre_p(n, 0.0, 0) * power_int_mn(kk, m, n);

        sum += dinc;

        aa *= alpha2;
        m1 *= -1.0;
        n += 2;


    }

    return sum;

}
