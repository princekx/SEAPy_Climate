#include <Stdio.h>
#include <stdlib.h>


/* function to evaluate part of the normalization factor of a non-isotropic 
   kernel function (power).                                                  */

double  legendre_p(int , double , int );


double power_int_mn(double kk, int mm, int nn)

{

    int i;
    int n2=2*nn;
    int ni[10];

    double ff=0.;
    double kki=1.0/kk;
    double km=kk;

    for(i=0; i < 2*mm + 1; i++) ni[i] = n2 - 2 * mm + 1 + 2 * i;
    for(i=0; i < mm-1; i++) km *= kk;


    switch(mm){
       case 0:
         ff = (legendre_p(nn-1, kki, 0) - legendre_p(nn+1, kki, 1)) / (float)ni[0];
         break;
       case 1:
         ff += legendre_p(nn-2, kki, 0) / (float)(ni[1] * ni[0]);
         ff -= 2.0 * legendre_p(nn, kki, 1) / (float)(ni[0] * ni[2]);
         ff += legendre_p(nn+2, kki, 1) / (float)(ni[1] * ni[2]);
         ff *= km;
         break;
       case 2:
         ff += legendre_p(nn-3, kki, 0) / (float)(ni[0] * ni[1]);
         ff -= 3.0 * legendre_p(nn-1, kki, 1) / (float)(ni[0] * ni[3]);
         ff += 3.0 * legendre_p(nn+1, kki, 1) / (float)(ni[1] * ni[4]);
         ff -= legendre_p(nn+3, kki, 1) / (float)(ni[3] * ni[4]);
         ff *= 2.0 * km / ni[2];
         break;
       case 3:
/* need to check this case */
         ff += legendre_p(nn-4, kki, 0) / (float)(ni[0] * ni[1] * ni[2] * ni[3]);
         ff -= 4.0 * legendre_p(nn-2, kki, 1) / (float)(ni[0] * ni[2] * ni[3] * ni[4]);
         ff += 6.0 * legendre_p(nn, kki, 1) / (float)(ni[1] * ni[2] * ni[4] * ni[5]);
         ff -= 4.0 * legendre_p(nn+2, kki, 1) / (float)(ni[2] * ni[3] * ni[4] * ni[6]);
         ff += legendre_p(nn+4, kki, 1) / (float)(ni[3] * ni[4] * ni[5] * ni[6]);
         ff *= 6.0 * km;
         break;
       case 4:
         ff += legendre_p(nn-5, kki, 0) / (float)(ni[0] * ni[1] * ni[2] * ni[3]);
         ff -= 5.0 * legendre_p(nn-3, kki, 1) / (float)(ni[0] * ni[2] * ni[3] * ni[5]);
         ff += 10.0 * legendre_p(nn-1, kki, 1) / (float)(ni[1] * ni[2] * ni[5] * ni[6]);
         ff -= 10.0 * legendre_p(nn+1, kki, 1) / (float)(ni[2] * ni[3] * ni[6] * ni[7]);
         ff += 5.0 * legendre_p(nn+3, kki, 1) / (float)(ni[3] * ni[5] * ni[6] * ni[8]);
         ff -= legendre_p(nn+5, kki, 1) / (float)(ni[5] * ni[6] * ni[7] * ni[8]);
         ff *= 24.0 * km / ni[4];
         break;
       default:
         printf("***ERROR***, no evaluation of integral for this index\r\n");
         printf("Error occured in %s\n\n", __FILE__);
         exit(1);

   } 
       




    return ff;

}


