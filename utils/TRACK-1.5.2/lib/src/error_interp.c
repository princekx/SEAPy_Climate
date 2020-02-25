#include <Stdio.h>
#include <stdlib.h>

/* function to report errors from smoopy.f */

void error_interp(int ier, char *fi, int line, int tmin)

{

    printf("       EEEEE  RRR   RRR      OO     RRR         \r\n"
           " ***** E      R  R  R  R   OO  OO   R  R *****  \r\n"
           " ***** EEE    RRR   RRR   O      O  RRR  *****  \r\n"
           " ***** E      R R   R R    OO  OO   R R  *****  \r\n"
           "       EEEEE  R  R  R  R     OO     R  R        \n\n"); 

    printf("the following report indicates a problem has been detected   \r\n"
           "with the surface fitting routine being used.                 \r\n"
           "Error occured in File %s at or near line %d\n\n", fi, line);

    switch(ier){

       case 1:

         printf("the required storage space exceeds the available space, \r\n"
                "specified by the parameters nxest and nyest( data init.)\r\n"
                "probably causes: s,nxest or nyest too small.              \n");
         break;

       case 2:

         printf("a theoretically impossible behaviour of the function    \r\n"
                "f(p) was found during the iteration process.            \r\n"
                "probably causes: tol too small.                           \n");
         break;

       case 3:

         printf("the maximum allowable number of iterations to find the  \r\n"
                "root of f(p)=s has been reached.                        \r\n"
                "probably causes: maxit or tol too small.                  \n");
         break;

       case 10:

         printf("some of the input data are invalid(see restrictions)      \n");

         break;

       case 20:

         printf("incompatable parameters kx and ky defined in            \r\n"
                "param_smpy.ptr and smoopy_setup.c                         \n");

         break;

       case 100:

         printf("incompatable parameters kx and ky defined for use in      \r\n"
                "sphery.f (see sphery.f) only bi-cubic available.          \n");

         break;

     }

     printf("\n");
     printf("*****ABNORMAL EXIT*****\n\n");
     printf("ERROR CODE = %d\n", ier);

     if(tmin) exit(1);

     return;

}
