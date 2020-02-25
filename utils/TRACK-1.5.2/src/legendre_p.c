#include <Stdio.h>

/* function to evaluate the Legendre polynomials of the first kind P_n(z) */


double legendre_p(int nn, double ki, int nw)

{

     int i;

     static int li=0;

     static double pnk;
     static double pn1;
     double pi=0.;

     if(li == 0 || !nw){li=1; pnk = 1.0; pn1 = ki;}

     if(nn < 0) nn = -1 * (nn + 1);

     if(nn == 0)  {li = 0; pnk = 1.0; pn1 = ki; return pnk;}

     else if(nn == 1) {li = 1; pnk = 1.0; pn1 = ki; return pn1;}

     else {

         if(nn >= li){

             if(nn == li) return pn1;

             for(i=li+1; i <= nn; i++){

                pi = (ki * (2.0 * (float)i - 1.0) * pn1 - ((float)i - 1.0) * pnk) / (float) i;
                pnk = pn1;
                pn1 = pi;

             }

             li = nn;

         }

         else {

             if(li - nn < nn){

                if(nn == li-1) return pnk;

                for(i=li-2; i >= nn; i--){

                   pi = (ki * (2.0 * (float) i + 3.0) * pnk - ((float) i + 2.0) * pn1) / ((float) i + 1.0);
                   pn1 = pnk;
                   pnk = pi;

                }

                li = nn+1;

             }

             else {

                pnk = 1.0;
                pn1 = ki;

                for(i=2; i <= nn; i++){

                   pi = (ki * (2.0 * (float)i - 1.0) * pn1 - ((float)i - 1.0) * pnk) / (float) i;
                   pnk = pn1;
                   pn1 = pi;

                }

                li = nn;

             }


         }


     }




     return pi;

}
