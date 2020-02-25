#include <Stdio.h>
#include <Math.h>

/* function to computed the sulution of an algebraic expression
   using the Secant method.                                     */

#define   TOLSEC    0.00001
#define   TOLFUN    0.000001
#define   MAXIT     30

float secant(float (*ff)(float ), float p1, float p2, float rhs)

{

    int j=0;

    float fl, f, dx=1.0;
    float swap, pl, rts;

    fl = (*ff)(p1) - rhs;
    f = (*ff)(p2) - rhs;

/* check bounds */

    if(fabs(fl) < fabs(f)){

       rts = p1;
       pl = p2;
       swap = fl;
       fl = f;
       f = swap;

     }

     else{rts = p2; pl = p1; }

/* secant loop */

     while(fabs(dx) > TOLSEC){

         dx = (pl - rts) * f / (f - fl);
         pl = rts;
         fl = f;
         rts += dx;
         f = (*ff)(rts) - rhs;

         if(++j > 30 ){

            printf("***ERROR***, Max. iterations exceeded for function secant, current value returned\n");
            break;

         }

         if(fabs(dx) < TOLSEC || fabs(f) <= TOLFUN) break;

      }

      return rts;

}
