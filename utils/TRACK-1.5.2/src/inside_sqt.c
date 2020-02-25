#include <Stdio.h>
#include "sqt.h"


int inside_sqt(VEC *pt, VEC *v1, VEC *v2, VEC *v3, char ior)

{

   int inn=0;

   double d1, d2, d3;

   VEC crs;

   if(ior == 'd') crosp(v1, v3, &crs)
   else crosp(v1, v2, &crs)
   d1 = dotp(pt, &crs);
   
   if(ior == 'd') {crosp(v3, v2, &crs);}
   else {crosp(v2, v3, &crs);}
   d2 = dotp(pt, &crs);

   if(ior == 'd') {crosp(v2, v1, &crs);}
   else {crosp(v3, v1, &crs);}
   d3 = dotp(pt, &crs);

   if(d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0) inn = 1;

   return inn;

}
