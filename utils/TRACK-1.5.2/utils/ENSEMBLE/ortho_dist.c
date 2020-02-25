#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "splice.h"
#include "m_values.h"
#include "sqt.h"


/* function to compute the along and across track seperation distances */

double fvec(VEC * , VEC * , VEC * , VEC * );

double ortho_dist(struct fet_pt_tr *fp1, struct fet_pt_tr *fp2, int nf2, int *ist, int ifnd, int nng, VEC *visec)

{
    int i=0;

    double dis=-1.0;
    double aa, bb, dotab;

    VEC a, b, c;
    VEC r;

    struct fet_pt_tr *ff1=NULL, *ff2=NULL;

    c.x = fp1->pp[0];
    c.y = fp1->pp[1];
    c.z = fp1->pp[2];

    if(ifnd) *ist = 1;
   
    for(i=*ist; i < nf2; i++){

       ff1 = fp2 + i - 1;
       a.x = ff1->pp[0];
       a.y = ff1->pp[1];
       a.z = ff1->pp[2];
       ff2 = fp2 + i;
       b.x = ff2->pp[0];
       b.y = ff2->pp[1];
       b.z = ff2->pp[2];

       dotab = fvec(&a, &b, &c, &r);

       aa = dotp(&r, &a);
       bb = dotp(&r, &b);

       if(dotab <= ((aa < bb) ? aa: bb)){
          *ist = i;
          dis = acos(dotp(&r, &c));
          memcpy(visec, &r, sizeof(VEC));

          break;
       }

       if(nng) *ist = -1;

    }

    return dis;
}

double fvec(VEC *a, VEC *b, VEC *c, VEC *r)
{

    VEC an, bn;

    double aa, bb, dotab;
    double norm=0.0;

    dotab = dotp(a, b);
    aa = dotp(a, c) - dotp(b, c) * dotab;
    bb = dotp(b, c) - dotp(a, c) * dotab;

    mulv(a, &an, aa);
    mulv(b, &bn, bb);

    addv(&an, &bn, r);

    norm = sqrt(dotp(r, r));
    normv(r, norm);

    return dotab;

}
