#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "m_values.h"
#include "mem_er.h"


struct dpt *grid_convert(struct tot_stat *st)
{

   int i;

   struct pt_stat *pst;
   struct dpt *stat=NULL, *dtt=NULL;

   double sn1, sn2, cn1, cn2;

   stat = (struct dpt *)calloc(st->ptnum, sizeof(struct dpt));
   mem_er((stat == NULL) ? 0 : 1, st->ptnum * sizeof(struct dpt));

   for(i=0; i < st->ptnum; i++){

       pst = st->ptst + i;
       dtt = stat + i;

       dtt->xdt = pst->xs * FP_PI;
       dtt->ydt = FP_PI2 - pst->ys * FP_PI;

       if(dtt->ydt < 0.) dtt->ydt = 0.;

       sincos(dtt->xdt, &sn1, &cn1);
       sincos(dtt->ydt, &sn2, &cn2);

       dtt->xdt = sn2 * cn1;
       dtt->ydt = sn2 * sn1;
       dtt->zdt = cn2;

   }

   return stat;
}
