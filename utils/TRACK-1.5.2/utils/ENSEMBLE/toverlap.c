#include <stdio.h>
#include "splice.h"

/* Check for overlap in time */

int toverlap(struct tot_tr *tr1, struct tot_tr *tr2, long int *is1, long int *is2)
{

/*    int i,j; */

    long int tr1t1, tr1t2, tr2t1, tr2t2;

    int iover=0;

    if(tr1->time && tr2->time){

       tr1t1 = tr1->trpt->time;
       tr1t2 = (tr1->trpt + tr1->num - 1)->time;

       tr2t1 = tr2->trpt->time;
       tr2t2 = (tr2->trpt + tr2->num - 1)->time;

    }

    else{

       tr1t1 = tr1->trpt->fr_id;
       tr1t2 = (tr1->trpt + tr1->num - 1)->fr_id;

       tr2t1 = tr2->trpt->fr_id;
       tr2t2 = (tr2->trpt + tr2->num - 1)->fr_id;

    }

    if(tr1t1 <= tr2t1 && tr1t2 >= tr2t2) {
      iover = 1;
      *is1 = tr2t1;
      *is2 = tr2t2;
    }
    else if(tr1t1 > tr2t1 && tr1t2 < tr2t2) {
      iover = 1;
      *is1 = tr1t1;
      *is2 = tr1t2;
    }
    else if(tr1t1 <= tr2t1 && tr2t1 <= tr1t2) {
      iover = 1;
      *is1 = tr2t1;
      *is2 = tr1t2;
    }
    else if(tr1t1 <= tr2t2 && tr1t2 >= tr2t2) {
      iover = 1;
      *is1 = tr1t1;
      *is2 = tr2t2;
    }

    return iover;

}
