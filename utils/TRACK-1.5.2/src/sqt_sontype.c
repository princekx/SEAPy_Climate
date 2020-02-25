#include <Stdio.h>
#include <stdlib.h>
#include "sqt.h"

int sqt_sontype(SQT *sq, int ilv, int mxl)
{

    int type=0;

    SQT *sqt=NULL;
    LEAF *lf=NULL;

    if(ilv == mxl){
       lf = (LEAF *)sq;
       sqt = lf->parent;
       if((LEAF *)(sqt->tt[0]) == lf) type = 1;
       else if((LEAF *)(sqt->tt[1]) == lf) type = 2;
       else if((LEAF *)(sqt->tt[2]) == lf) type = 3;
       else if((LEAF *)(sqt->tt[3]) == lf) type = 4;
       else {
         printf("***ERROR***: file \"%s\", line %d\n", __FILE__, __LINE__);
         exit(1);
       }
    }

    else{

       sqt = sq->parent;
       if(sqt->tt[0] == sq) type = 1;
       else if(sqt->tt[1] == sq) type = 2;
       else if(sqt->tt[2] == sq) type = 3;
       else if(sqt->tt[3] == sq) type = 4;
       else {
         printf("***ERROR***: file \"%s\", line %d\n", __FILE__, __LINE__);
         exit(1);
       }

    }


    return type;

}
