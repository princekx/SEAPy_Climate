#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include "statistic.h"
#include "sqt.h"
#include "mem_er.h"
#include "m_values.h"

int inside_sqt(VEC * , VEC * , VEC * , VEC * , char );

void data_partition(SQT **sq, VEC *gv, int nfc, int nl, int dtype, struct dpt *dt , int ndt)
{

    int i, j, k, m;

    struct dpt *dtt=NULL;

    VEC pp, v1, v2, v3;
    SQT *sqt=NULL, *sqt1=NULL;
    LEAF *lf;

    for(i=0; i < ndt; i++){

       dtt = dt + i;
       pp.x = dtt->xdt;
       pp.y = dtt->ydt;
       pp.z = dtt->zdt;


       for(j=0; j < nfc; j++){

           sqt = *sq + j;
           v1 = *(gv + sqt->ivec[0]);
           v2 = *(gv + sqt->ivec[1]);
           v3 = *(gv + sqt->ivec[2]);

           if(inside_sqt(&pp, &v1, &v2, &v3, sqt->ud)){

              for(k=1; k < nl-1; k++){

                  for(m=0; m < 4; m++){
                     sqt1 = sqt->tt[m];
                     v1 = *(gv + sqt1->ivec[0]);
                     v2 = *(gv + sqt1->ivec[1]);
                     v3 = *(gv + sqt1->ivec[2]);
                     if(inside_sqt(&pp, &v1, &v2, &v3, sqt1->ud)) break;
                  }

                  sqt = sqt1;

              }

              for(m=0; m < 4; m++){
                 lf = (LEAF *)(sqt->tt[m]);

                 v1 = *(gv + lf->ivec[0]);
                 v2 = *(gv + lf->ivec[1]);
                 v3 = *(gv + lf->ivec[2]);

                 if(inside_sqt(&pp, &v1, &v2, &v3, lf->ud)) {

                    if(!dtype){

                       ++(lf->ng);

                       if(!(lf->lgrid)){
                         lf->lgrid = (int *)malloc_initl(sizeof(int));
                         mem_er((lf->lgrid == NULL) ? 0 : 1, sizeof(int)); 
                       }

                       else{

                         lf->lgrid = (int *)realloc_n(lf->lgrid, (lf->ng) * sizeof(int));
                         mem_er((lf->lgrid == NULL) ? 0 : 1, (lf->ng) * sizeof(int)); 
                       }

                       *(lf->lgrid + lf->ng - 1) = i;
                       

                    }

                    else {

                       ++(lf->ndata);

                       if(!(lf->ldata)){
                         lf->ldata = (int *)malloc_initl(sizeof(int));
                         mem_er((lf->ldata == NULL) ? 0 : 1, sizeof(int)); 
                       }

                       else{

                         lf->ldata = (int *)realloc_n(lf->ldata, (lf->ndata) * sizeof(int));
                         mem_er((lf->ldata == NULL) ? 0 : 1, (lf->ndata) * sizeof(int)); 
                       }

                       *(lf->ldata + lf->ndata - 1) = i;

                    }

                    break;

                 }


              }
          
              break;
           }

       }

    }

    return;

}
