#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include <string.h>
#include "sqt.h"
#include "mem_er.h"

/* function to generate the SQT hierarchy  */

double area_sp_triangle(SQT * , VEC * , int );

VEC *create_sqt(SQT **sqth, VEC *gv, int nfc, int nl, int *nvec, int mxnl)
{

    int i, j, k;
    int nelm=nfc;
    int nt;
    int ifnd=0;
    int imap[3];
    int infc=0;
    int llev=0;

    SQT *sq0=NULL, *sqtmp=NULL, *sq1t=NULL;
    LEAF *lf=NULL;
    VEC pp[3];
    VEC *pa=NULL, *pb=NULL;

    double nrm;

    if(nl == (mxnl - 1)) llev = 1;

    nelm *= 4;
    sq0 = sqth[nl-1];

    for(i=0; i<nfc; i++){

        sqtmp = sq0 + i;

        pa = gv + sqtmp->ivec[0];
        pb = gv + sqtmp->ivec[1];

        addv(pa, pb, &pp[0]);
        nrm = sqrt(dotp(&pp[0], &pp[0]));
        normv(&pp[0], nrm);

        pa = gv + sqtmp->ivec[1];
        pb = gv + sqtmp->ivec[2];

        addv(pa, pb, &pp[1]);
        nrm = sqrt(dotp(&pp[1], &pp[1]));
        normv(&pp[1], nrm);

        pa = gv + sqtmp->ivec[2];
        pb = gv + sqtmp->ivec[0];

        addv(pa, pb, &pp[2]);     
        nrm = sqrt(dotp(&pp[2], &pp[2]));
        normv(&pp[2], nrm);

        nt = *nvec;

        for(k=0; k<3; k++){
            ifnd = 0;
            for(j=0; j< *nvec; j++){
               if(fabs(dotp(&pp[k], gv+j) - 1.0) < TOLVEC){
                  imap[k] = j;
                  ifnd = 1;
                  break;
               }
            }
            if(!ifnd){
               imap[k] = nt;  
               ++nt;
               gv = (VEC *)realloc_n(gv, nt * sizeof(VEC));
               mem_er((gv == NULL) ? 0 : 1, nt * sizeof(VEC));
               *(gv + nt - 1) = pp[k];
            }

        }


        *nvec = nt;


        if(!llev) {
          sq1t = sqth[nl] + infc;
          sqtmp->tt[0] = sq1t;
          sq1t->ud = (sqtmp->ud == 'u') ? 'u' : 'd';
          sq1t->parent = sqtmp;
          sq1t->ivec[0] = sqtmp->ivec[0];
          sq1t->ivec[1] = imap[0];
          sq1t->ivec[2] = imap[2];
          ++infc;

          sq1t = sqth[nl] + infc;
          sqtmp->tt[1] = sq1t;
          sq1t->ud = (sqtmp->ud == 'u') ? 'u' : 'd';
          sq1t->parent = sqtmp;
          sq1t->ivec[0] = imap[0];
          sq1t->ivec[1] = sqtmp->ivec[1];
          sq1t->ivec[2] = imap[1];
          ++infc; 

          sq1t = sqth[nl] + infc;
          sqtmp->tt[2] = sq1t;
          sq1t->ud = (sqtmp->ud == 'u') ? 'd' : 'u';
          sq1t->parent = sqtmp;
          sq1t->ivec[0] = imap[1];
          sq1t->ivec[1] = imap[0];
          sq1t->ivec[2] = imap[2];
          ++infc;

          sq1t = sqth[nl] + infc;
          sqtmp->tt[3] = sq1t;
          sq1t->ud = (sqtmp->ud == 'u') ? 'u' : 'd';
          sq1t->parent = sqtmp;
          sq1t->ivec[0] = imap[2];
          sq1t->ivec[1] = imap[1];
          sq1t->ivec[2] = sqtmp->ivec[2];
          ++infc;

        }
        else {
          lf = (LEAF *)sqth[nl] + infc;
          sqtmp->tt[0] = (SQT *)lf;
          lf->ud = (sqtmp->ud == 'u') ? 'u' : 'd';
          lf->parent = sqtmp;
          lf->ivec[0] = sqtmp->ivec[0];
          lf->ivec[1] = imap[0];
          lf->ivec[2] = imap[2];
          lf->area = area_sp_triangle((SQT *)lf, gv, 1);
          ++infc;

          lf = (LEAF *)sqth[nl] + infc;
          sqtmp->tt[1] = (SQT *)lf;
          lf->ud = (sqtmp->ud == 'u') ? 'u' : 'd';
          lf->parent = sqtmp;
          lf->ivec[0] = imap[0];
          lf->ivec[1] = sqtmp->ivec[1];
          lf->ivec[2] = imap[1];
          lf->area = area_sp_triangle((SQT *)lf, gv, 1);
          ++infc; 

          lf = (LEAF *)sqth[nl] + infc;
          sqtmp->tt[2] = (SQT *)lf;
          lf->ud = (sqtmp->ud == 'u') ? 'd' : 'u';
          lf->parent = sqtmp;
          lf->ivec[0] = imap[1];
          lf->ivec[1] = imap[0];
          lf->ivec[2] = imap[2];
          lf->area = area_sp_triangle((SQT *)lf, gv, 1);
          ++infc;

          lf = (LEAF *)sqth[nl] + infc;
          sqtmp->tt[3] = (SQT *)lf;
          lf->ud = (sqtmp->ud == 'u') ? 'u' : 'd';
          lf->parent = sqtmp;
          lf->ivec[0] = imap[2];
          lf->ivec[1] = imap[1];
          lf->ivec[2] = sqtmp->ivec[2];
          lf->area = area_sp_triangle((SQT *)lf, gv, 1);
          ++infc;
        }

    }

    return gv;

}
