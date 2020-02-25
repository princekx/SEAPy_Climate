#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "sqt.h"
#include "mem_er.h"
#include "m_values.h"

#define  TOLARC  1.0e-5  /* tolerance for tests */ 

/* function to find all neighbouring leaf triangles that overlap with a
   kernel returning a list of leaf triangle pointers.                     */

LEAF **sqt_leaf_sample(LEAF **ilf, struct dpt *dtt, LEAF *lf, VEC *gv, double cn, int *nghlf, int *maxlnf, int iret)
{

    int isec[3]={0, 0, 0};

    VEC *pa=NULL, *pb=NULL;
    VEC cpv, cpv1, cpv2;
    VEC dpp;

    double norm, vdd;
    double aa, bb, cc;
    double tolarc=TOLARC;

    if(cn < 0.){
       printf("***ERROR***, negative smoothing parameter in %s\n", __FILE__);
       exit(1);
    }

    if(lf->ingh) return ilf;

    ++(*nghlf);

    if(!ilf){
       ilf = (LEAF **)malloc_initl(sizeof(LEAF *));
       mem_er((ilf == NULL) ? 0 : 1, sizeof(LEAF *));
       *maxlnf = 1;
    }
    else if(*nghlf > *maxlnf){
       ilf = (LEAF **)realloc_n(ilf, (*nghlf) * sizeof(LEAF *));
       mem_er((ilf == NULL) ? 0 : 1, (*nghlf) * sizeof(LEAF *));
       *maxlnf = *nghlf;
    }

    *(ilf + *nghlf - 1) = lf;
    lf->ingh = 1;

    if(iret) return ilf;


/* check for neighbours */

    dpp.x = dtt->xdt;
    dpp.y = dtt->ydt;
    dpp.z = dtt->zdt;


/* simple check, verticies in support region */

    pa = gv + lf->ivec[0];

    vdd = cn * dotp(&dpp, pa) - 1.0 + tolarc;
    if(vdd >= 0.0){

       isec[0] = 1;
       isec[2] = 1;
       ilf = sqt_leaf_sample(ilf, dtt, lf->lft, gv, cn, nghlf, maxlnf, 0);
       ilf = sqt_leaf_sample(ilf, dtt, lf->rgt, gv, cn, nghlf, maxlnf, 0);

    }

    pa = gv + lf->ivec[1];

    vdd = cn * dotp(&dpp, pa) - 1.0 + tolarc;
    if(vdd >= 0.0){

       isec[0] = 1;
       isec[1] = 1;
       ilf = sqt_leaf_sample(ilf, dtt, lf->lft, gv, cn, nghlf, maxlnf, 0);
       ilf = sqt_leaf_sample(ilf, dtt, lf->udw, gv, cn, nghlf, maxlnf, 0);

    }

    pa = gv + lf->ivec[2];

    vdd = cn * dotp(&dpp, pa) - 1.0 + tolarc;
    if(vdd >= 0.0){

       isec[1] = 1;
       isec[2] = 1;
       ilf = sqt_leaf_sample(ilf, dtt, lf->rgt, gv, cn, nghlf, maxlnf, 0);
       ilf = sqt_leaf_sample(ilf, dtt, lf->udw, gv, cn, nghlf, maxlnf, 0);

    }  

/* Left (Side Id 1) */

    if(!isec[0]){

       pa = gv + lf->ivec[0];
       pb = gv + lf->ivec[1];

       crosp(pa, pb, &cpv);
       crosp(&dpp, &cpv, &cpv1);
       crosp(&cpv, &cpv1, &cpv2);
       norm = sqrt(dotp(&cpv2, &cpv2));
       normv(&cpv2, norm);
       if(cn * dotp(&cpv2, &dpp) - 1.0 + tolarc >= 0.0){

          aa = dotp(pa, pb);
          bb = dotp(pa, &cpv2);
          cc = dotp(pb, &cpv2);

          if(aa <= ((bb < cc) ? bb : cc))
             ilf = sqt_leaf_sample(ilf, dtt, lf->lft, gv, cn, nghlf, maxlnf, 0);
       }


    }


/* Verticle (Side Id 2) */

    if(!isec[1]){

       pa = gv + lf->ivec[1];
       pb = gv + lf->ivec[2];

       crosp(pa, pb, &cpv);
       crosp(&dpp, &cpv, &cpv1);
       crosp(&cpv, &cpv1, &cpv2);
       norm = sqrt(dotp(&cpv2, &cpv2));
       normv(&cpv2, norm);
       if(cn * dotp(&cpv2, &dpp) - 1.0 + tolarc >= 0.0){
          aa = dotp(pa, pb);
          bb = dotp(pa, &cpv2);
          cc = dotp(pb, &cpv2);

          if(aa <= ((bb < cc) ? bb : cc))
             ilf = sqt_leaf_sample(ilf, dtt, lf->udw, gv, cn, nghlf, maxlnf, 0);

       }

    }


/*  Right (Side Id 3) */

    if(!isec[2]){


       pa = gv + lf->ivec[2];
       pb = gv + lf->ivec[0];


       crosp(pa, pb, &cpv);
       crosp(&dpp, &cpv, &cpv1);
       crosp(&cpv, &cpv1, &cpv2);
       norm = sqrt(dotp(&cpv2, &cpv2));
       normv(&cpv2, norm);
       if(cn * dotp(&cpv2, &dpp) - 1.0 + tolarc >= 0.0){
          aa = dotp(pa, pb);
          bb = dotp(pa, &cpv2);
          cc = dotp(pb, &cpv2);

          if(aa <= ((bb < cc) ? bb : cc))
             ilf = sqt_leaf_sample(ilf, dtt, lf->rgt, gv, cn, nghlf, maxlnf, 0);

       }

    }

    return ilf;
}
