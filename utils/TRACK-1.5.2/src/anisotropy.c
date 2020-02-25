#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "grid.h"
#include "statistic.h"
#include "m_values.h"
#include "proj.h"
#include "p_vecs.h"
#include "boundary.h"

#define  DISTOL   0.00001
#define  MAXDIFF  1.0e+10
#define  TOLEIGN  1.0e-8

/* function to compute the local anisotropy around a feature point. */

void assign_pt(struct dpt * , float * , float * );
struct dpt *search_aniso(struct dpt * , int * , float * , float , int , int ,int , int , int , float , float * );
void *shape_vecs(struct dpt * , int, struct cvecs * , int );
double  invtan(double , double );
void shape_setup(struct boundary_cntl * , int );
int boundary_find(struct object * , struct boundary_cntl * , int * , int * , int );
float sys_area(struct boundary_pt * , int , int , int , float , int * , int );
int inside(int , struct boundary_pt * , float , float );
float area_int(struct boundary_pt * , int , float * , float * , int , int , int , int , float * , float , float , int );
void getdxdy(float * , float * , int , int , int );

float bilinear_intrp(float * , float , float , float , int );

int iian=0;
int *xan=NULL, *yan=NULL;
int an_plt='n';

extern GRID *gr, *gt;
extern int x1u, y1u;
extern int pb;

extern int tf;
extern int tom;
extern int dfil;

extern float period;

extern PROJ *pp;

extern float *ap;

void anisotropy(struct frame_objs *fo, int msg)

{

   static int noa=0;
   static int vtype='z';
   static int reg_th=1;
   static int direc=0;
   static int ssmth='n';
   static int aar='n';
   static int iff=1;
   static int atype=0;
   static int bnop=0;

   int i,j,k,l;
   int npt=0;
   int t_dfil=0;

   int pxmx, pxmn, pymx, pymn;
   int ydim, xdim, dim;
   int itt=0;
   int iian_tmp=0;
   int nn=0;
   int iwarn=0;
   int iover = 0;

   int nerx=0, nery=0;
   int ntx, nty;
   
   int intrp=0;

   int *iorder=NULL;

   static float ffrac=1.0, athr=0.0;
   float fpth, apth;
   float *aa=NULL, *ab=NULL;
   float nrx, nry, nrxo, nryo;
   float mrxo, mryo;
   float obx, oby;
   float d1, d2;

   float mindis, dis1, dis2;
   float msk=-2.0;
   float mxstr=0.0;
   float *forder=NULL;
   float vv;

   double snp, csp, snt, cst;

   struct object *ob;
   struct object ob_tmp, ob_area;
   struct point *ptt;
   struct feature_pts *fpts=NULL, *fptt=NULL;
   struct boundary_pt *btt=NULL;

   static struct boundary_cntl bcntl = {0, 0, 0, 'n', 0};

   struct dpt *pta=NULL, pdir, pbnd;
   struct boundary_pt *bpt=NULL, *bbpt=NULL;
   struct cvecs eigp;
   
   iian = 0;

   t_dfil = dfil;

   if(!msg){

/* re-initialize */

       noa=0;
       vtype='z';
       reg_th=1;
       direc=0;
       ssmth='n';
       aar='n';
       iff=1;
       atype=0;
       memset(&bcntl, 0L, sizeof(struct boundary_cntl));
       bcntl.ofill = 'n'; 

      if(tom != 'g'){
 
         printf("***WARNING***, shape and area computation only available on the sphere\r\n"
                "               no anisotropy measure computed.                    \n\n");

         noa = 1;

      }

      if(noa == 1) {
         printf("****WARNING****, no shape or area calculation performed\n\n");
         return;

      }


      printf("Do you want consistent orientation vectors. For none or   \r\n"
             "zonal orientation vectors or meriodinal orientation       \r\n"
             "vectors. Input 'n', 'z' or 'm'.                           \n\n");

      scanf("\n");
      vtype = getchar();

      printf("Input a region threshold for shape determination in point numbers,\r\n"
             "any region with fewer points will be treated as isotropic.        \n\n");

      scanf("%d", &reg_th);

      printf("Do you want shape determined centered on the direction of the     \r\n"
             "feature point or centered on the direction of the mean            \r\n"
             "(principle) vector direction, '0' for feature point, '1' for mean.\n\n");

      scanf("%d", &direc);

      printf("Do you want to allow object overlap, 'y' or 'n'\n\n");
      scanf("\n");
      if(getchar() == 'y') iover = 1;


      if(reg_th < 1) {

         printf("****WARNING****, region thresholding for shape determination must be >=1\r\n"
                "                 re-setting to 1.                                    \n\n");
         reg_th = 1;

      }

      printf("Do you want to smooth the region shape before computing the anisotropy descriptors. 'y' or 'n'.\n\n");
      scanf("\n");
      if((ssmth = getchar()) == 'y') {
         printf("To compute anisotropy from a sub-region , what fraction of a feature points  \r\n"
                "value defines the shape associated with the feature point?                   \n\n");

         scanf("%f", &ffrac);

         shape_setup(&bcntl, 0);

      }

      printf("Will you want to plot the regions from which the shape is computed? 'y' or 'n'\n");
      scanf("\n");
      an_plt = getchar();

      printf("Do you want the region areas computed, 'y' or 'n'\n\n");
      scanf("\n");
      aar = getchar();

      if(aar == 'y'){

         printf("Set the fraction threshold for area calculation\n\n");
         scanf("%f", &athr);

         printf("Do you want just area,                       '0' or, \r\n"
                "Area multiplied by feature point itensity,   '1' or, \r\n"
		"Area integrated field                        '2'     \n\n");
         scanf("%d", &atype);
	 
	 if(atype < 0 || atype > 2){
	    printf("****ERROR****, wrong value for option.\n\n");
	    exit(1);
	 }

         if(gr->prty) {
	    printf("Compute area on the projection rather than the sphere, '0' for no or '1' for yes.\n\n");
	    scanf("%d", &bnop);
	    if(bnop < 0 || bnop > 1){
	       printf("****ERROR****, wrong value for option.\n\n");
	       exit(1);
	    }	    
	 
	 }


      }


   }

   if(noa == 1) return;  

   intrp = 0;

   for(i=0; i < fo->obj_num; i++){

      ob = (fo->objs) + i;

      if(! ob->ext){

        pxmx = pxmn = ob->pt->x;
        pymx = pymn = ob->pt->y;

        for(j=0; j < ob->point_num; j++){

            ptt = (ob->pt) + j;

            if(ptt->x > pxmx) pxmx = ptt->x;
            else if(ptt->x < pxmn) pxmn = ptt->x;

            if(ptt->y > pymx) pymx = ptt->y;
            else if(ptt->y < pymn) pymn = ptt->y;

        }


        ob->ext = (struct extent * )malloc_initl(sizeof(struct extent));
        mem_er((ob->ext == NULL) ? 0 : 1, sizeof(struct extent));

        ob->ext->x1 = pxmn;
        ob->ext->x2 = pxmx;
        ob->ext->y1 = pymn;
        ob->ext->y2 = pymx;

      }

      else {

         pxmn = ob->ext->x1;
         pxmx = ob->ext->x2;
         pymn = ob->ext->y1;
         pymx = ob->ext->y2;

      }


      xdim = pxmx - pxmn + 3;
      ydim = pymx - pymn + 3;

      dim = xdim*ydim;

      aa = (float *)calloc(dim, sizeof(float));
      mem_er((aa == NULL) ? 0 : 1, dim * sizeof(float));

      ab = (float *)calloc(dim, sizeof(float));
      mem_er((ab == NULL) ? 0 : 1, dim * sizeof(float));


      for(j=0; j < dim; j++) *(ab + j) = -1.0;

      for(j=0; j < ob->point_num; j++){

         ptt = (ob->pt) + j;
         *(ab + (ptt->y - pymn + 1) * xdim + ptt->x - pxmn + 1) = ptt->val;

      }

      memcpy(aa, ab, dim * sizeof(float));

/* determine feature point intensity order */

      iorder = (int *)calloc(ob->fet->feature_num, sizeof(int));
      mem_er((iorder == NULL) ? 0 : 1, (ob->fet->feature_num) * sizeof(int));

      forder = (float *)calloc(ob->fet->feature_num, sizeof(float));
      mem_er((forder == NULL) ? 0 : 1, (ob->fet->feature_num) * sizeof(float));

      for(j=0; j < ob->fet->feature_num; j++){
          fpts = (ob->fet->fpt) + j;
          *(forder + j) = fpts->ostr;          
      }

      for(j=0; j < ob->fet->feature_num; j++){

          mxstr = 0.0;
          itt = 0;

          for(k=0; k < ob->fet->feature_num; k++){

              if(*(forder + k) > mxstr) {mxstr = *(forder + k); itt = k;}
          }

          *(iorder + j) = itt;

          *(forder + itt) = 0.0; 


      }

      free(forder);


      for(j=0; j < ob->fet->feature_num; j++){ 

          if(iover) memcpy(aa, ab, dim * sizeof(float));


/* ------- test code ------ */

/*for(k=0; k < ydim; k++){
for(l=0; l < xdim; l++) printf("%+6.4f ", *(aa + k*xdim +l));
printf("\n");
}
printf("\n\n"); */

/* ----------------------- */


          fpts = (ob->fet->fpt) + *(iorder + j); 

          if(fpts->str < CHECK_PT) {
            fpts->r_sh = fpts->ornt[0] = fpts->ornt[1] = fpts->area = 0.;
            continue;
          }

/* check ordering of old points */

          if(tf == 7 || tf == 8){

             nrx = (fpts->x).xy;
             nry = (fpts->y).xy;
             nrxo = (fpts->ox).xy;
             nryo = (fpts->oy).xy;

             for(k=0; k < ob->fet->feature_num; k++){

                 fptt = ob->fet->fpt + k;

                 if(k != *(iorder + j) && fptt->str < CHECK_PT){ 
                    mrxo = (fptt->ox).xy;
                    mryo = (fptt->oy).xy; 
                    d1 = nrx - nrxo;
                    d2 = nry - nryo;
                    dis1 = d1*d1 + d2*d2;
                    d1 = nrx - mrxo;
                    d2 = nry - mryo;
                    dis2 = d1*d1 + d2*d2;
  
                    if(dis2 < dis1){

                      nrxo = mrxo;
                      nryo = mryo;
                      vv = fpts->ostr;
                      fpts->ostr = fptt->ostr;
                      (fpts->ox).xy = mrxo;
                      (fpts->oy).xy = mryo;
                      (fptt->ox).xy = (fpts->ox).xy;
                      (fptt->oy).xy = (fpts->oy).xy;
                      fptt->ostr = vv;
                      
                    }
                 
                 }
             }

             nrx = nrxo;
             nry = nryo;
             (fpts->ox).xy = nrx;
             (fpts->oy).xy = nry;

          }

          else {

             if(tf == 3){
               nrx = obj_xreal((fpts->ox).ixy + x1u - 2);
               nry = *(gr->ygrid + (fpts->oy).ixy + y1u -2);
             }

             else{
               nrx = (fpts->ox).xy;
               nry = (fpts->oy).xy;
             }

          }

          fpth = fpts->str * ffrac;
          apth = fpts->str * athr;

          if(fpth < 0.) continue;

          npt = 0;


/* search for nearest maximum grid point */

          mindis = -1.0;

          for(k=0; k < ob->point_num; k++){

               ptt = (ob->pt) + k;

               obx = obj_xreal(ptt->x + x1u - 2);

               oby = *(gr->ygrid + ptt->y + y1u - 2);


               d1 = nrx - obx;
               d2 = nry - oby;
  

               dis1 = d1*d1 + d2*d2;


               if(mindis < 0.0) {

                 mindis = dis1; 
                 nerx = ptt->x - pxmn + 1;
                 nery = ptt->y - pymn + 1;

               }

               else if(dis1 < mindis){

                   mindis = dis1;
                   nerx = ptt->x - pxmn + 1;
                   nery = ptt->y - pymn + 1;
               }

          }

          vv = *(ab + nery * xdim + nerx);
          ntx = nerx;
          nty = nery;

          for(k=-1; k <= 1; k++){
             for(l=-1; l <= 1; l++){
                 if(vv < *(ab + (nery + k) * xdim + nerx + l)){
                    vv = *(ab + (nery + k) * xdim + nerx + l);
                    ntx = nerx + l;
                    nty = nery + k;
                 }
             }
          }

          nerx = ntx;
          nery = nty;

          msk = -2.0;
          iian_tmp = iian;

          pta = search_aniso(pta, &npt, aa, apth, nerx, nery, xdim, pxmn, pymn, vv, &msk);


/* ------- test code ------ */

/*for(k=0; k < ydim; k++){
for(l=0; l < xdim; l++){
    if(*(aa + k*xdim +l)>0)printf("%3d ", 1);
    else if(*(aa + k*xdim +l) < -1) printf("%3d ", (int)(*(aa + k*xdim +l)));
    else printf("%3d ", (int)(*(aa + k*xdim +l)));
}
printf("\n");
}
printf("\n\n"); */

/* ----------------------- */


/* compute area of a sub-region */


          fpts->area = 0.0;

          if(aar == 'y' && npt > 1 ){

             memset(&ob_area, 0L, sizeof(struct object));

             ob_area.pt = (struct point *)calloc(npt, sizeof(struct point));
             mem_er((ob_area.pt == NULL) ? 0 : 1, npt * sizeof(struct point));

             itt = 0;

             for(k=0; k < ydim; k++){

                for(l=0; l < xdim; l++){

                   if(*(aa + k*xdim + l) < -1) {
                     (ob_area.pt + itt)->x = l + pxmn - 1;
                     (ob_area.pt + itt)->y = k + pymn - 1;
                     (ob_area.pt + itt)->val = 1.0;

                     ++itt;
 
                   }

                }

             }

             ob_area.point_num = itt;

             boundary_find(&ob_area, NULL, NULL, NULL, 1);


             for(k=0; k < ob_area.bound_num; k++){

                btt = ob_area.bound + k;

                (btt->x).xy = obj_xreal((btt->x).ixy + x1u - 2);

                (btt->y).xy = *(gr->ygrid + (btt->y).ixy + y1u - 2);   

             }


             if(!bnop){

                bbpt = ob_area.bound;
                nn = ob_area.bound_num;

                for(k=0; k < ob_area.bound_num; k++){
                    bpt = ob_area.bound + k;
                    pbnd.ic = 0;                
                    assign_pt(&pbnd, &(bpt->x).xy, &(bpt->y).xy);
                    if(!pbnd.ic){
                      (bbpt->x).xy = pbnd.lng;
                      (bbpt->y).xy = pbnd.lat;

                      ++bbpt;
                    }
                    else --nn;
                 
                }

                ob_area.bound_num = nn;
	     
	     }

             fpts->area = sys_area(ob_area.bound, ob_area.bound_num, iff, atype, fpts->str, &intrp, bnop);

             iff = 0;

             free(ob_area.pt);
             free(ob_area.bound);


          }


/* Get sub-region and smooth boundary */

          if(ssmth == 'y' && npt >= reg_th){

             memset(&ob_tmp, 0L, sizeof(struct object));

             ob_tmp.pt = (struct point *)calloc(npt, sizeof(struct point));
             mem_er((ob_tmp.pt == NULL) ? 0 : 1, npt * sizeof(struct point));

             itt = 0;

             for(k=0; k < ydim; k++){

                 for(l=0; l < xdim; l++){

                    if(*(aa + k*xdim + l) < -1 && *(ab + k*xdim + l) >= fpth) {
                       (ob_tmp.pt + itt)->x = l + pxmn - 1;
                       (ob_tmp.pt + itt)->y = k + pymn - 1;
                       (ob_tmp.pt + itt)->val = 1.0;
 
                       ++itt;

                    }

                 }

             }
 
             ob_tmp.point_num = itt;

             boundary_find(&ob_tmp, &bcntl, NULL, NULL, 0);

             free(pta);


             npt = ob_tmp.point_num;

             pta = (struct dpt *)calloc(npt, sizeof(struct dpt));
             mem_er((pta == NULL) ? 0 : 1, npt * sizeof(struct dpt));

             iian = iian_tmp;

             if(an_plt == 'y'){

                xan = (int *)realloc_n(xan, (iian + npt) * sizeof(int));
                mem_er((xan == NULL) ? 0 : 1, (iian + npt) * sizeof(int));

                yan = (int *)realloc_n(yan, (iian + npt) * sizeof(int));
                mem_er((yan == NULL) ? 0 : 1, (iian + npt) * sizeof(int));


                for(k=0; k < npt; k++){

                   ptt = ob_tmp.pt + k;

                   *(xan + iian_tmp + k) = ptt->x + x1u - 1;

                   if(*(xan + iian_tmp + k) < 1) *(xan + iian_tmp + k) += gr->ix;
                   else if(*(xan + iian_tmp + k) > gr->ix) *(xan + iian_tmp + k) -= gr->ix;

                   *(yan + iian_tmp + k) = ptt->y + y1u - 1;
                   ++iian;

                }

             }

             for(k=0; k < npt; k++){

                ptt = ob_tmp.pt + k;

                obx = obj_xreal(ptt->x + x1u - 2);

                oby = *(gr->ygrid + ptt->y + y1u - 2);

                assign_pt(pta+k, &obx, &oby);

             }

             free(ob_tmp.pt);
             free(ob_tmp.bound);


          }

          if(npt >= reg_th) {

            assign_pt(&pdir, &nrx, &nry); 

            eigp.p1[0] = pdir.xdt;
            eigp.p1[1] = pdir.ydt;
            eigp.p1[2] = pdir.zdt;


            shape_vecs(pta, npt, &eigp, direc);

            if(eigp.p1[3] < 0. || eigp.p2[3] < 0. || eigp.p3[3] < 0.){

               if(!iwarn)
                  printf("****WARNING****, negative eigenvalues for shape determination.\r\n"
                         "                 e1 = %f, e2 = %f, e3 = %f\n\n", eigp.p1[3], eigp.p2[3], eigp.p3[3]);

               iwarn = 1;

               if(eigp.p1[3] < 0.){
                  if(fabs(eigp.p1[3]) < TOLEIGN) eigp.p1[3] = 0.0;
               }
               if(eigp.p2[3] < 0.){
                  if(fabs(eigp.p2[3]) < TOLEIGN) eigp.p2[3] = 0.0;
               }
               if(eigp.p3[3] < 0.){
                  if(fabs(eigp.p3[3]) < TOLEIGN) eigp.p3[3] = 0.0;
               }

            }


            if(direc){

               fpts->r_sh = 1.0 - (eigp.p1[3]/eigp.p2[3]);

               if(pdir.xdt * eigp.p3[0] + pdir.ydt * eigp.p3[1] + pdir.zdt * eigp.p3[2] < 0.)
                 {eigp.p3[0] *= -1.; eigp.p3[1] *= -1.; eigp.p3[2] *= -1.;}

               nrx = invtan(eigp.p3[1], eigp.p3[0]);

               nry = acos(eigp.p3[2]);

               sincos(nrx, &snp, &csp);
               sincos(nry, &snt, &cst);

               fpts->ornt[0] = -eigp.p2[0] * snp + eigp.p2[1] * csp;
               fpts->ornt[1] = -eigp.p2[0] * cst * csp - eigp.p2[1] * cst * snp + eigp.p2[2] * snt;



            }


            else{

               fpts->r_sh = 1.0 - (eigp.p2[3]/eigp.p3[3]);


/* compute orientation vector */

               sincos(nrx, &snp, &csp);
               sincos(nry, &snt, &cst);

               fpts->ornt[0] = -eigp.p3[0] * snp + eigp.p3[1] * csp;
               fpts->ornt[1] = -eigp.p3[0] * cst * csp - eigp.p3[1] * cst * snp + eigp.p3[2] * snt;

            }


/*            fpts->ornt[0] *= fpts->r_sh;
            fpts->ornt[1] *= fpts->r_sh;    */


            if(vtype != 'n'){

               if(((fpts->ornt[0] <= 0. && fpts->ornt[1] < 0.)  || 
                  (fpts->ornt[0] < 0. && fpts->ornt[1] >= 0. )) &&
                   vtype == 'z'                                    ){

                 fpts->ornt[0] *= -1.;  
                 fpts->ornt[1] *= -1.;

               } 

               else if(((fpts->ornt[0] < 0. && fpts->ornt[1] >= 0.)  || 
                       (fpts->ornt[0] >= 0. && fpts->ornt[1] > 0. )) &&
                        vtype == 'm'                                    ){

                 fpts->ornt[0] *= -1.;  
                 fpts->ornt[1] *= -1.;

               }

            }
              

          }

          else {

             fpts->r_sh = fpts->area = 0.0;
             fpts->ornt[0] = fpts->ornt[1] = 0.0;

          }

/* ------- test code ------ */


/* printf("%f %f %f\n", fpts->r_sh, fpts->ornt[0], fpts->ornt[1]); */

/* ----------------------- */
 
          free(pta);

          if(!iover){
             for(k=0; k < dim; k++) {
                 if(*(aa + k) < -1) *(aa + k) = -1.0;
             }
          }

      }

      free(iorder);
      free(aa);
      free(ab);


   }

   dfil = t_dfil;


   return;

}


void assign_pt(struct dpt *pta, float *xx, float *yy)

{
   int ic=1;

   double sn1, sn2, cn1, cn2;


/* push point back into domain */

    if(tf != 7 && pb =='y'){
       if(*xx < *(gr->xgrid)) *xx += period;
       else if(*xx > *(gr->xgrid + gr->ix - 1)) *xx -= period;
    }

/* project back to the sphere */

    if(gr->prty) {

       switch(gr->prgr){
           case 0:
             (pp->prj2)(xx, 'x', 0);
             (pp->prj2)(yy, 'y', 0);
             break;
           case 1:
             ic=azimuthal(yy, xx, gr->alat, gr->alng, 0, gr->prty);
             break;
       }

    }

    if(!ic) {pta->ic = 1; return;}

    pta->lng = *xx;
    pta->lat = *yy;

/* convert to cartesian co-ordinates (x, y, z) */

    *xx *= FP_PI;
    *yy = FP_PI2 - *yy * FP_PI;

    if(*yy < 0.) *yy = 0.;

    sincos(*xx, &sn1, &cn1);
    sincos(*yy, &sn2, &cn2);

    pta->xdt = sn2 * cn1;
    pta->ydt = sn2 * sn1;
    pta->zdt = cn2;

    return;

}

/* search_aniso, does a 4-connected search for neighbours with value 
   within the threshold.                                              */

struct dpt *search_aniso(struct dpt *pta, int *npt, float *aa, float fpth, int nerx, int nery, int xdim, int pxmn, int pymn, float val, float *msk)

{

    float obxx=0., obyy=0.;

    float app;

    struct dpt *ptaa=NULL;

    app = *(aa + nery * xdim + nerx);

    if(app >= fpth && app <= val){

       ++(*npt);

      *(aa + nery * xdim + nerx) = *msk;

       *msk -= 1.0;

       obxx = obj_xreal(nerx + pxmn + x1u - 3);

       obyy = *(gr->ygrid + nery + pymn + y1u - 3);

       if(an_plt == 'y'){

          ++iian;

          if(!xan){

             xan = (int *)malloc_initl(sizeof(int));
             mem_er((xan == NULL) ? 0 : 1, sizeof(int));

             yan = (int *)malloc_initl(sizeof(int));
             mem_er((yan == NULL) ? 0 : 1, sizeof(int));

           }

           else {

             xan = (int *)realloc_n(xan, iian * sizeof(int));
             mem_er((xan == NULL) ? 0 : 1, iian * sizeof(int));

             yan = (int *)realloc_n(yan, iian * sizeof(int));
             mem_er((yan == NULL) ? 0 : 1, iian * sizeof(int));

           }

           *(xan + iian - 1) = nerx + pxmn + x1u - 2;

           if(*(xan + iian - 1) < 1) *(xan + iian - 1) += gr->ix;
           else if(*(xan + iian - 1) > gr->ix) *(xan + iian - 1) -= gr->ix;

           *(yan + iian - 1) = nery + pymn + y1u - 2;


       }


       if(!pta){

          pta = (struct dpt *)malloc_initl(sizeof(struct dpt));
          mem_er((pta == NULL) ? 0 : 1, sizeof(struct dpt));

       }

       else {

          pta = (struct dpt *)realloc_n(pta, *npt * sizeof(struct dpt));
          mem_er((pta == NULL) ? 0 : 1, *npt * sizeof(struct dpt));

       }

       ptaa = pta + *npt - 1;
       ptaa->ic = 0;

       assign_pt(ptaa, &obxx, &obyy);

/* eight connected search */

       pta = search_aniso(pta, npt, aa, fpth, nerx+1, nery, xdim, pxmn, pymn, app, msk);
       pta = search_aniso(pta, npt, aa, fpth, nerx, nery+1, xdim, pxmn, pymn, app, msk);
       pta = search_aniso(pta, npt, aa, fpth, nerx-1, nery, xdim, pxmn, pymn, app, msk);
       pta = search_aniso(pta, npt, aa, fpth, nerx, nery-1, xdim, pxmn, pymn, app, msk);

       pta = search_aniso(pta, npt, aa, fpth, nerx+1, nery+1, xdim, pxmn, pymn, app, msk);
       pta = search_aniso(pta, npt, aa, fpth, nerx-1, nery+1, xdim, pxmn, pymn, app, msk);
       pta = search_aniso(pta, npt, aa, fpth, nerx+1, nery-1, xdim, pxmn, pymn, app, msk);
       pta = search_aniso(pta, npt, aa, fpth, nerx-1, nery-1, xdim, pxmn, pymn, app, msk);

    }


    return pta;

}

float sys_area(struct boundary_pt *bpt, int nbp, int iff, int atype, float fptstr, int *intrp, int bnop)
{

   int i=0, j=0, iper=0;
   int ibf=0, nt=0, ibb=0;
   static int iref=0;
   int dim=0;
   static int ifirst=0;
   static int nx=0, ny=0;

   float x1, x2, y1, y2;
   float diff=0.0;
   static float dx=0.0, dy=0.0;
   static float *xg=NULL, *yg=NULL;  /* refined grid positions */
   static float *fld=NULL;           /* refined grid values    */
   float area=0.0;
   float period2=0.5*period;

   struct boundary_pt *btmp1=NULL, *btmp2=NULL, *bbt=NULL, *btmp=NULL;
   
   GRID *gtemp=NULL;

   if(iff > 0) ifirst = 0; /* reset initialization flag    */
   else if (iff < 0 && iref) {     /* free memory for refined grid */
      free(xg);
      free(yg);
      if(atype == 2) free(fld);
      ifirst = 0;
      return 0.0;
   }
   else if(iff < 0 && !iref) return 0.0;


   if(!ifirst){

      printf("****WARNING****, the calculation of the sub-region areas   \r\n"
             "                 expects the sub-regions to be much smaller\r\n"
             "                 than the domain.                          \n\n");  

      printf("****INFORMATION****, current region is X -- (%f, %f) \r\n"
             "                                       Y -- (%f, %f) \n\n",
             *(gt->xgrid), *(gt->xgrid + gt->ix - 1), *(gt->ygrid), 
             *(gt->ygrid + gt->iy - 1));
	     
      printf("Use current grid or use grid refinment, '0' for current or '1' for refined\n\n");
      scanf("%d", &iref);
      
      if(iref < 0 || iref > 1){
         printf("****ERROR****, wrong input for option.\n\n");
	 exit(1);
      }
      
      if(!iref){
         xg = gt->xgrid;
	 yg = gt->ygrid;
	 nx = gt->ix;
	 ny = gt->iy;
	 
	 if(!(gt->ixfrm)) dx = (*(xg + nx - 1) - *xg) / (nx - 1);
	 if(!(gt->iyfrm)) dy = (*(yg + ny - 1) - *yg) / (ny - 1);
	 
	 if(atype == 2) fld = ap;
      }
      else{

         printf("Specify the grid refinement required for computing area\r\n"
                "in terms of the region and the  number of grid nodes   \r\n"
                "in X and Y.                                            \r\n"
                "Input full region, x1, x2, y1, y2.                     \n\n");

         scanf("%f %f %f %f", &x1, &x2, &y1, &y2);
	 
	 if(x1 < *(gt->xgrid)) x1 = *(gt->xgrid);
	 if(x2 > *(gt->xgrid + gt->ix - 1)) x2 = *(gt->xgrid + gt->ix - 1);
	 if(y1 < *(gt->ygrid)) x1 = *(gt->ygrid);
	 if(y2 > *(gt->ygrid + gt->iy - 1)) y2 = *(gt->ygrid + gt->iy - 1);	 

         printf("Input the number of points required in the refined grid\r\n"
                "in X and Y.                                            \n\n");
      
         scanf("%d %d", &nx, &ny);

         dx = (x2 - x1) / (nx - 1);
         dy = (y2 - y1) / (ny - 1);


         xg = (float * )calloc(nx, sizeof(float));
         mem_er((xg == NULL) ? 0 : 1, nx * sizeof(float));
         yg = (float * )calloc(ny, sizeof(float));
         mem_er((yg == NULL) ? 0 : 1, ny * sizeof(float));


         for(i=0; i < nx; i++) *(xg + i) = x1 + dx * (float) i;
         for(i=0; i < ny; i++) *(yg + i) = y1 + dy * (float) i;
	 
	 
	 if(atype == 2) {
	    dim = nx * ny;
	    fld = (float *)calloc(dim, sizeof(float));
	    mem_er((fld == NULL) ? 0 : 1, dim * sizeof(float));
	 }

/*
         if(*xg < x1) *xg = x1;
         if(*(xg + nx - 1) > x2) *(xg + nx - 1) = x2;
         if(*yg < y1) *yg = y1;
         if(*(yg + ny - 1) > y2) *(yg + ny - 1) = y2;  
*/

      }  

      ifirst = 1;
      
      dx *= FP_PI;
      dy *= FP_PI;

   }

   if(!bpt){
      printf("****WARNING****, no boundary information for determining area\n\n");
      return 0.0;

   }
   
   
   if(iref && atype == 2 && !(*intrp)){
   
      gtemp = gr;
      gr = gt;
      
      for(i=0; i < ny; i++){
          for(j=0; j < nx; j++){
              *(fld + i * nx + j) =  bilinear_intrp(ap, *(xg + j), *(yg + i), MVAL, 2);
	  }      
      }
   
      gr = gtemp;       
      gtemp = NULL;
   
      *intrp = 1;
   
   }


/* Find extent of object boundary */

/* Does object straddel a periodic boundary */

   for(i=0; i < nbp-1; i++){

       btmp1 = bpt + i;
       btmp2 = btmp1 + 1;

       diff = fabs((btmp2->x).xy - (btmp1->x).xy);
       if(diff > period2) {iper = 1; break;}

   }

   if(!iper) area = area_int(bpt, nbp, xg, yg, nx, ny, atype, iref, fld, dx, dy, bnop);

   else {

       btmp = (struct boundary_pt * )calloc(nbp+4, sizeof(struct boundary_pt));
       mem_er((btmp == NULL) ? 0 : 1, (nbp+4) * sizeof(struct boundary_pt));

       ibf = 0;

/* find first boundary cross-over point */

       for(i=0; i < nbp-1; i++){

          btmp1 = bpt + i;
          btmp2 = btmp1 + 1;

          diff = fabs((btmp2->x).xy - (btmp1->x).xy);
	  
          if(diff > period2) {ibf = i+1; break;}

       }

       ibb = ibf;

       while(ibf >= 0){

          if(ibf == nbp - 1) break;
	  
          btmp1 = bpt + ibf;
          btmp2 = btmp1 + 1;


          if(fabs((btmp2->x).xy - (btmp1->x).xy) > period2) { ++ibf; continue;}

          nt = 1;

          if((btmp1->x).xy < period2) {

             (btmp->x).xy = 0.0;
             (btmp->y).xy = (btmp1->y).xy;

          }

          else {

             (btmp->x).xy = period;
             (btmp->y).xy = (btmp1->y).xy;
          }


          while(fabs((btmp2->x).xy - (btmp1->x).xy) < period2){
             ++nt;
             bbt = btmp + nt - 1;
             (bbt->x).xy = (btmp1->x).xy;
             (bbt->y).xy = (btmp1->y).xy;

             ++ibf;
             if(ibf == nbp-1) ibf = 0;
             btmp1 = bpt + ibf;
             btmp2 = btmp1 + 1;

          }

          ++ibf;
          nt += 3;
          ++bbt;

          (bbt->x).xy = (btmp1->x).xy;
          (bbt->y).xy = (btmp1->y).xy;

          ++bbt;
          if((btmp1->x).xy < period2){
            (bbt->x).xy = 0.0;
            (bbt->y).xy = (btmp1->y).xy;
          }
          else {
            (bbt->x).xy = period;
            (bbt->y).xy = (btmp1->y).xy;

          }

          ++bbt;
          (bbt->x).xy = (btmp->x).xy;
          (bbt->y).xy = (btmp->y).xy;

          area += area_int(btmp, nt, xg, yg, nx, ny, atype, iref, fld, dx, dy, bnop);

          if(ibf == ibb) break;

       }

       free(btmp);

   }

/* printf("area %e\n", area); */

   if(atype == 1) area = fptstr * area;
   

   return area;

}


float area_int(struct boundary_pt *bpt, int nbp, float *xg, float *yg, int nx, int ny, int atype, int iref, float *fld, float dx, float dy, int bnop)
{

    int i,j;
    int i1=0, i2=0, j1=0, j2=0;

    float mxlng, mnlng;
    float mxlat, mnlat;
    float area = 0.0;
    float d1, d2, diff;
    float coslat=0.0;

    struct boundary_pt *btmp=NULL;

    mxlng = mnlng = (bpt->x).xy;
    mxlat = mnlat = (bpt->y).xy;

    for(i=0; i < nbp; i++){
 
        btmp = bpt + i;

        if((btmp->x).xy > mxlng) mxlng = (btmp->x).xy;
        if((btmp->x).xy < mnlng) mnlng = (btmp->x).xy;

        if((btmp->y).xy > mxlat) mxlat = (btmp->y).xy;
        if((btmp->y).xy < mnlat) mnlat = (btmp->y).xy;

    }

/* nearest grid points */

    d1 = d2 = MAXDIFF;

    for(i=0; i < nx; i++){

        diff = fabs(mnlng - *(xg + i));
        if(diff < d1) {d1 = diff; j1 = i;}
        diff = fabs(mxlng - *(xg + i));
        if(diff < d2) {d2 = diff; j2 = i;}

    }
    
    if(fabs((*xg + period) - *(xg + j2)) < TOLGRID) --j2;

    d1 = d2 = MAXDIFF;

    for(i=0; i < ny; i++){

        diff = fabs(mnlat - *(yg + i));
        if(diff < d1) {d1 = diff; i1 = i;}
        diff = fabs(mxlat - *(yg + i));
        if(diff < d2) {d2 = diff; i2 = i;}

    }


    if(atype == 2){
       for(i=i1; i <= i2; i++){
           if(!bnop) coslat = cos(FP_PI * *(yg + i));
           for(j=j1; j <= j2; j++){
               if(inside(nbp, bpt, *(xg + j), *(yg + i))){
	          if(!iref) {
		     if(gt->ixfrm) getdxdy(xg, &dx, nx, j, 0);
		     if(gt->iyfrm) getdxdy(yg, &dy, ny, i, 1);
		  }	       
                  if(!bnop) {area += dx * dy * coslat * *(fld + i * nx + j);}
		  else {area += dx * dy;}
               }
  
           }

       }    
    }
    else{
    
       for(i=i1; i <= i2; i++){
           if(!bnop) coslat = cos(FP_PI * *(yg + i));
           for(j=j1; j <= j2; j++){
               if(inside(nbp, bpt, *(xg + j), *(yg + i))){
	          if(!iref) {
		     if(gt->ixfrm) getdxdy(xg, &dx, nx, j, 0);
		     if(gt->iyfrm) getdxdy(yg, &dy, ny, i, 1);		  
		  }
                  if(!bnop) {area += dx * dy * coslat;}
		  else {area += dx * dy;}
               }
  
           }

       }
    
    }


    return area;

}

void getdxdy(float *xy, float *dxy, int nxy, int ixy, int idr)
{

   int i1, i2;
   
   float dxy1=0.0, dxy2=0.0;
   
   i1 = ixy - 1;
   i2 = ixy + 1;
   
   if(!idr){   /* x-direction */

     if(i1 >= 0) dxy1 = *(xy + ixy) - *(xy + i1);
     else if(i1 < 0 && pb == 'y') dxy1 = *(xy + ixy) - (*(xy + nxy - 2) - period);
     else dxy1 = *(xy + i2) - *(xy + ixy);

     if(i2 < nxy) dxy2 = *(xy + i2) - *(xy + ixy);
     else if(i2 >= nxy && pb == 'y') dxy2 =  *(xy + nxy - 1) - *(xy + ixy) + *(xy + 1);
     else dxy2 = *(xy + ixy) - *(xy + i1);
     
   }
   else if(idr){  /* y-direction */

     if(i1 >= 0) dxy1 = *(xy + ixy) - *(xy + i1);
     else dxy1 = *(xy + i2) - *(xy + ixy);

     if(i2 < nxy) dxy2 = *(xy + i2) - *(xy + ixy);
     else dxy2 = *(xy + ixy) - *(xy + i1);
   
   }
   
   else{
     printf("****ERROR****, incorrect direction specifier.\n\n");
     exit(1);
   }

   *dxy = 0.5 * (dxy1 + dxy2) * FP_PI;

}

