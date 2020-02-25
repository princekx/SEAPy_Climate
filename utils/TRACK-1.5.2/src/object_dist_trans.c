#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_im.h"
#include "grid.h"
#include "boundary.h"



/* function to compute feature points based on the distance transform.
   Object shapes can be smoothed using the distance transform.
   Feature points are determined from the maxima in the transform?     

   See Rosenfeld and Kak, Digital Picture Processing.                   */

void inv_dist_trans(int * , int , int , int );
void dist_trans(int * , int , int , int );
int *segment_maxima(int * , int , int );
void hierarc_segment(struct image ** , struct image * , int , struct frame_objs * , int , int , int );
int powi(int , int );
float obj_xreal(int );
void shape_setup(struct boundary_cntl * , int );
void frame_boundaries(struct frame_objs * , struct boundary_cntl * , int * , int * );

extern int cc;
extern int x1u, y1u;
extern int dfil;
extern int pb;
extern int tf;

extern GRID *gr;


void object_dist_trans(struct frame_objs *ff, struct frame_objs *fp, int *llb, int *rrb, int fst, int nfpt)

{

    static int osm='n';
    static int mth=0;
    static int ex=0;
    static int rdob='n';
    static int oset=1;
    static int ifs='n';

    static struct boundary_cntl bcntl = {0, 0, 0, 'n', 0};

    int i, j, k;
    int pxmx, pxmn, pymx, pymn, ydim, xdim;
    int dim;
    int nb = 1;
    int nobpt=0;
    int level, lmax, p2;
    int offx, offy;
    int oobptn=0;
    int ix, iy;
    int nly, nry;

    int *aa=NULL, ab;

    float *ob_field=NULL;

    float sx, sy, sstr;


    struct object *ob, *obb;
    struct point *ptt;
    struct feature_pts *fpts;
    struct image *ia=NULL, **hierarchy=NULL;


    if(fst){

      if(nfpt) {

         printf("What object shape smoothing is required, either distance transform smoothing, or \r\n"
                "Fourier smoothing with object fill.                                               \n\n");

         printf("Do you want to smooth the object shapes with distance transforms, 'y' or 'n'\n\n");

         scanf("\n");
         osm = getchar();

         if(osm == 'n'){

            printf("Do you want to smooth the object shapes with fourier descriptors, 'y' or 'n'. \n\n");
            scanf("\n");
            ifs = getchar();

            if(ifs == 'y') shape_setup(&bcntl, 0);

         }

      }

      else {

         printf("What type of object smoothing is required, Distance Transform (d) or Fourier (f).\n\n");
         scanf("\n");
         if(getchar() == 'f'){
            osm = 'n';
            ifs = 'y';
            shape_setup(&bcntl, 0);
         }

         else osm = 'y';

      }

      if(osm == 'y'){

         printf("What is the expansion factor >= 0?\n\n");
         scanf("%d", &ex);

         if(ex < 0) ex = 0;

         if(ex){

            printf("What is the maximum suppresion threshold?\n\n");
            scanf("%d", &mth);

            if(mth != ex) printf("***WARNING***, using this value of maximum   \r\n"
                                "               suppresion may lead to a      \r\n"
                                "               gross change in object shape. \n\n");

         } 

         oset += ex;

      }

      printf("do you want the object redefined by the distance transform\r\n"
             "'y' or 'n'\n\n");

      scanf("\n");
      rdob = getchar();


    }

    if(ifs == 'y')
      frame_boundaries(ff, &bcntl, llb, rrb);

    if(rdob == 'y') {

      dfil = 0;

/* re-initialize periodic boundarys if in use */

      if(pb == 'y'){

         for(j=0; j < gr->iy; j++) *(llb + j) = *(rrb + j) = 0;

      }

    }



    for(i=0; i < ff->obj_num; i++){


        ob = (ff->objs) + i;

        pxmx = pxmn = ob->pt->x;
        pymx = pymn = ob->pt->y;

        for(j=0; j < ob->point_num; j++){

            ptt = (ob->pt) + j;

            if(ptt->x > pxmx) pxmx = ptt->x;
            else if(ptt->x < pxmn) pxmn = ptt->x;

            if(ptt->y > pymx) pymx = ptt->y;
            else if(ptt->y < pymn) pymn = ptt->y;

        }

        xdim = pxmx - pxmn + 3 + 2*ex;
        ydim = pymx - pymn + 3 + 2*ex;

        dim = xdim*ydim;

        aa = (int *)calloc(dim, sizeof(int));
        mem_er((aa == NULL) ? 0 : 1, dim * sizeof(int));

        ob_field = (float *)calloc(dim, sizeof(float));
        mem_er((ob_field == NULL) ? 0 : 1, dim * sizeof(float));


/* initialize */


        if(osm == 'y' && ex){
           for(j=0; j < dim; j++) *(aa + j) = 255;
           nb = 0;
        }
        else{
           for(j=0; j < dim; j++) *(aa + j) = 0;
           nb = 1;
        }        

        for(j=0; j < ob->point_num; j++){

           ptt = (ob->pt) + j;

           *(aa+(ptt->y - pymn+oset)*xdim + (ptt->x - pxmn +oset)) = nb;
           *(ob_field+(ptt->y - pymn+oset)*xdim + (ptt->x - pxmn +oset)) = ptt->val;

        }


/* object smoothing using distance transforms */

        if(osm == 'y' && ex) {
           inv_dist_trans(aa, xdim, ydim, mth);

           dist_trans(aa, xdim, ydim, mth); 
        }
        else dist_trans(aa, xdim, ydim, 0);


/* redefine object */

        if(rdob == 'y'){

           ob->b_or_i = '\0';

           nly = nry = 0;

           nobpt = 0;

           for(j=0; j<dim; j++){if(*(aa + j)) ++nobpt; }

/* check for domain boundary and re-allocate memory */

           oobptn = ob->point_num;

           ob->point_num = nobpt;

           ob->mem_pt_size = nobpt;

           if(pb == 'y'){

              for(j=oset; j<ydim-oset; j++){

                 iy = j + pymn - oset;

                 for(k=oset; k <xdim-oset; k++){

                     ix = k + pxmn - oset;

                     ab = *(aa + j*xdim + k);
  
                     if(ab && ix == 1){
                       if(ob->b_or_i == '\0') {ob->b_or_i = 'l';} 
                       else if(ob->b_or_i == 'r') {ob->b_or_i = 'b';}
                       *(llb + iy - 1) = ob->lab;
                       ++nly;
                     }

                     else if(ab && ix == gr->ix){
                       if(ob->b_or_i == '\0') {ob->b_or_i = 'r';} 
                       else if(ob->b_or_i == 'l') {ob->b_or_i = 'b';}
                       *(rrb + iy - 1) = ob->lab;
                       ++nry;
                     }
                  
                 }

              }

              if(ob->b_or_i == 'l') {ob->mem_pt_size += nly;}
              else if (ob->b_or_i == 'r'){ob->mem_pt_size += nry;}
              else if (ob->b_or_i == 'b'){ob->mem_pt_size = (nly > nry) ? nly : nry;}

           }            



           if(ob->mem_pt_size > oobptn){

              ob->pt = (struct point * )realloc_n(ob->pt, ob->mem_pt_size*sizeof(struct point));
              mem_er((ob->pt == NULL) ? 0 : 1, ob->mem_pt_size*sizeof(struct point));

           }

           dfil += ob->mem_pt_size;

           ptt = ob->pt;

           for(j=oset; j<ydim-oset; j++){

              for(k=oset; k <xdim-oset; k++){

                  ab = *(aa + j*xdim + k);

                  if(ab){

                     ptt->x = k + pxmn - oset;
                     ptt->y = j + pymn - oset;
                     ptt->val = (float) ab;

                     ++ptt; 

                  }


              }

           }

           if(cc == 'v') {ptt=ob->pt; for(j=0; j <ob->point_num; j++) (ptt++)->val *= 0.5;}

        }

        if(nfpt){


           aa = segment_maxima(aa, xdim, ydim);


/* connected component search of local maxima */

           level = 1;
           lmax = (xdim >= ydim) ? xdim : ydim;
           p2 = powi(2, level);
 
           while(p2 <= lmax){
                ++level;
                p2 = powi(2, level);
           }
 
           if(p2 != lmax) ++level;

           hierarchy = (struct image ** )calloc(level, sizeof(struct image *));
           mem_er((hierarchy == NULL) ? 0 : 1, level * sizeof(struct image *));

           ia = (struct image *)calloc(p2*p2, sizeof(struct image));
           mem_er((ia == NULL) ? 0 : 1, p2*p2 * sizeof(struct image));
 
           for(j=0; j < ydim; j++){
               for(k=0; k < xdim; k++) {
                   ab = *(aa + j*xdim + k);
                   (ia + j*p2 + k)->pval = (ab) ? 1 : 0;
               }
           }


           hierarc_segment(hierarchy, ia, level, fp, p2, p2, 'v');

           free(ia);

           ff->tot_f_f_num += fp->obj_num;

           ob->fet = (struct features *)malloc_initl(sizeof(struct features));
           mem_er((ob->fet == NULL) ? 0 : 1, sizeof(struct features));

           ob->fet->feature_num = fp->obj_num;

           ob->fet->fpt = (struct feature_pts *)calloc(fp->obj_num, sizeof(struct feature_pts));
           mem_er((ob->fet->fpt == NULL) ? 0 : 1, fp->obj_num * sizeof(struct feature_pts));

           offx = pxmn + x1u - oset - 3;
           offy = pymn + y1u - oset - 3;

           for(j=0; j < fp->obj_num; j++){

               obb = (fp->objs) + j;

               sx = sy = sstr = 0.0;

               for(k=0; k < obb->point_num; k++){

                   sx += obj_xreal((((obb->pt)+k)->x) + offx);
                   sy += *(gr->ygrid + (((obb->pt)+k)->y) + offy);
                   sstr += *(ob_field + ((((obb->pt)+k)->y) - 1)*xdim + (((obb->pt)+k)->x) - 1);

               }

               fpts = (ob->fet->fpt)+j;

               (fpts->x).xy = sx/obb->point_num;
               (fpts->y).xy = sy/obb->point_num;
               if(rdob == 'n') 
                 fpts->str = sstr / obb->point_num;
               else
                 fpts->str = *(aa + ((obb->pt->y) - 1)*xdim + (obb->pt->x) - 1);
               if(cc == 'v') fpts->str *= 0.5;

               free(obb->pt);

           }

           fp->obj_num = 0;
           fp->tot_f_f_num = 0;
           free(fp->objs);
           fp->objs = NULL;
           fp->obj_num = 0;
           free(hierarchy);

        }
   
        free(aa);
        free(ob_field);

        ob->ext = (struct extent * )malloc_initl(sizeof(struct extent));
        mem_er((ob->ext == NULL) ? 0 : 1, sizeof(struct extent));

        if(rdob == 'y'){

          pxmx = pxmn = ob->pt->x;
          pymx = pymn = ob->pt->y;

          for(j=0; j < ob->point_num; j++){

              ptt = (ob->pt) + j;

              if(ptt->x > pxmx) pxmx = ptt->x;
              else if(ptt->x < pxmn) pxmn = ptt->x;

              if(ptt->y > pymx) pymx = ptt->y;
              else if(ptt->y < pymn) pymn = ptt->y;

          }

        }

        ob->ext->x1 = pxmn;
        ob->ext->x2 = pxmx;
        ob->ext->y1 = pymn;
        ob->ext->y2 = pymx;        

    }

    return;

}
