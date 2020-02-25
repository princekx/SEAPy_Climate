#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "zones.h"


/* function to filter feature points with no chance of connection 
   to other feature points based on the upper bound on the possible
   displacements. */


extern int izz;
extern int nf, nfld;

float measure(struct feature_pts * , struct feature_pts * );
float dmaxx(ZONE * , struct feature_pts * , struct feature_pts * , float );

struct frame_objs *feature_pt_filter(ZONE **zone, struct frame_objs *fo, float *ddmax, int frame_num, int npar)
{

    int i, j, k, l, m, fdel;
    int iim1=0, iim2=0;

    float ndmax=*ddmax;
    float dmax=*ddmax;

    float dist;
    
    struct frame_objs *foo, *fo1=NULL, *fo2=NULL;
    struct object *ob, *ob1;
    struct feature_pts *fpts, *fpts1;
    

    for(i=0; i < frame_num; i++){

       foo = fo + i;
       iim1 = (foo->nmiss > npar - 1) ? npar - 1 : foo->nmiss;

       if(i > 0) {
         fo1 = foo - 1;
         iim2 = (fo1->nmiss > npar - 1) ? npar - 1 : fo1->nmiss;
       }
       if(i < frame_num-1) fo2 = foo + 1;

       for(j=0; j < foo->obj_num; j++){

          ob = (foo->objs) + j;

          for(k=0; k < ob->fet->feature_num; k++){

reset:

             fpts = (ob->fet->fpt) + k;

             fdel = 0;

             if(i < frame_num-1){

                dmax = ndmax = *(ddmax + iim1);

                for(l=0; l < fo2->obj_num; l++){

                    ob1 = (fo2->objs) + l;

                    for(m=0; m < ob1->fet->feature_num; m++){

                        fpts1 = (ob1->fet->fpt) + m;
                        dist = measure(fpts, fpts1);

                        if(izz) ndmax = dmaxx(zone[iim1], fpts, fpts1, dmax);

                        if(dist <= ndmax && fdel == 0) {fdel = 1; break;}
                    }

                    if(fdel == 1) break;

                 }

             }

             if(i > 0 && fdel == 0){

                dmax = ndmax = *(ddmax + iim2);

                for(l=0; l < fo1->obj_num; l++){

                    ob1 = (fo1->objs) + l;

                    for(m=0; m < ob1->fet->feature_num; m++){

                        fpts1 = (ob1->fet->fpt) + m;

                        dist = measure(fpts, fpts1);

                        if(izz) ndmax = dmaxx(zone[iim2], fpts, fpts1, dmax);

                        if(dist <= ndmax && fdel == 0) {fdel=1; break;}
                    }

                    if(fdel == 1) break;

                 }

             }

             if(fpts->str < CHECK_PT) fdel = 0;


             if(fdel == 0){

                l = 0;

                while(l+k < ob->fet->feature_num-1){

                     if(!l && nf) free((fpts+l)->add_fld);
                     *(fpts+l) = *(fpts + l + 1);
		     
                     ++l;

                }

                --(ob->fet->feature_num);
                --(foo->tot_f_f_num);

                ob->fet->fpt = (struct feature_pts * )realloc_n(ob->fet->fpt, ob->fet->feature_num * sizeof(struct feature_pts));
                mem_er((ob->fet->fpt == NULL) ? 0 : 1, ob->fet->feature_num * sizeof(struct feature_pts));

             }

             if(fdel == 0 && k < ob->fet->feature_num) goto reset;

          }

        }

    }


    return fo;

}
