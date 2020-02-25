#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include "st_obj.h"
#include "st_fo.h"
#include "mem_er.h"
#include "m_values.h"
#include "grid.h"
#include "proj.h"

/* function to remove close by maxima */

float measure(struct feature_pts * , struct feature_pts * );
void geo_conv(struct feature_pts * );

extern int tom;
extern GRID *gr;
extern int x1u, y1u;
extern int geo_init;
extern PROJ *pp;

struct frame_objs *mfilt_obj_dist(struct frame_objs *fo, float rads, int tf)
{

   int j, k, m;
   int imx=0;
   int nstr=0;
   int nfet=0;
   int ift=0;

   float fmax=0.0;
   float *strtmp=NULL;
   float xa, ya;

   struct features *fet=NULL;
   struct object *ob=NULL;
   struct feature_pts *fpt=NULL, *fptm=NULL, *fptn=NULL;
   struct feature_pts fp1;


   if(tom == 'g') geo_init = 1;

   fo->tot_f_f_num = 0;

   for(j=0; j < fo->obj_num; j++){

      ob = fo->objs + j;
      fet = ob->fet;

      if(fet->feature_num > 1){

/* put intensities into tempory array */

         nstr = fet->feature_num;
         nfet = 0;

         strtmp = (float *)calloc(fet->feature_num, sizeof(float));
         mem_er((strtmp == NULL) ? 0 : 1, (fet->feature_num)*sizeof(float));

         for(k=0; k < fet->feature_num; k++){
            fpt = fet->fpt + k;
            *(strtmp + k) = fpt->str;
         }

/* convert to cartesians for distance measure */

         if(gr->prty && tom == 'g'){

            for(k=0; k < fet->feature_num; k++){
                fpt = fet->fpt + k;
                memcpy(&fp1, fpt, sizeof(struct feature_pts));
                  if(tf == 3){
                     xa = *(gr->xgrid + ((fp1.x).ixy)+x1u-2);
                     ya = *(gr->ygrid + ((fp1.y).ixy)+y1u-2);
                     switch(gr->prgr){
                       case 0:
                          (pp->prj2)(&xa, 'x', 0);
                          (pp->prj2)(&ya, 'y', 0);
                          break;
                       case 1:
                          azimuthal(&ya, &xa, gr->alat, gr->alng, 0, gr->prty);
                          break;
		       case 2:
		          conic(&ya, &xa, gr->alat, gr->alng, gr->sp1, gr->sp2, 0, gr->prty, &ift);
			  break;
		       case 3:
		          rotated_cyln(&ya, &xa, gr->alat, gr->alng, 0, &ift);
			  break;
                     }
                     (fp1.x).xy = xa;
                     (fp1.y).xy = ya;
                  }
                  else {
                     switch(gr->prgr){
                       case 0:
                         (pp->prj2)(&((fp1.x).xy), 'x', 0);
                         (pp->prj2)(&((fp1.y).xy), 'y', 0);
                         break;
                       case 1:
                         azimuthal(&((fp1.y).xy), &((fp1.x).xy), gr->alat, gr->alng, 0, gr->prty);
                         break;
		       case 2:
		         conic(&((fp1.y).xy), &((fp1.x).xy), gr->alat, gr->alng, gr->sp1, gr->sp2, 0, gr->prty, &ift);
			 break;
		       case 3:
		         rotated_cyln(&((fp1.y).xy), &((fp1.x).xy), gr->alat, gr->alng, 0, &ift);
			 break;
                     }
                  }

                  geo_conv(&fp1);
                  fpt->gwk[0] = fp1.gwk[0];
                  fpt->gwk[1] = fp1.gwk[1];
                  fpt->gwk[2] = fp1.gwk[2];
                  fpt->gwky = 1;
            }

         }

         else if(tom == 'g'){
            for(k=0; k < fet->feature_num; k++){
                fpt = fet->fpt + k;
                geo_conv(fpt);
            }
         }

         while(nstr){

/* find maxima in order largest first */

            fmax = CHECK_PT;
            for(k=0; k < fet->feature_num; k++){
                if(*(strtmp + k) > fmax){fmax = *(strtmp + k); imx = k;}
            }

            if(fmax > CHECK_PT){

               fptm = fet->fpt + imx;

               for(k=0; k < fet->feature_num; k++){
                  if(k == imx) continue;
                  fpt = fet->fpt + k;
                  if(*(strtmp + k) > CHECK_PT){

                     if(measure(fptm, fpt) < rads) {
                        *(strtmp + k) = DUFF_PT; 
                        fpt->str = DUFF_PT;
                        --nstr;
                     }
                  }
               }

               *(strtmp + imx) = DUFF_PT;
            }
            --nstr;
         }

         for(k=0; k < fet->feature_num; k++){
             fpt = fet->fpt + k;
             if(fpt->str < CHECK_PT){
                fptn = fpt;
                fptm = fpt + 1;
                for(m=k; m < fet->feature_num - 1; m++){
                    memcpy(fptn, fptm, sizeof(struct feature_pts));
                    ++fptn;
                    ++fptm;
                }
                --k;
                --(fet->feature_num);
             }
             else ++nfet;
         }

         fet->fpt = (struct feature_pts * )realloc_n(fet->fpt, nfet * sizeof(struct feature_pts));
         mem_er((fet->fpt == NULL) ? 0 : 1, nfet*sizeof(struct feature_pts));

         fet->feature_num = nfet;
         fo->tot_f_f_num += nfet;

         free(strtmp);


      }

      else if(fet->feature_num > 0) ++(fo->tot_f_f_num);

   }

   if(tom == 'g') geo_init = 0;

   return fo;

}
