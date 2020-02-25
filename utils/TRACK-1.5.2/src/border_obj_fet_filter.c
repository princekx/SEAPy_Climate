#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_im.h"
#include "st_obj.h"
#include "st_fo.h"

#define  ALL  'e'     /* 'a' to filter all feature points outside of region,
                         anything else to filter according to the amount of 
                         mass within region.                                */

/* function to filter out feature points for border objects according to
   the ratio of mass inside to the total mass of the object  .           */

extern int tf;
extern int xdm, ydm;
extern float fbr;
extern float xmn, ymn, xmx, ymx;

void border_obj_fet_filter(struct object *ob , double **fuz, struct frame_objs *fo )

{

   int i, j, k, l=0, m, al=ALL, ftag;
   int iptx, ipty, numpt;

   float ptx, pty;

   double divv, mply, intm, totm, *fz=NULL, **fzz=NULL;

   struct point *ptt;
   struct feature_pts *fpts, *fpts1, *fpts2;

   fzz = (double **)calloc(ob->fet->feature_num, sizeof(double *));
   mem_er((fzz == NULL) ? 0 : 1, (ob->fet->feature_num * sizeof(double *)));

   numpt = ob->fet->feature_num;

   if(al != 'a' && (tf == 1 || tf == 5 || tf == 6)) for(i=0; i < ob->fet->feature_num; i++) fzz[i] = *(fuz+i);

   for(i=0; i < ob->fet->feature_num; i++){

       m = 0;
       fpts = (ob->fet->fpt) + i;

reset:

       ftag = 0;

       if(tf == 3){

          iptx = ((fpts->x).ixy - 1)*(xdm - (fpts->x).ixy);
          ipty = ((fpts->y).ixy - 1)*(ydm - (fpts->y).ixy);

          if(iptx < 0 || ipty < 0) {ftag = 1; ++l;}

       }

       else if(al == 'a'){

          ptx = ((fpts->x).xy - xmn)*(xmx - (fpts->x).xy);
          pty = ((fpts->y).xy - ymn)*(ymx - (fpts->y).xy);

          if(ptx < 0 || pty < 0) {ftag = 1; ++l;}

       }

       else if(tf == 1 || tf == 5 || tf == 6){

          intm = totm = 0.0;

          for(j=0; j < ob->point_num; j++){

              fz = fzz[i] + j;

              ptt = (ob->pt)+j;           
              iptx = (ptt->x - 1)*(xdm - ptt->x);
              ipty = (ptt->y - 1)*(ydm - ptt->y);
              mply = (*fz) * ptt->val;

              if(iptx >= 0 && ipty >= 0) intm += mply;

              totm += mply;

           }

           divv = intm * 100.0/totm;

           if(divv < fbr) {ftag = 1; ++l;}

       }

       if(ftag == 1){

          ++m;
          k = 0;

          while(i+k < ob->fet->feature_num-m){

                fpts1 = fpts + k;
                fpts2 = fpts1 + 1; 
               
                (fpts1)->x = (fpts2)->x;
                (fpts1)->y = (fpts2)->y;
                (fpts1)->str = (fpts2)->str;

                if(al != 'a' && (tf == 1 || tf == 5 || tf == 6)) fzz[i+k] = fzz[i+k+1];

                ++k;

          }

          if(i+m < ob->fet->feature_num) goto reset;

       }

       if(m > 0) ob->fet->feature_num -= m;

   }

   if(l == numpt) {free(ob->fet->fpt);}

   else if(l > 0) {

      ob->fet->fpt = (struct feature_pts * )realloc_n(ob->fet->fpt, ob->fet->feature_num * sizeof(struct feature_pts));
      mem_er((ob->fet->fpt == NULL) ? 0 : 1, (ob->fet->feature_num * sizeof(struct feature_pts)));

   }

   fo->tot_f_f_num -= l;

   free(fzz);

   return;

} 
