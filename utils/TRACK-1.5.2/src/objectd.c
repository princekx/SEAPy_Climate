#include <Stdio.h>
#include <stdlib.h>
#include "st_fo.h"
#include "st_obj.h"
#include "boundary.h"
#include "mem_er.h"


/* function to write object data for each frame to an object dump file */

extern int objwr, fptwrt;
extern int aniso;

extern int nf, nfld;
extern int *nfwpos;

void objectd(struct frame_objs *ff, FILE *fobjo, int fnum, int tf)

{

     int i, j, k;
     struct object *ob=NULL;
     struct point *ptt=NULL;
     struct boundary_pt *btt=NULL;
     struct feature_pts *fp=NULL;

     fprintf(fobjo, "FRAME  %d  NMISS %d\n", fnum, ff->nmiss);
     fprintf(fobjo, "OBJ_NUM  %d\n", ff->obj_num);
     fprintf(fobjo, "BOUNDARY_STATE  %d\n", ff->b_state);
     fprintf(fobjo, "FRAME_FEATURE_NUM %d\n", ff->tot_f_f_num);

     if(abs(ff->b_state) > 1){

        printf("****ERROR****, boundary state not recognised. \r\n"
               "               no object data written to file.\n\n");
        exit(1);

     }

     for(i=0; i < ff->obj_num; i++){

         ob=(ff->objs)+i;

         fprintf(fobjo, "OBJECT  %d\n", i+1);

         if(objwr){
            fprintf(fobjo, "POINT_NUM  %d\n", ob->point_num);

            for(j=0; j < ob->point_num; j++){

                ptt=(ob->pt)+j;

                fprintf(fobjo, "%d %d %e \n", ptt->x, ptt->y, ptt->val);

            }


            fprintf(fobjo, "BOUNDARY_NUM  %d\n", ob->bound_num); 

            if(!(ff->b_state)){

               for(j=0; j < ob->bound_num; j++){

                  btt=(ob->bound)+j;

                  fprintf(fobjo, "%d %d %e \n", (btt->x).ixy, (btt->y).ixy, btt->val);

               }

            } 

            else {


               for(j=0; j < ob->bound_num; j++){

                  btt=(ob->bound)+j;

                  fprintf(fobjo, "%f %f %e \n", (btt->x).xy, (btt->y).xy, btt->val);
               }

            }  

         }

         else {
            fprintf(fobjo, "POINT_NUM  %d\n", 0);
            fprintf(fobjo, "BOUNDARY_NUM  %d\n", 0);
         }


         if(tf && ob->fet){

            fprintf(fobjo, "OBJECT_FEATURE_NUM  %d\n", ob->fet->feature_num);

            if(aniso == 'y'){

               for(j=0; j < ob->fet->feature_num; j++){

                  fp = (ob->fet->fpt)+j;

/* initialize track_id's */

                  fp->track_id = 0;

                  fprintf(fobjo, "%f %f %e %f %f %f %e ", (fp->x).xy, (fp->y).xy, fp->str, fp->r_sh, fp->ornt[0], fp->ornt[1], fp->area);
                  if(nf) {
                    fprintf(fobjo, "& ");
                    for(k=0; k < nfld; k++) fprintf(fobjo, "%e & ", *(fp->add_fld + k));
                  }
                  fprintf(fobjo, "\n");

               }

            }

            else {

               for(j=0; j < ob->fet->feature_num; j++){

                  fp = (ob->fet->fpt)+j;

/* initialize track_id's */

                  fp->track_id = 0;

                  if(tf == 3) fprintf(fobjo, "%d %d %e ", (fp->x).ixy, (fp->y).ixy, fp->str);

                  else fprintf(fobjo, "%f %f %e ", (fp->x).xy, (fp->y).xy, fp->str);
                  if(nf) {
                    fprintf(fobjo, "& ");
                    for(k=0; k < nfld; k++) fprintf(fobjo, "%e & ", *(fp->add_fld + k));
                  }
                  fprintf(fobjo, "\n");


               }

            }

         }

         else {

            fprintf(fobjo, "OBJECT_FEATURE_NUM  %d\n", 0);

         }

         if(fptwrt) {

            if(nf){
               for(j=0; j < ob->fet->feature_num; j++) {fp = (ob->fet->fpt)+j; free(fp->add_fld);}
            }
            free(ob->fet->fpt);
            free(ob->fet);
            free(ob->ext);


         }

         ob->point_num = 0;
         free(ob->pt);    
         ob->bound_num = 0;
         free(ob->bound);         

      }

      if(!tf) free(ff->objs);
      if(fptwrt) free(ff->objs);


      fflush(fobjo);


      return;

}
