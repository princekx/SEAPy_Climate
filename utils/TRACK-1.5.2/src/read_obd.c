#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "st_obj.h"
#include "st_fo.h"
#include "boundary.h"
#include "mem_er.h"
#include "file_handle.h"

#define  STRLEN   200

extern int tf, fruser;
extern int objwr;
extern int ffirst, fint, frterm;
extern int aniso;
extern int x1u, x2u, y1u, y2u;
extern int irmiss;

extern int nf, nfld;
extern int *nfwpos;

int wc(char * );


/* function to read in existing object and feature data for re-use */

struct frame_objs *read_obd(FILE *fobjo, float *alat, float *alng, int *pgr, int *pty)

{


    int xx1=x1u, xx2=x2u, yy1=y1u, yy2=y2u;
    int i=0, j=0, k=0;
    int fr_id, ob_id, nmiss=0;
    int ofl=0;
    int nf0, nfld0;
    int iadd=0;
    
    int *nfwpos0=NULL;

    off_t plt, pltn;

    char tex[STRLEN], tex2[STRLEN];
    char *fldt=NULL;

    struct frame_objs *foo=NULL, *fp=NULL;
    struct object *ob=NULL;
    struct point *ptt=NULL;
    struct boundary_pt *btt=NULL;
    struct feature_pts *fpts=NULL;


    fruser = 0;

    fgets(tex, STRLEN, fobjo);

    sscanf(tex, "%d %d %d", &tf, &aniso, &irmiss);
    if(aniso == 1) aniso = 'y';
    else aniso = 'n';

    if(!tf){

       printf("***WARNING***, this file contains threshold data only\r\n"
              "               there are no feature points.          \n\n");
    }

    fgets(tex, STRLEN, fobjo);
    fgets(tex, STRLEN, fobjo);
    sscanf(tex, "%d %d", pgr, pty);
    if(*pgr) {
      fgets(tex, STRLEN, fobjo);
      sscanf(tex, "%f %f", alat, alng);
    }

    plt = ftello(fobjo);
    fgets(tex, STRLEN, fobjo);
    if(strstr(tex, "REGION_DETAILS")) {
       fgets(tex2, STRLEN, fobjo);
       sscanf(tex2, "%d %d %d %d", &xx1, &xx2, &yy1, &yy2);
       pltn = ftello(fobjo);
       fgets(tex, STRLEN, fobjo);
       if(strstr(tex, "ADDITIONAL_FIELDS")) iadd = 1;
       else fseeko(fobjo, pltn, FSTART);
    }
    else if(strstr(tex, "ADDITIONAL_FIELDS")) iadd = 1;
    else fseeko(fobjo, plt, FSTART);
    
    if(iadd){
       fgets(tex2, STRLEN, fobjo);
       sscanf(tex2, "%d %d", &nf0, &nfld0);
       nfwpos0 = (int *)calloc(nf0, sizeof(int));
       mem_er((nfwpos0 == NULL) ? 0 : 1, nf0 * sizeof(int));
       fldt = strchr(tex2, '&') + 1;
       for(i=0; i < nf0; i++) {sscanf(tex2, "%1d", nfwpos0 + i); ++fldt;}
       
       if(nf){
          if(nf != nf0 || nfld != nfld0){
             printf("****ERROR****, additional field information does not match in %s at line %d.\n\n", __FILE__, __LINE__);
             exit(1);
          }
          for(i=0; i < nf; i++){
              if(*(nfwpos + i) != *(nfwpos0 + i)){
                 printf("****ERROR****, additional field information does not match in %s at line %d.\n\n", __FILE__, __LINE__);
                 exit(1);
              }
          } 
       }
       else {
          nf = nf0;
          nfld = nfld0;
          nfwpos = (int *)calloc(nf, sizeof(int));
          mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));
          memcpy(nfwpos, nfwpos0, nf0 * sizeof(int));
       }
    }

    if(((xx1 != x1u) || (xx2 != x2u)) || ((yy1 != y1u) || (yy2 != y2u))){

       printf("****WARNING****, chosen region does not match original region.\r\n"
              "                 use the same input parameters that created   \r\n"
              "                 this object file.                            \n\n");

    }


    while(fgets(tex2, STRLEN, fobjo) != NULL && !feof(fobjo)){

          if(!strstr(tex2, "FRAME")){
             printf("****ERROR****, reading object file,       \r\n"
                    "               no FRAME tag, at Frame %d, \r\n"
                    "               exiting.\n\n", fr_id + 1);
             exit(1);

          }

          sscanf(tex2, "%*s %d %*s %d", &fr_id, &nmiss);

          ++fruser;

          if(fruser == 1){

             foo = (struct frame_objs * )malloc_initl(sizeof(struct frame_objs));
             mem_er((foo == NULL) ? 0 : 1, sizeof(struct frame_objs));

             fp = foo;
             ffirst = fr_id;

          }

          else{

             foo = (struct frame_objs * )realloc_n(foo, fruser*sizeof(struct frame_objs));
             mem_er((foo == NULL) ? 0 : 1, fruser*sizeof(struct frame_objs));

             fp = foo + fruser - 1;

          }

          fp->b_state = -1;
          fp->frame_id = fr_id;
          fp->nmiss = nmiss;

          fgets(tex2, STRLEN, fobjo);

          if(!strstr(tex2, "OBJ_NUM")){

             printf("****ERROR****, reading object file,           \r\n"
                    "               no OBJ_NUM tag, at Frame %d,\r\n"
                    "               exiting.\n\n", fr_id);
             exit(1);

          }

          sscanf(tex2, "%*s %d", &fp->obj_num);
         
          plt = ftello(fobjo);
          fgets(tex2, STRLEN, fobjo);

          if(strstr(tex2, "BOUNDARY_STATE")) sscanf(tex2, "%s %d", tex, &(fp->b_state));

          else fseeko(fobjo, plt, FSTART);

          if(abs(fp->b_state) > 1){

             printf("****ERROR****, boundary state not recognised. \r\n"
                    "               no object data read from file.\n\n");
             exit(1);

          }

          fgets(tex2, STRLEN, fobjo);

          if(!strstr(tex2, "FRAME_FEATURE_NUM")){

             printf("****ERROR****, reading object file,                  \r\n"
                    "               no FRAME_FEATURE_NUM tag, at Frame %d,\r\n"
                    "               exiting.\n\n", fr_id);
             exit(1);

          }

          sscanf(tex2, "%*s %d", &fp->tot_f_f_num);

          fp->objs = (struct object * )calloc(fp->obj_num, sizeof(struct object));
          mem_er((fp->objs == NULL) ? 0 : 1, fp->obj_num * sizeof(struct object));

          for(i=0; i < fp->obj_num; i++){

              ob = (fp->objs) + i;
              ob->bound_num = 0;
              ob->point_num = 0;
              ob->bound = NULL;
              ob->fet = NULL;
              ob->ext = NULL;

              ob->fet = (struct features * )malloc_initl(sizeof(struct features));
              mem_er((ob->fet == NULL) ? 0 : 1, sizeof(struct features));

              fgets(tex2, STRLEN, fobjo);
              sscanf(tex2, "%s %d", tex, &ob_id);

              fgets(tex2, STRLEN, fobjo);
              sscanf(tex2, "%s %d", tex, &ob->point_num);

              if(objwr) {

                 ob->pt = (struct point * )calloc(ob->point_num, sizeof(struct point));
                 mem_er((ob->pt == NULL) ? 0 : 1, ob->point_num * sizeof(struct point));

                 for(j=0; j < ob->point_num; j++){

                     ptt = (ob->pt) + j;
                     fgets(tex2, STRLEN, fobjo);
                     sscanf(tex2, "%d %d %e", &ptt->x, &ptt->y, &ptt->val);
                 }

              }

              else{

                 for(j=0; j < ob->point_num; j++){

                     ptt = (ob->pt) + j;
                     fgets(tex2, STRLEN, fobjo);

                 }

                 ob->point_num = 0;

              }

              plt = ftello(fobjo);
              fgets(tex2, STRLEN, fobjo);

              if(strstr(tex2, "BOUNDARY_NUM")){

                 sscanf(tex2, "%s %d", tex, &(ob->bound_num));

                 if(objwr && ob->bound_num > 0){

                    ob->bound = (struct boundary_pt * )calloc(ob->bound_num, sizeof(struct boundary_pt));
                    mem_er((ob->bound == NULL) ? 0 : 1, ob->bound_num * sizeof(struct boundary_pt));

                    if(!(fp->b_state)){

                       for(j=0; j < ob->bound_num; j++){

                           btt = (ob->bound) + j;
                           fgets(tex2, STRLEN, fobjo);
                           sscanf(tex2, "%d %d %e", &(btt->x).ixy, &(btt->y).ixy, &(btt->val));


                       }

                    }

                    else{

                       for(j=0; j < ob->bound_num; j++){

                           btt = (ob->bound) + j;
                           fgets(tex2, STRLEN, fobjo);
                           sscanf(tex2, "%f %f %e", &(btt->x).xy, &(btt->y).xy, &(btt->val));

                       }


                    }

                 }

                 else{

                    if(!(fp->b_state)){

                       for(j=0; j < ob->bound_num; j++){

                           btt = (ob->bound) + j;
                           fgets(tex2, STRLEN, fobjo);

                       }

                    }

                    else{

                       for(j=0; j < ob->bound_num; j++){

                           btt = (ob->bound) + j;
                           fgets(tex2, STRLEN, fobjo);

                       }


                    }

                    ob->bound_num = 0;


                 }

              }

              else fseeko(fobjo, plt, FSTART);

              fgets(tex, STRLEN, fobjo);
              sscanf(tex, "%*s %d", &ob->fet->feature_num);

              ob->fet->fpt = (struct feature_pts * )calloc(ob->fet->feature_num, sizeof(struct feature_pts));
              mem_er((ob->fet->fpt == NULL) ? 0 : 1, ob->fet->feature_num * sizeof(struct feature_pts));

              if(aniso == 'y'){

                 for(j=0; j < ob->fet->feature_num; j++){

                    fpts = (ob->fet->fpt) + j;

                    fgets(tex2, STRLEN, fobjo);

                    ofl = sscanf(tex2, "%f %f %e %f %f %f %e", &(fpts->x).xy, &(fpts->y).xy, &(fpts->str), &(fpts->r_sh), &(fpts->ornt[0]), &(fpts->ornt[1]), &(fpts->area));

                    if(nf){
                    
                       fpts->add_fld = (float *)calloc(nfld, sizeof(float));
                       mem_er((fpts->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                       fldt = strchr(tex2, '&') + 1;
                       for(k=0; k < nfld; k++) {
                           sscanf(fldt, "%e", fpts->add_fld + k);
                           fldt = strchr(fldt, '&') + 1;
                       }
                    }

                 }
              }

              else {

                 for(j=0; j < ob->fet->feature_num; j++){

                    fpts = (ob->fet->fpt) + j;

                    if(tf == 3){
                       fgets(tex2, STRLEN, fobjo);
                       sscanf(tex2, "%d %d %e", &(fpts->x).ixy, &(fpts->y).ixy, &(fpts->str));

                    }
                    else {
                       fgets(tex2, STRLEN, fobjo);
                       sscanf(tex2, "%f %f %e", &(fpts->x).xy, &(fpts->y).xy, &(fpts->str));
                       
                       if(nf){
                          fpts->add_fld = (float *)calloc(nfld, sizeof(float));
                          mem_er((fpts->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                          fldt = strchr(tex2, '&') + 1;
                          for(k=0; k < nfld; k++) {
                             sscanf(fldt, "%e", fpts->add_fld + k);
                             fldt = strchr(fldt, '&') + 1;
                          }
                       }

                    }

                 }

              }
 
          }

    }

    if(aniso == 'y' && ofl != 7){

       printf("****WARNING****, object file is an older version with no area data.\n\n");

    }

    frterm = fr_id;
    if(fruser > 1){
      if((frterm - ffirst) % (fruser - 1)){
         printf("***WARNING***, reading object data file, incorrect sequence of frames\r\n");

      }
      fint = (frterm - ffirst) / (fruser - 1);
    }
    else fint = 1;
    
    if(nfwpos0) free(nfwpos0);


    return foo;

}
              



          
