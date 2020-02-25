#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include <sys/types.h>
#include "file_cat_out.h"
#include "proj.h"
#include "mem_er.h"
#include "st_obj.h"
#include "st_fo.h"
#include "st_track.h"
#include "file_handle.h"
#include "grid.h"
#include "splice.h"

#define  TOL      1.0e-6
#define  UTOL     1.0
#define  BDIST    1000.0
#define  MMSTR    1.0e+10

#define  LPERC     0.5    /* percentage of track allowed outside chosen frame range */

#define  IDEV     0       /* switch for outputing the track deviation as a function */
                          /* displacement for calibrating tracking parameters from  */
                          /* sub-sampled tracks.                                    */

/* function to read in track data from different files and combine them */

float measure(struct feature_pts * , struct feature_pts * );
struct frame_objs *read_obd(FILE * , float * , float * , int * , int * );
void meantrd(FILE * , struct tot_tr * , int , int );
struct track_ind *read_td(FILE * , int);
struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * );
void splice_plot(struct tot_tr * , int , int ,int , int , int );
void disp_filter(struct tot_tr * , int , int * );
void prop_filt(struct tot_tr * , int , int * );
void speed_filt(struct tot_tr * , int , int * );
void track_sample_dev(struct tot_tr* , int );
void track_merge(struct tot_tr * , int , int );
void write_track_netcdf(struct tot_tr *  , int , char *, int );

/* void fourier_match(struct tot_tr * , int , int * ); */

char *trfil=NULL;

int trtyp='s';

extern int tf, fruser, track_num;
extern int x1u, y1u, x2u, y2u;
extern GRID *gr, *gr1;
extern CNTRY *cm, *cm1;

extern float xmn, ymn, xmx, ymx;

extern int aniso;

extern char *fext;

extern int iext;

extern int nfld, nf;
extern int *nfwpos;

struct tot_tr *splice(int *trackn)

{

    int izm, mnpt, mxpt, ff=0;
    int i, j, k, tt=0;
    int tot_fet=0, tr_count=0;
    int dnum, dnm=0;
    int tplt, fr_tot=0;
    int asmm=0, am, jf, rp=0;
    int ex;
    int tpr, gpr;
    int tn=0;
    int lframe=0;
    int anis_t='n';
    int strsw;
    int ntrnum=0;
    int ifr1, ifr2;
    int ifb, ife, nmiss, iout;
    int id1=0, id2=0, idd=0;
    int numt=0;
    int strty=0;
    int iadd=0;
    int igrwth=0;
    int iff=0, ifloc=0;
    int imiss=0;
    int lpos=0, lneg=0, iextm=0;
    int ld1=0, ld2=0;
    int tstep=0;
    
    int tnum=0, trty=0;

    off_t pl=0;

    float alat=gr->alat, alng=gr->alng;

    float strth1, strth2, str, str1, str2;

    float dd, ddf;
    float mmstr=0.0;

    double diff;

    char charin[MAXCHR];
    char trout[MAXCHR];
    char grtrs[MAXCHR];

    FILE *fobj=NULL, *ftr=NULL, *tsf=NULL;

    struct fet_pt_tr *atr=NULL, *atr2=NULL, *at=NULL;
    struct tot_tr *all_tr=NULL, *altr=NULL;
    struct track_points *tps=NULL;
    struct feature_pts *fpts=NULL, fpts1, fpts2;
    struct object *ob=NULL;

    struct frame_objs *fo, *fp=NULL;          /* pointer to feature points */
    struct track_ind *tind=NULL, *tr=NULL;    /* pointer to track points   */

    trtyp = 's';

    printf("do you want to use an existing combined track data file, 'y' or 'n'\n");

    scanf("\n");
    ex = getchar();

reread:

    if(ex == 'y'){

        trfil = (char *)malloc_initl(MAXCHR*sizeof(char));
        mem_er((trfil == NULL) ? 0 : 1, MAXCHR*sizeof(char));

/* read combined data from file */
        strncpy(charin, FPTTRS, MAXCHR);
        if(iext) strcpy(strstr(charin, EXTENSION), fext);

       printf("Do you want to use a track file different from the default, 'y' or 'n'.\n\n");
       scanf("\n");
       if(getchar() == 'y'){
          printf("what is the name of the combined file to be read\n");

          scanf("%s", charin);

          printf("What is the track file type, i.e. 's' or 'v'\n\n");
          scanf("\n");
          trtyp = getchar();

          if(!(trtyp == 's' || trtyp == 'v')){
             printf("****WARNING****, incorrect track file type defaulting to 's'.\n\n");
             trtyp = 's';
          }

       }

       strncpy(trfil, charin, MAXCHR);

       tsf = open_file(charin, "r");

       all_tr = read_tracks(tsf, &tr_count, &gpr, &tpr, trtyp, &gr->alat, &gr->alng);

       close_file(tsf, charin);

       if(gpr){alat = gr1->alat; alng = gr1->alng;}
       
       if(!(gr->iaz)){
          if((gpr != gr->prgr) || (tpr != gr->prty)  ||
             (fabs(gr->alat - alat) > TOLPROJ)       ||
             (fabs(gr->alng - alng) > TOLPROJ)         ) {

             printf("***WARNING***, map projection and data do not match\r\n"
                    "               grid will be corrected\n\n");


             proj_group(gpr, tpr);


          }
       }

       goto plot;

    }


    printf("do you want track data set addition or matching,\r\n"
           "input '0' for addition and '1' for matching\n");
    scanf("%d", &am);

    if(am == 1){

      printf("**************************************************************\n\n");
      printf("   W     W     W    A     RRR   N  N  III  N  N  GG           \n");
      printf("    W   W W   W    A A    R  R  NN N   I   NN N G             \n");
      printf("     w w   w w    AAAAA   RRR   N NN   I   N NN G  GG         \n");
      printf("      W     w    A     A  R  R  N  N  III  N  N  GG           \n\n");
      printf("**************************************************************\n\n");
      printf("data sets must overlap by at least one frame for the matching\r\n"
             "to be carried out correctly. This is the users responsibility!!\n");
      printf("**************************************************************\n\n");

      printf("***Warning*** Make shure that data sets are read in\r\n"
             "              chronological order if matching is required!!\n\n");

    }
    
    if(nf){
       nf = 0;
       nfld = 0;
       free(nfwpos);
    }

    printf("How many sets of data are to be combined\n\n");
    scanf("%d", &dnum);

    while(++dnm <= dnum){

        printf("What is the %d object file to be read\n\n", dnm);

        scanf("%s", charin);

        fobj = open_file(charin, "r");


        fo = read_obd(fobj, &gr->alat, &gr->alng, &gpr, &tpr);

        close_file(fobj, charin);

        if(dnm > 1 && anis_t != aniso){

           printf("****WARNING****, object files have incompatable shape descriptor tag. \n\n");

        }

        anis_t = aniso;

        if(gpr){alat = gr1->alat; alng = gr1->alng;}

        if((gpr != gr->prgr) || (tpr != gr->prty) ||
          (fabs(gr->alat - alat) > TOLPROJ)  ||
          (fabs(gr->alng - alng) > TOLPROJ)  ) {

           printf("***WARNING***, map projection and data do not match\r\n"
                  "               grid will be corrected\n\n");


           proj_group(gpr, tpr);

           if(gpr) { 

              gr = gr1;
              cm = cm1;

              x1u = 1;
              x2u = gr->ix;
              y1u = 1;
              y2u = gr->iy;

              xmn = *(gr->xgrid + x1u-1);
              ymn = *(gr->ygrid + y1u-1);
              xmx = *(gr->xgrid + x2u-1);
              ymx = *(gr->ygrid + y2u-1);

              printf("***WARNING***, region extent must be the same as when\r\n"
                     "               data file was created, check original \r\n"
                     "               input. Otherwise errors may result.  \n\n");


           }


        }

        if(dnm == 1) tt = tf;

        else if(tt != tf){

          printf("***Warning***, data in these files may not be compatable,\r\n"
                 "               different feature point detection methods \r\n"
                 "               has been used.\n\n");

        }

/* read in corresponding track data */

        printf("What is the %d track file to be read: MUST CORRESPOND TO OBJECT DATA FILE\n\n", dnm);
        scanf("%s", charin);

        ftr = open_file(charin, "r");

        tind = read_td(ftr, fruser);

        close_file(ftr, charin);

        printf("What is the first frame (> 0) for this data set\n");
        scanf("%d", &jf);

retry:

        iff = 0;

        if(jf > fruser) {

           printf("***error***, illegal value for starting frame, number\r\n"
                  "             of frames for this data set is %d, try again\n", fruser);
           iff = 1;
        }

        if(iff == 1) goto retry;

/* assign storage for combined tracks */

        for(i=0; i < track_num; i++){

           tr = tind + i;

           rp = 0;

           for(j=0; j < jf-1; j++)

              if(((tr->tp)+j)->feature_id > 0) ++rp;

/* what are the total number of feature points */

            if(dnm == 1){

              ++tr_count;
              asmm = 2;

              if(tr_count == 1){

                 all_tr = (struct tot_tr * )malloc_initl(sizeof(struct tot_tr));
                 mem_er((all_tr == NULL) ? 0 : 1, sizeof(struct tot_tr));

                 all_tr->trpt = (struct fet_pt_tr * )calloc(tr->num_point_track-rp, sizeof(struct fet_pt_tr));
                 mem_er((all_tr->trpt == NULL) ? 0 : 1, tr->num_point_track-rp * sizeof(struct fet_pt_tr));

                 all_tr->num = tr->num_point_track-rp;

                 altr = all_tr;
                 atr = all_tr->trpt;


              }

              else asmm = 1;

            }

            else if(am == 1){

              asmm = 1;

              tps = (tr->tp) + jf - 1;
              fp = fo + (tps->frame_id)-1;
              ob = (fp->objs) + (tps->object_id) - 1;

              if(tps->feature_id != -1){

                 fpts = (ob->fet->fpt) + (tps->feature_id) - 1;

                 for(j=0; j < tr_count; j++){

                     altr = all_tr + j;

                     diff = UTOL;

                     if(altr->eofs && altr->num > 0){

                        atr = (altr->trpt) + altr->num - 1;
                        if(atr->fr_id != fr_tot) continue;
                        (fpts1.x).xy = atr->xf;
                        (fpts1.y).xy = atr->yf;

                        if(tf == 3){

                           tf = 4;
                           (fpts2.x).xy = *(gr->xgrid + ((fpts->x).ixy)+x1u-2);
                           (fpts2.y).xy = *(gr->ygrid + ((fpts->y).ixy)+y1u-2);
                           diff = measure(&fpts2, &fpts1);
                           tf = 3;     

                        }

                        else diff = measure(fpts, &fpts1);

                     }

                     else if(altr->eofs && altr->num <= 0) {

                           printf("****WARNING****, possible track miss-alignment for track %d in %s\n", j+1, __FILE__);


                     }

                     if(diff < TOL) {

                        altr->trpt = (struct fet_pt_tr * )realloc_n(altr->trpt, (altr->num + tr->num_point_track - rp -1)*sizeof(struct fet_pt_tr));
                         mem_er((altr->trpt == NULL) ? 0 : 1, (altr->num + tr->num_point_track - rp -1)*sizeof(struct fet_pt_tr));
                        asmm = 0;

                        break;

                     }
                   
                 }

              }

              if(asmm == 1) ++tr_count;
              
            }

            else{asmm = 1; ++tr_count;}

            if(asmm == 1){

              all_tr = (struct tot_tr * )realloc_n(all_tr, tr_count*sizeof(struct tot_tr));
              mem_er((all_tr == NULL) ? 0 : 1, tr_count * sizeof(struct tot_tr));


              (all_tr+tr_count-1)->trpt = (struct fet_pt_tr * )calloc(tr->num_point_track - rp, sizeof(struct fet_pt_tr));
              mem_er(((all_tr+tr_count-1)->trpt == NULL) ? 0 : 1, (tr->num_point_track - rp) * sizeof(struct fet_pt_tr));

              altr = all_tr + tr_count - 1;
              altr->awt = 0;
              altr->time = 0;
              altr->num = tr->num_point_track - rp;
              atr = altr->trpt;

            }

            else if(asmm == 0) {

              atr = (altr->trpt) + altr->num;
              altr->num += tr->num_point_track - rp -1;

            }

            ff = (am == 0) ? jf -1 : ((am*asmm > 0) ? jf-1 : jf);

            lframe = 0;

            for(j=ff; j < fruser; j++){

                tps = (tr->tp) + j;

                if(tps->feature_id != -1){

                   at = atr++;

                   lframe = tps->frame_id;

                   fp = fo + (tps->frame_id)-1;
                   ob = (fp->objs) + (tps->object_id) - 1;
                   fpts = (ob->fet->fpt) + (tps->feature_id) - 1;

                   at->fr_id = fr_tot + tps->frame_id - ((am == 0 || dnm == 1) ? jf-1 : jf);

                   at->obj_id = tps->object_id;

                   at->zf = fpts->str;
                   at->nfm = tps->nmpt;

                   at->sh_an = fpts->r_sh;
                   at->or_vec[0] = fpts->ornt[0];
                   at->or_vec[1] = fpts->ornt[1];
                   at->area = fpts->area;
		   
		   if(nf){
                      at->add_fld = (float *)calloc(nfld, sizeof(float));
                      mem_er((at->add_fld == NULL) ? 0 : 1, nfld * sizeof(float));
                      memcpy(at->add_fld, fpts->add_fld, nfld * sizeof(float)); 
                   }
                   else at->add_fld = NULL;

                   if(tf == 3){

                      at->xf = *(gr->xgrid + ((fpts->x).ixy)+x1u-2);
                      at->yf = *(gr->ygrid + ((fpts->y).ixy)+y1u-2);

                    }

                    else {

                      at->xf = (fpts->x).xy;
                      at->yf = (fpts->y).xy;

                    }

                }

            }


            altr->eofs = 0;

            if(lframe == fruser) altr->eofs = dnm;

        }

        if(dnm > 1){

           for(i=0; i < tr_count; i++){
              altr = all_tr + i;
 
              if(altr->eofs < dnm) altr->eofs = 0;
           }

        }

        fr_tot += ((am == 0 || dnm == 1) ? fruser - jf + 1 : fruser - jf);

        for(i=0; i < fruser; i++){

             fp = fo + i;

             for(j=0; j < fp->obj_num; j++){

                 ob = (fp->objs) + j;
                 free(ob->pt);
	         if(nf){
                    for(k=0; k < ob->fet->feature_num; k++){fpts = (ob->fet->fpt)+k; free(fpts->add_fld);}
                 } 
                 free(ob->fet->fpt);
                 free(ob->fet);
                 free(ob->ext);
                 free(ob->bound);

             }

             free(fp->objs);

         }

         free(fo);

         for(i=0; i < track_num; i++){

            tr = tind + i;
            free(tr->tp);

         }
         free(tind);

    }
    
/* set track id's */
    
    for(i=0; i < tr_count; i++){
        altr = all_tr + i;
	altr->trid = i + 1;
    }

    printf("Do you want to check for split tracks based on common object Id's and merge them, 'y' or 'n'?\n\n");
    scanf("\n");
    if(getchar() == 'y') track_merge(all_tr, tr_count, 's');

/* write combined data to file */

    strncpy(trout, FPTTRS, MAXCHR);
    if(iext) strcpy(strstr(trout, EXTENSION), fext);

    tsf = open_file(trout, "w");

    meantrd(tsf, all_tr, tr_count, 's');  
 
    close_file(tsf, trout);

    write_track_netcdf(all_tr, tr_count, trout, 's');

plot:

    printf("***INFORMATION***, you can now filter the tracks according  \r\n"
           "                   to lifetime and/or displacement distance,\r\n"
           "                   and/or maximum strength (i.e. system must\r\n"
           "                   reach this strength to be included. Also,\r\n"
           "                   it is possible to filter according to    \r\n"
           "                   mean propgation direction.               \n\n");

    printf( "input the min. and max. number of points in a track\n\n");

    scanf("%d %d", &mnpt, &mxpt);

    if(mnpt < 0 || mxpt < 0){
       printf("****ERROR****, min. and max. points must be positive,\r\n"
              "               try again.                            \n\n");
       goto plot;
    }

/* check each track for lifetime */

    tot_fet = 0;

    tn = 0;

    for(i=0; i < tr_count; i++){

        altr = all_tr + i;

        numt = altr->num;

        for(j=0; j < altr->num; j++) numt += (altr->trpt + j)->nfm;

        if(numt < mnpt || numt > mxpt){

	   for(j=0; j < altr->num; j++) {if(altr->trpt) free((altr->trpt + j)->add_fld);}
           altr->num = 0;
           free(altr->trpt);
           altr->trpt = NULL;

        }

        else ++tn;

        tot_fet += altr->num;

    }

    printf("***INFORMATION***, current number of tracks is %d\n\n", tn);

    printf("Do you want to filter according to distance, 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y') disp_filter(all_tr, tr_count, &tot_fet);

    printf("Do you want to filter according to system strength, 'y' or 'n'.\n\n");

    scanf("\n");
    if(getchar() == 'y'){

       if(nf && trtyp != 'v'){
            printf("There are %d additional fields available, specify which one to use \r\n"
                   "for intensity or input '0' for default.                            \n\n", nf);
            scanf("%d", &iadd);
            if(iadd < 0 || iadd > nf) {
               printf("****ERROR****, additional field %d identifier not known, using default.\n\n", iadd);
               iadd = 0;
            }
            if(iadd){
               iff = 0;
               for(i=0; i < iadd; i++){
                  if(*(nfwpos + i)) iff += 3;
                  else iff += 1;
               }
               --iff;
            }

       }

       printf("Do you want to use tendency or growth rate as the intensity value              \r\n"
              "to filter on. input '0' for no, '1' for tendency or '2' for growth rate.       \n\n");
       scanf("%d", &igrwth);

       if(igrwth < 0 || igrwth > 2){
          printf("****ERROR****, intensity type %d not valid.\n\n", igrwth);
          exit(1);
       }

       if(trtyp != 'v' && igrwth == 2){
          printf("What is the time step for computing growth rates?\n\n");
          scanf("%d", &tstep);
       }

       printf("Do you want, the systems must attain value in range, input '0'\r\n"
              "             system minimum must be in range,        input '1'\r\n"
              "             system maximum must be in range,        input '2'\n\n");

       scanf("%d", &strty);

       printf("What are the lower and upper strength thresholds?\n\n");
       scanf("%f %f", &strth1, &strth2);

       printf("Check for missing location or missing value, 'y' or 'n'\n\n");
       scanf("\n");
       if(getchar() == 'y'){
         imiss = 1;
         if(iadd){
            ifloc = 0;
            if(nfwpos[iadd - 1]){
               ifloc = iff - 1;
            }
         }
       }

       if(strty){
          printf("Do you want to constrain lifecycle either side of the extremum, 'y' or 'n'\n\n");
          scanf("\n");
          if(getchar() == 'y'){
             printf("Input the length of track in time steps either side of the extremum, lneg, lpos.\n\n");
             scanf("%d %d", &lneg, &lpos);
          }

       }


       tot_fet = 0;

       tn = 0;

       for(i=0; i < tr_count; i++){

           altr = all_tr + i;

           strsw = 1;

           if(!strty){

              for(j=0; j<altr->num; j++){

                 atr = altr->trpt + j;
                 if(igrwth) {
                    if(trtyp == 'v') {
                       if(igrwth == 1) str = atr->tend;
                       else str = atr->gwthr;
                    }
                    else {
                       if(j >= altr->num-1)break;
                       atr2 = atr + 1;
                       str1 = (iadd) ? atr->add_fld[iff]: atr->zf;
                       str2 = (iadd) ? atr2->add_fld[iff]: atr2->zf;
                       if(igrwth == 1) str = str2 - str1;
                       else  str = (str2 - str1) / (float)tstep;           
                    }
                 }
                 else str = (iadd) ? atr->add_fld[iff]: atr->zf;
                 if(imiss){
                    if(str > ADD_CHECK) continue;
                    else if(iadd && ifloc){
                       if(atr->add_fld[ifloc] > ADD_CHECK) continue;
                    }
                 }
                 if((str - strth1) * (strth2 - str) > 0.0) {
                    strsw = 0;
                    break;

                 }

              }

           }

           else if (strty == 1){

              mmstr = MMSTR;

              for(j=0; j<altr->num; j++){

                 atr = altr->trpt + j;
                 if(igrwth) {
                    if(trtyp == 'v') {
                       if(igrwth == 1) str = atr->tend;
                       else str = atr->gwthr;
                    }
                    else {
                       if(j >= altr->num-1)break;
                       atr2 = atr + 1;
                       str1 = (iadd) ? atr->add_fld[iff]: atr->zf;
                       str2 = (iadd) ? atr2->add_fld[iff]: atr2->zf;
                       if(igrwth == 1) str = str2 - str1;
                       else  str = (str2 - str1) / (float)tstep;           
                    }
                 }
                 else str = (iadd) ? atr->add_fld[iff]: atr->zf;
                 if(imiss){

                    if(str > ADD_CHECK) continue;
                    else if(iadd && ifloc){
                       if(atr->add_fld[ifloc] > ADD_CHECK) continue;
                    }
                 }
                 if(str < mmstr) {mmstr = str; iextm = j;}

              }

              if((mmstr - strth1) * (strth2 - mmstr) > 0.0) strsw = 0;
           }

           else if (strty == 2){

              mmstr = -MMSTR;

              for(j=0; j<altr->num; j++){

                 atr = altr->trpt + j;
                 if(igrwth) {
                    if(trtyp == 'v') {
                       if(igrwth == 1) str = atr->tend;
                       else str = atr->gwthr;
                    }
                    else {
                       if(j >= altr->num-1)break;
                       atr2 = atr + 1;
                       str1 = (iadd) ? atr->add_fld[iff]: atr->zf;
                       str2 = (iadd) ? atr2->add_fld[iff]: atr2->zf;
                       if(igrwth == 1) str = str2 - str1;
                       else  str = (str2 - str1) / (float)tstep;           
                    }
                 }
                 else str = (iadd) ? atr->add_fld[iff]: atr->zf;
                 if(imiss){

                    if(str > ADD_CHECK) continue;
                    else if(iadd && ifloc){
                       if(atr->add_fld[ifloc] > ADD_CHECK) continue;
                    }
                 }
                 if(str > mmstr) {mmstr = str; ; iextm = j;}

              }
 
              if((mmstr - strth1) * (strth2 - mmstr) > 0.0) strsw = 0;

           }


           if(strty && !strsw){
              strsw = 1;
              ld1 = (altr->trpt + iextm)->fr_id - altr->trpt->fr_id;
              ld2 = (altr->trpt + altr->num - 1)->fr_id - (altr->trpt + iextm)->fr_id;
              if(ld1 >= lneg && ld2 >= lpos) strsw = 0;

           }

           if(strsw){

	      for(j=0; j < altr->num; j++) {if(altr->trpt) free((altr->trpt + j)->add_fld);}
              altr->num = 0;
              free(altr->trpt);
              altr->trpt = NULL;

           }
           else {
              ++tn;
              tot_fet += altr->num;
           }

       }    


       printf("***INFORMATION***, current number of tracks is %d\n\n", tn);   

    }

    printf("Do you want to filter according to propogation direction,\r\n"
           "and/or speed, 'y' or 'n'.                                    \n\n");
    scanf("\n");
    if(getchar() == 'y') {
       printf("Use propogation direction filter, 'y' or 'n'.\n\n");
       scanf("\n");
       if(getchar() == 'y') prop_filt(all_tr, tr_count, &tot_fet);
       printf("Use propogation speed filter, 'y' or 'n'.\n\n");
       scanf("\n");
       if(getchar() == 'y') speed_filt(all_tr, tr_count, &tot_fet);      
    }

    printf("Do you want to restrict the frame range, 'y' or 'n'\r\n"
           "This will identify tracks that overlap with the    \r\n"
           "frame range by at least %5.2f%%                    \n\n", 100*LPERC );
    scanf("\n");
    if(getchar() == 'y'){

       tot_fet = 0;
       tn = 0;

       printf("Input the frame range required, f1, f2\n\n");
       scanf("%d %d", &ifr1, &ifr2);

       if(ifr2 <= ifr1){
          printf("****ERROR****, incorrect range chosen.\n\n");
          exit(1);

       }

       for(i=0; i < tr_count; i++){

         altr = all_tr + i;
         atr = altr->trpt;

         if(!altr->num) continue;

         nmiss = 0;
         for(j=0; j<altr->num; j++) nmiss += (atr+j)->nfm;

         ifb =  atr->fr_id;
         ife =  (atr + (altr->num - 1))->fr_id;

         id1 = id2 = idd = iout = 0;

         if(ife < ifr1 && ifb < ifr1) iout = 1;
         else if(ife > ifr2 && ifb > ifr2) iout = 1;
         else {
            idd = ife - ifb + nmiss;
            if(ifb < ifr1 && ife >= ifr1) id1 = ifr1 - ifb;
            if(ife > ifr2 && ifb <= ifr2) id2 = ife - ifr2;


            if((float)(id2 + id1)/(float)idd > LPERC) iout = 1;


         }

         if(iout){
            if(altr->trpt){
	       for(j=0; j < altr->num; j++) free((altr->trpt + j)->add_fld);
               altr->num = 0;
               free(altr->trpt);
               altr->trpt = NULL;
            }

         }

         else {tot_fet += altr->num;  ++tn;}

       }

       printf("***INFORMATION***, current number of tracks is %d\n\n", tn);

    }
    
/* select individual tracks or use all tracks. */
   
    printf("Do you want to use all tracks or a selected track, 'a' for all and 's' for selected.\n\n");
    scanf("\n");
    if(getchar() == 's'){
       printf("Use track number or track Id. to select track, input '0' for track number and '1' for track Id.\n\n");
       scanf("%d", &trty);
       if(trty < 0 || trty > 1){
          printf("****ERROR****, incorrect specifier for track identification type.\n\n");
       	  exit(1);
       }
       if(!trty){
          printf("What track number is required?\n\n");
	  scanf("%d", &tnum);
	  if(tnum < 1 || tnum > tr_count){
	     printf("****ERROR****, track number %d not available.\n\n", tnum);
	     exit(1);
	  }

          for(i=0; i < tr_count; i++){
	     altr = all_tr + i;
	     if(!(i+1 == tnum)){
	        for(j=0; j < altr->num; j++) free((altr->trpt + j)->add_fld);
                free(altr->trpt);
	        altr->num = 0;
	     } 
	    
	  }	  
	  
       }
       else {
          printf("What track Id. is required?\n\n");
	  scanf("%d", &tnum);
	 
	  for(i=0; i < tr_count; i++){
	     altr = all_tr + i;
	     if(!(altr->trid == tnum)){
	        for(j=0; j < altr->num; j++) free((altr->trpt + j)->add_fld);
                free(altr->trpt);
	        altr->num = 0;
	     } 
	    
	  }
       }
    }
    else {
       printf("****INFORMATION****, using all available tracks.\n\n");
    }

/* track deviations */

    if(IDEV) track_sample_dev(all_tr, tr_count);

/* write filtered data to file */

    strncpy(trout, FILTRS, MAXCHR);
    if(iext) strcpy(strstr(trout, EXTENSION), fext);

    tsf = open_file(trout, "w");

    meantrd(tsf, all_tr, tr_count, trtyp);  
 
    close_file(tsf, trout);

    write_track_netcdf(all_tr, tr_count, trout, trtyp); 

/* free track data and re-read to get correct number of tracks */

    for(i=0; i < tr_count; i++) {
       altr = all_tr + i;
       if(altr->num){
          for(j=0; j < altr->num; j++) {if(altr->trpt) free((altr->trpt + j)->add_fld);}
          free(altr->trpt);  
       }    
    }

    free(all_tr);

    tsf = open_file(trout, "r");

    all_tr = read_tracks(tsf, &tr_count, &gpr, &tpr, trtyp, &gr->alat, &gr->alng);

    close_file(tsf, trout);

    tot_fet = 0;
    for(i=0; i < tr_count; i++) tot_fet += (all_tr+i)->num;

    printf("Do you want nearest grid point positions, 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y'){

       ntrnum=0;

       for(i=0; i < tr_count; i++){

           altr = all_tr + i;
           atr = altr->trpt;

           for(j=0; j < altr->num; j++){

               at = atr + j;
               at->ix = 0;
               at->iy = 0;
               dd = BDIST;

               for(k=0; k < gr->ix; k++) {

                   if((ddf=fabs(*(gr->xgrid + k) - at->xf)) < dd){dd = ddf; at->ix = k+1;}

               }

               dd = BDIST;

               for(k=0; k < gr->iy; k++) {

                   if((ddf=fabs(*(gr->ygrid + k) - at->yf)) < dd){dd = ddf; at->iy = k+1;}

               }

           }


       }

       strncpy(grtrs, GRTRS, MAXCHR);
       if(iext) strcpy(strstr(grtrs, EXTENSION), fext);

       tsf = open_file(grtrs, "w");

       pl = ftello(tsf);

       fprintf(tsf, "TRACK_NUM  %8d\n", tr_count);

       for(i=0; i < tr_count; i++){

           altr = all_tr + i;

           if(altr->num){

              ++ntrnum;

              fprintf(tsf, "TRACK_ID  %d\n", i+1);
              fprintf(tsf, "POINT_NUM  %d\n", altr->num);

              for(j=0; j < altr->num; j++){

                 atr = altr->trpt + j;
                 fprintf(tsf, "%d %d %d\n", atr->fr_id, atr->ix, atr->iy);

              }

           }

       }

       fseeko(tsf, pl, FSTART);
       fprintf(tsf, "TRACK_NUM  %8d", ntrnum);
 
       close_file(tsf, grtrs);


    }

/*    printf("Do you want to compute principle tracks? 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y') fourier_match(all_tr, tr_count, &tot_fet); */

/* produce a plot of tracks if required */

replot:

    if(!(gr->prgr) && !(gr->prty)){
       printf("do you want the combined track data plotted, y or n\n\n");

       scanf("\n");
       tplt = getchar();

       if(tplt == 'y'){

          printf("if a plot is required of the combined track data\r\n"
                 "do you want area of interest (input 1) or whole area (input 0)\n\n");
          scanf("%d", &izm);

          splice_plot(all_tr, tot_fet, tr_count, izm, 0, trtyp);


       }

       printf("do you want to change the minimum lifetime and re-plot\r\n"
              "or analyse a different combined data set\r\n"
              "or change the area of interest, 'y' or 'n'\n\n");

       scanf("\n");
       if(getchar() == 'y') {

          printf("***INFORMATION***, to change the minimum lifetime data must be reread\n\n");
          printf("do you want to change the minimum lifetime or choose a different data set, 'y' or 'n'\n\n");

          scanf("\n");
          if(getchar() == 'y'){ ex = 'y'; goto reread; }

          printf("do you wish to change the area of interest and replot, 'y' or 'n'\n\n");

          scanf("\n");
          if(getchar() == 'y'){

             printf("input new area of interest\n\n");

             printf("X grid numbers =\n");
             scanf("%d %d", &x1u, &x2u);

             printf("Y grid numbers =\n");
             scanf("%d %d", &y1u, &y2u);

             xmn = *(gr->xgrid + x1u - 1);
             xmx = *(gr->xgrid + x2u - 1);
             ymn = *(gr->ygrid + y1u - 1);
             ymx = *(gr->ygrid + y2u - 1);

             goto replot;


          }

       }
    
    }

    *trackn = tr_count;

    return all_tr;

}
