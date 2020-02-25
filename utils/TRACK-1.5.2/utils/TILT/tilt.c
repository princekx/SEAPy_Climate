#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"
#include "m_values.h"
#include "vec.h"

#define  NBIN     201
#define  PLARGE   1.0e+10
#define  NLARGE  -1.0e+10


/* combine sets of tr_trs files */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);
double fvec(VEC * , VEC * , VEC * , VEC * );
void centre( struct tot_tr * , int , int , int );
void sincos(double , double * , double * );

extern int aniso;
extern int nfld, nff;
extern int *nfwpos;


int main(void)

{
    int i, j, k;
    int iu, iv;
    int idirv=0, idirt=0;
    int iuf=0, ivf=0;
    int iuvpx=0, iuvpy=0;
    int tr_count=0;
    int gpr, ipr;
    int nlev=0;
    int *ilev=NULL;
    int *ipx=NULL;
    int wnc='n';
    int iavg=0, iref=0;
    int ifint=0, ifft=0;
    int lfnd=0;

    int icent=NBIN/2;
    int **ntlt=NULL;
    int nloc[NBIN];

    char fncout[100], stub[100];

    float alat, alng;
    float xx, yy;
    float uu, vv;
    float norm;
    float arcl;
    float loct=0.0;

    float **tavg=NULL;
    float *lev=NULL;
    float loc[NBIN];

    double sn1, cn1, sn2, cn2;
    double aa, bb, dotab;

    char filnamin[100];

    FILE *fin=NULL;
    FILE *fout=NULL;
    FILE *fnc=NULL;
    FILE *favg=NULL;

    struct tot_tr *all_tr=NULL, *altr=NULL;
    struct fet_pt_tr *atr=NULL;

    VEC rad, tang, rx, ry;
    VEC levv, nwv;

    printf("What is the track file to read?\n\n");
    scanf("%s", filnamin);

    fin=fopen(filnamin, "r");
    if(!fin){
       printf("****ERROR****, can't open file %s\n", filnamin);
       exit(1);
    }

    all_tr = read_tracks(fin, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);

    printf("This program computes storm tilts relative to the flow direction,\r\n"
           "or the storm propagation direction.                              \r\n"
           "For this it requires field values at different levels with storm \r\n"
           "center positions. It also requires the flow direction provided by\r\n"
           "large scale wind field if tilt with respect to flow is required. \n\n");

    printf("Do you want the tilt with respect to the flow or propagation? \r\n"
           "Input '0' for flow or '1' for propagation.                    \n\n");
    scanf("%d", &idirt);

    if(!idirt){

       printf("There are %d additional fields available, which ones are the U and V\r\n"
              "components of the large scale flow?                                 \n\n", nff);
       scanf("%d %d", &iu, &iv);
       if((iu < 1 || iu > nff) || (iv < 1 || iv > nff)){
          printf("****ERROR****, incorrect Id's for U or V.\n\n");
          exit(1);
       }

       printf("What field ID should be used for locating the directional\r\n"
              "vector? Input '0' for default.                           \n\n");
       scanf("%d", &idirv);
       if(! *(nfwpos + idirv - 1)){
         printf("****ERROR****, additional field Id., has no positional information.\n\n");
         exit(1);
       }

    }
    else {
       printf("****INFORMATION****, note that the directional vector will be computed from the default positions.\n\n");
    }

    printf("How many other levels are to be used?\n\n");
    scanf("%d", &nlev);
    ilev = (int *)calloc(nlev, sizeof(int));
    mem_er((ilev == NULL) ? 0 : 1, nlev*sizeof(int));

    ipx = (int *)calloc(nlev, sizeof(int));
    mem_er((ipx == NULL) ? 0 : 1, nlev*sizeof(int));

    lev = (float *)calloc(nlev, sizeof(float));
    mem_er((lev == NULL) ? 0 : 1, nlev*sizeof(float));

    printf("What are the other field Id's required and their level values.\n\n");
    for(i=0; i < nlev; i++){
       printf("Input next Id and level ");
       scanf("%d %f", ilev + i, lev + i);
       if(*(ilev + i) < 0 || *(ilev + i) > nff){
         printf("****ERROR****, additional field Id %d greater than available number of fields %d,\n\n", *(ilev+i), nff);
         exit(1);
       }
       if(!(*(nfwpos + *(ilev + i) - 1))){
          printf("****ERROR****, this additional field Id has no positional information.\n\n");
          exit(1);
       }
       for(j=0; j < *(ilev + i); j++){
          if(*(nfwpos + j)) *(ipx + i) += 3;
          else *(ipx + i) += 1;
       }
       --(*(ipx + i));

    }
    printf("\n\n");

    printf("Write individual tracks as seperate netcdf files, 'y' or 'n'\n\n");
    scanf("\n");
    wnc = getchar();

    printf("Do you want to average the tilts, centered on the maximum or minimum intensity? 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y'){
       iavg = 1;
       printf("Center tracks wrt '0' max. intensity     \r\n"
              "                  '1' min. intensity     \n\n");
       scanf("%d", &iref);
       printf("Which of the additional fields do you want to use for intensity, input '0' for default?\n\n");
       scanf("%d", &ifint);
       if(ifint < 0 || ifint > nff){
          printf("****ERROR****, not a valid field Id.\n\n");
          exit(1);
       }
       ifft = 0;
       if(ifint){
          for(i=0; i < ifint; i++){
              if(*(nfwpos + i)) ifft += 3;
              else ifft += 1;
          }
          --ifft;
       }

       tavg = (float **)calloc(nlev, sizeof(float *));
       mem_er((tavg == NULL) ? 0 : 1, nlev*sizeof(float *));
       ntlt = (int **)calloc(nlev, sizeof(int *));
       mem_er((ntlt == NULL) ? 0 : 1, nlev*sizeof(int *));
       for(i=0; i < nlev; i++){
           *(tavg + i) = (float *)calloc(NBIN, sizeof(float));
           mem_er((*(tavg + i) == NULL) ? 0 : 1, NBIN*sizeof(float));
           *(ntlt + i) = (int *)calloc(NBIN, sizeof(int));
           mem_er((*(ntlt + i) == NULL) ? 0 : 1, NBIN*sizeof(int));
           for(j=0; j < NBIN; j++){*(*(tavg + i) + j) = 0.0; *(*(ntlt + i) + j) = 0; loc[j] = 0.0; nloc[j] = 0;}
       }

    }

    if(!idirt){
    
       for(i=0; i < iu; i++){
           if(*(nfwpos + i)) iuf += 3;
           else iuf += 1;
       }

       --iuf;

       for(i=0; i < iv; i++){
           if(*(nfwpos + i)) ivf += 3;
           else ivf += 1;
       }

       --ivf;

       if(idirv){
          for(i=0; i < idirv; i++){
              if(*(nfwpos + i)) iuvpx += 3;
              else iuvpx += 1;
          }

          iuvpx -= 3;

          iuvpy = iuvpx + 1;

       }

    }

    fout = fopen("tilt.dat", "w");
    if(!fout){
       printf("****ERROR****, can't open file %s\n", "tilt.dat");
       exit(1);
    }

    fprintf(fout, "NTRACK %d NFIELD %d\n", tr_count, nlev);
    fprintf(fout, "LEVELS ");
    for(i=0; i < nlev; i++) fprintf(fout, "%f ", *(lev + i));
    fprintf(fout, "\n");    

    for(i=0; i < tr_count; i++){

       altr = all_tr + i;

       printf("Processing track %d\n", i+1);

       if(iavg) centre(altr, iref, ifint, ifft);

       fprintf(fout, "TRACK_NO %d NUMPT %d\n", i+1, altr->num);

       if(wnc == 'y'){

          sprintf(fncout, "tr%04d.nc", i+1);
          sprintf(stub, "tr%04d", i+1);
          
          fnc = fopen(fncout, "w");
          if(!fnc){
             printf("****ERROR****, can't open file %s\n", fncout);
             exit(1);
          }
          fprintf(fnc, "netcdf %s{\n", stub);
          fprintf(fnc, "dimensions:\n\t time\t = %d;\n\t level\t = %d;\n", altr->num, nlev);
          fprintf(fnc, "variables:\n\t float time(time);\n\t float level(level);\n\t float tilt(time, level);\n\t       tilt:missing_value = %e;\n", ADD_UNDEF);
          fprintf(fnc, "data:\n\t time = ");
          if(altr->time){
             for(j=0; j < altr->num; j++){
                 atr = (altr->trpt) + j;
                 fprintf(fnc, "%ld, ", atr->time);
             }
          }
          else {
             for(j=0; j < altr->num; j++){
                 atr = (altr->trpt) + j;
                 fprintf(fnc, "%d, ", atr->fr_id);
             }
          }
          fseek(fnc, -2, SEEK_CUR);
          fprintf(fnc, ";\n");
          fprintf(fnc, "\t level = ");
          for(j=0; j < nlev; j++) fprintf(fnc, "%5.1f, ", *(lev + j));
          fseek(fnc, -2, SEEK_CUR);
          fprintf(fnc, ";\n");           
          fprintf(fnc, "\t tilt = ");

       }

       for(j=0; j < altr->num; j++){

          atr = (altr->trpt) + j;

          if(altr->time) fprintf(fout, "%ld ", atr->time);
          else fprintf(fout, "%d ", atr->fr_id);

/* compute radial and tangential vectors */


          if(!idirt){

             if(!idirv){
                xx = atr->xf;
                yy = atr->yf;
             }
             else {
                xx = *(atr->add_fld + iuvpx);
                yy = *(atr->add_fld + iuvpy);
             }

             xx *= FP_PI;
             yy = FP_PI2 - yy * FP_PI;
             if(yy < 0.) yy = 0.;
             sincos(xx, &sn1, &cn1);
             sincos(yy, &sn2, &cn2);
             rad.x = sn2 * cn1;
             rad.y = sn2 * sn1;
             rad.z = cn2;

             rx.x = -sn1;
             rx.y = cn1;
             rx.z = 0.0;

             ry.x = cn2 * cn1;
             ry.y = cn2 * sn1;
             ry.z = -sn2;

             uu = *(atr->add_fld + iuf);
             vv = *(atr->add_fld + ivf);

             mulv(&rx, &rx, uu);
             mulv(&ry, &ry, vv);

             addv(&rx, &ry, &tang);
             norm = sqrt(dotp(&tang, &tang));
             normv(&tang, norm);

          }

          else {
             if(!j){
                xx = atr->xf;
                yy = atr->yf;
                xx *= FP_PI;
                yy = FP_PI2 - yy * FP_PI;
                if(yy < 0.) yy = 0.;
                sincos(xx, &sn1, &cn1);
                sincos(yy, &sn2, &cn2);
                rad.x = sn2 * cn1;
                rad.y = sn2 * sn1;
                rad.z = cn2;
             }

             else {
                rad = rx;
             }

             if(j < altr->num - 1){
                xx = (atr+1)->xf;
                yy = (atr+1)->yf;
                xx *= FP_PI;
                yy = FP_PI2 - yy * FP_PI;
                if(yy < 0.) yy = 0.;
                sincos(xx, &sn1, &cn1);
                sincos(yy, &sn2, &cn2);
                rx.x = sn2 * cn1;
                rx.y = sn2 * sn1;
                rx.z = cn2;

                crosp(&rad, &rx, &ry);
                crosp(&ry, &rad, &tang);
                norm = sqrt(dotp(&tang, &tang));
                normv(&tang, norm);

             }

          }

/*          printf("%f %f %f\n", rad.x, rad.y, rad.z);
          printf("%f %f %f\n", tang.x, tang.y, tang.z);
          printf("%f\n", dotp(&tang, &rad)); */

          lfnd = 0;

          for(k=0; k < nlev; k++){
              xx = *(atr->add_fld + *(ipx + k) - 2);
              yy = *(atr->add_fld + *(ipx + k) - 1);

              xx *= FP_PI;
              yy = FP_PI2 - yy * FP_PI;
              if(yy < 0.) yy = 0.;
              sincos(xx, &sn1, &cn1);
              sincos(yy, &sn2, &cn2);
              levv.x = sn2 * cn1;
              levv.y = sn2 * sn1;
              levv.z = cn2;

              dotab = fvec(&rad, &tang, &levv, &nwv);

              aa = dotp(&rad, &nwv);
              bb = dotp(&tang, &nwv);

              arcl = acos(aa) / FP_PI;
              if(dotab <= ((aa < bb) ? aa: bb)) arcl = arcl;
              else arcl = -arcl;

              if(xx > ADD_CHECK) arcl = ADD_UNDEF;
              else if(iavg){
                 if(!lfnd){
                    loct = arcl;
                    loc[icent + atr->fr_id] += arcl;
                    nloc[icent + atr->fr_id] += 1;
                    lfnd = 1;
                 }
                 tavg[k][icent + atr->fr_id] += arcl - loct;
                 ntlt[k][icent + atr->fr_id] += 1;

              }
              fprintf(fout, "%e ", arcl);
              if(wnc == 'y') fprintf(fnc, "%e, ", arcl);

          }

          fprintf(fout, "\n");


       }


       if(wnc == 'y'){
          fseek(fnc, -2, SEEK_CUR);
          fprintf(fnc, ";}\n");
          fclose(fnc);
          fnc = NULL;
       }

    }

    if(iavg){
       for(i=0; i < NBIN; i++) {
           if(nloc[i]) loc[i] /= nloc[i];
       }
       for(i=0; i < nlev; i++){
           for(j=0; j < NBIN; j++){
/*printf("%f %d %f\n", tavg[i][j], ntlt[i][j], loc[i]);*/
               if(ntlt[i][j]) {
                  tavg[i][j] /= ntlt[i][j];
                  tavg[i][j] += loc[i];
               }
 
               else tavg[i][j] = ADD_UNDEF;
/*printf("%f\n", tavg[i][j]);*/
           }
       }
    }

    favg = fopen("tilt_avg.nc", "w");
    if(!favg){
       printf("****ERROR****, can't open file %s\n", fncout);
       exit(1);
    }

    fprintf(favg, "netcdf tilt_avg{\n");
    fprintf(favg, "dimensions:\n\t time\t = %d;\n\t level\t = %d;\n", NBIN, nlev);
    fprintf(favg, "variables:\n\t float time(time);\n\t float level(level);\n\t float tilt(time, level);\n\t       tilt:missing_value = %e;\n", ADD_UNDEF);
    fprintf(favg, "data:\n\t time = ");
    for(i=0; i < NBIN; i++){
       fprintf(favg, "%d, ", -icent + i);
    }
    fseek(favg, -2, SEEK_CUR);
    fprintf(favg, ";\n");
    fprintf(favg, "\t level = ");
    for(i=0; i < nlev; i++) fprintf(favg, "%5.1f, ", *(lev + i));
    fseek(favg, -2, SEEK_CUR);
    fprintf(favg, ";\n");           
    fprintf(favg, "\t tilt = ");
    for(i=0; i < NBIN; i++){
        for(j=0; j < nlev; j++){
            fprintf(favg, "%e, ", tavg[j][i]);
        }
    }
    fseek(favg, -2, SEEK_CUR);
    fprintf(favg, ";}\n");
    fclose(favg);


    if(tavg){
       for(i=0; i < nlev; i++) {free(*(tavg + i)); free(*(ntlt + i));}
       free(tavg);
       free(ntlt);
    }

    fclose(fout);

    return 0;

}


double fvec(VEC *a, VEC *b, VEC *c, VEC *r)
{

    VEC an, bn;

    double aa, bb, dotab;
    double norm=0.0;

    dotab = dotp(a, b);

    aa = dotp(a, c) - dotp(b, c) * dotab;
    bb = dotp(b, c) - dotp(a, c) * dotab;

    mulv(a, &an, aa);
    mulv(b, &bn, bb);

    addv(&an, &bn, r);

    norm = sqrt(dotp(r, r));
    normv(r, norm);

    return dotab;

}


void centre(struct tot_tr *trr, int iref, int ifr, int iffr)
{
    int i;
    int ipref=0;
    int frid=0;

    float ref, str;

    struct fet_pt_tr *fp=NULL;

    if(!iref) ref = NLARGE;
    else ref = PLARGE;


/* find reference point */
    
    for(i=0; i < trr->num; i++){
       fp = trr->trpt + i;
       str = (ifr) ? fp->add_fld[iffr] : fp->zf;

       if(!iref){
          if(str > ref) {ref = str; ipref = i;}
       }
       else{
          if(str < ref) {ref = str; ipref = i;}
       }

    }

    frid = (trr->trpt + ipref)->fr_id;

    for(i=0; i < trr->num; i++){
       fp = trr->trpt + i;
       fp->fr_id -= frid;
    }


    return;
}
