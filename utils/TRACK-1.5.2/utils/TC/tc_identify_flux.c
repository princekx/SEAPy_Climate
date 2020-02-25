#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"
#include "vec.h"
#include "geo_values.h"
#include "st_obj.h"
#include "st_fo.h"

#define  SMALL -1.0e+12
#define  LARGE 1.0e+12
#define  GTOL  1.0e-4

/* program to identify tropical cyclones and produce composites, now extended to extra-tropical cyclones 
   and energetics                                                                                        */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
void write_fields(char * , double ** , float * , float * , int , int , int );
void convert_track(struct tot_tr * , int );
double anomaly(float * , double * , int , int , int );
void sincos(double , double * , double * );
float *write_fields_netcdf(char * , double ** , float * , float * , float * , int , int , int , int , int , int , int , long int );
void write_intflux_netcdf(char * , double ** , float * , float * , int , int , int , int , int , long int );
void write_profile_netcdf(char * , double * , float * , int , int , int , long int );

void radial_extreme(FILE * , FILE * , FILE * , float * , float * , float, int , int , off_t * , off_t * , off_t * );
double dirangle(struct tot_tr * , double * , double * , int , int );
float *leveldata(int );
double prop_speed(struct tot_tr * , int );
void line_area_integ(double * , double * , double * , double , double , double * , double , int , int);
FILE *open_radial_file(char * , int * , long int * , int * , int * , int * , int * , int * , int * , int * , int * );
double line_anomaly(float * , int , int );
double vert_integ(double * , double * , int , int , int );
void setup_rotmatrix(double , double , double , double , VEC * );
void setup_radvec(VEC * , VEC * , VEC * , VEC * , VEC * , VEC * , VEC * , VEC * , double , double , int , int );
void find_last(struct fet_pt_tr * , struct feature_pts * , int * , int * , int , int );
double fvec(VEC * , VEC * , VEC * , VEC * );
void read_model_levs(char * , float ** , float ** , int * );
void model_level_diff(float * , float * , float * , float * , int , int );

int orog_test(float * , float * , int , int , int * , float , float );

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int main(void )
{

    int i, j, k;
    int trnum=0;
    int gpr=0, ipr=0;
    int istc=0;
    int nth=1;
    int ift=0, if1=0, if2=0, ifs=0;
    int iff1=0, iff2=0, ifft=0;
    int iv=0, nlp=0, nl1=0, nl2=0;
    int nn=0;
    int ntrack=0;
    int nvr=0;
    int ireg=0;
    int inreg=0;
    int ilm=0;
    int lglng, lglat;
    int ilng=0, ilat=0, igen=0;
    int ianom=0, iwind=0, irel=0;
    int itv='n';
    int igg=0;
    int nvert=0;
    int imiss=0;
    int istd=0;
    int iflx=0, ifanom=0, iffl=0, irif=0;
    int irr=0, iranom=0, iflxrd=0;
    int iaddf=0, ipnorm=0;
    int ndsmth=0, ifcnt=0, itpadd=0, iprojd=0;
    int iadv=0, idevazm=0;
    int toff=0;
    int iselt=0, itrst=0, itren=0;
    int *nalv=NULL;
    int yinteg=0;
    int ikecb=0, ike_cb_anom=0;
    int lev_typ=0, lyes=0, nlab=0;

    int irog=0;

    int irad=0, irdim=0;
    int irtrn=0, irnr=0, irnth=0, irnf=0, irnft=0;
    int iext=0, imxmn=1, ifld=0;
    int nll=0, nll1=0, nll2=0;

    int frst=0, last=0;
    int ispd=0, ifsp=0;
    int isprec=0, ifsprec=0;
    int nwtr=0, nctr=0;
    int icomp_order=0;

    int *ilms=NULL;
    int idir=0;

    int *ifncnt=NULL, *ifncntp=NULL;

    long int iptnum;
    long int tim1, tim2, tim=0;
    long int ctime=0;
    off_t place1=0, place2=0, blklen=0;
    off_t **rmposx=NULL, **rmposn=NULL;
    off_t **rrpos=NULL;

    float alat=0.0, alng=0.0;
    float thresh=0.0, fth=0.0;
    float str=0.0;
    float tend=0.0, str2=0.0;
    float rthresh=0.0;

    float *lms=NULL;
    float *glng=NULL, *glat=NULL;

    float lng1, lng2, lat1, lat2;
    float fmax, fmin, lngm=0.0, latm=0.0;
    float fff=0.0, fdif=0.0, fmxmn=0.0, fdifmax=0.0;
    float sumprc=0.0, sumpp=0.0;

    float glat1, glat2;
    float **ncount=NULL;

    float *sradf1=NULL, *sradf2=NULL;
    float *sadd=NULL, *sadd2=NULL;;
    float *tmp1=NULL, *tmp2=NULL;
    float *sgrdf1=NULL, *sgrdf2=NULL;
    float *slng=NULL, *slat=NULL;
    float *slng2=NULL, *slat2=NULL;
    float rad=0.0, thet=0.0, sq=0.0, rr=0.0;
    float *levdat=NULL;
    float ke=0.0, nrad=0.0, anom_rad=0.0;

    float *sspres=NULL;

    float *tmpdat=NULL;

    float *aa=NULL, *bb=NULL;
    float ad=0.0, bd=0.0;

    double pdi=0.0;
    double s1=0.0, s2=0.0, c1=0.0, c2=0.0;
    double ss1=0.0, cc1=0.0;
    double norm=0.0, arot=0.0;
    double pdirl=0.0, tstp=0.0;
    double ang=0.0;
    double lscl=0.0, ascl=0.0, asclr=0.0;
    double scladv=0.0;
    double xx, yy;

    double uu=0.0, vv=0.0;
    double gx=0.0, gy=0.0;
    double vphi=0.0;

    double lavg1=0.0, lavg2=0.0, aravg=0.0;

    double sumarc=0.0;
    double arcll=0.0, arclt=0.0;
    double vlintg1=0.0, vaintg1=0.0, vlintg2=0.0, vaintg2=0.0;
    double vkeintg=0.0;

    double *coslat=NULL;
    double *rintg1=NULL, *rintg2=NULL, *rintg3=NULL;
    double *vint=NULL, *vke=NULL;
    double **linteg1=NULL, **linteg2=NULL, **linteg3=NULL; 
    double *lvtmp=NULL;
    double *alv=NULL;

    double *vsavg1=NULL, *vsavg2=NULL;

    double **aavg=NULL, **aavg1=NULL, **aavg2=NULL;
    double **savg1=NULL, **savg2=NULL, **savg3=NULL;
    double **astd=NULL, **astd1=NULL, **astd2=NULL;

    VEC ptt={0.0, 0.0, 0.0}, ww={0.0, 0.0, 0.0};
    VEC erd={0.0, 0.0, 0.0}, etn={0.0, 0.0, 0.0};
    VEC vel={0.0, 0.0, 0.0}, vx={0.0, 0.0, 0.0}, vy={0.0, 0.0, 0.0};

    VEC pt1,  pt2, vtmp;
    VEC vt1, vt2, vt3, vn;

    VEC *vecg=NULL, *vt=NULL;
    VEC rot[3]={{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    VEC vsdir={0.0, 0.0, 0.0}, tvec={0.0, 0.0, 0.0};
    VEC *newv=NULL, *etht=NULL, *ephi=NULL, *erad=NULL, *etan=NULL;
    VEC *erdd=NULL, *etnn=NULL;

    FILE *fin=NULL, *fout=NULL, *flm=NULL;
    FILE *freg=NULL, *fprc=NULL;
    FILE *frad1=NULL, *frad2=NULL, *fadd=NULL;
    FILE *ftim=NULL;
    FILE *fdout=NULL;
    FILE *frmax=NULL, *frmin=NULL, *frange=NULL;
    FILE *fthresh=NULL;
    FILE *fgrdx=NULL, *fgrdy=NULL;
    FILE *fdir=NULL;
    FILE *flev=NULL;
    FILE *fsfp=NULL;

    char filin[MAXCHR];
    char filout[MAXCHR];
    char filrad[MAXCHR];
    char filsfp[MAXCHR];
    char buff[MAXCHR];
    char flevnam[MAXCHR];
    char hurr[]="track.dat";

    struct tot_tr *tracks=NULL, *tracks2=NULL, *atr=NULL, *atr2=NULL;
    struct fet_pt_tr *fp=NULL, *fpc=NULL;
    struct feature_pts fpt;

    freg = fopen("region.dat", "r");
    if(freg){
       ireg = 1;
       fscanf(freg, "%f %f", &lng1, &lng2);
       fscanf(freg, "%f %f", &lat1, &lat2);
       fscanf(freg, "%ld %ld", &tim1, &tim2);
       printf("****INFORMATION****, using region information from region.dat file   \r\n"
              "                     region is (%f, %f) longitude, (%f, %f) latitude,\r\n"
              "                     (%ld, %ld) time.                                \n\n",
              lng1, lng2, lat1, lat2, tim1, tim2);
       fclose(freg);
    }
    else {
       printf("****INFORMATION****, no region information available.\n\n");
    }


    flm = fopen("lmask.dat", "r");
    if(flm){
       printf("****INFORMATION****, using land mask file lmask.dat\n\n");
       ilm = 1;
       fgets(buff, MAXCHR, flm);
       sscanf(buff, "%d %d", &lglng, &lglat);
       glng = (float *)calloc(lglng, sizeof(float));
       glat = (float *)calloc(lglat, sizeof(float));
       lms = (float *)calloc(lglng*lglat, sizeof(float));
       ilms = (int *)calloc(lglng*lglat, sizeof(int));
       for(i=0; i < lglng; i++) fscanf(flm, "%f", glng + i); 
       for(i=0; i < lglat; i++) fscanf(flm, "%f", glat + i);
       fgets(buff, MAXCHR, flm);
       fgets(buff, MAXCHR, flm);
       fread(lms, lglng*lglat*sizeof(float), 1, flm);
       fclose(flm);
       for(i=0; i < lglng*lglat; i++) *(ilms + i) = (int)(*(lms + i)); 

    }

    else {
       printf("Do you want genesis region filtering, used mainly for TC's.\n\n");
       scanf("\n");
       if(getchar() == 'y'){
          igg = 1;
          printf("What are the latitudinal range for the genesis filter, lat1 and lat2?\n\n");
          scanf("%f %f", &glat1, &glat2);
          if(glat2 < glat1){
             printf("****ERROR****, latitudinal range incorrect, lat1 > lat2, should be lat1 < lat2.\n\n");
             exit(1);
          }
       }
    }

    printf("What is the track file to read?\n\n");
    scanf("%s", filin);

    fin = fopen(filin, "r");
    if(!fin){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
       exit(1);
    }

    tracks = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    itrst = 0;
    itren = trnum;

    fseeko(fin, 0L, SEEK_SET);

    tracks2 = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);

    printf("Are there along track radial field files available, 'y' or 'n'\n\n");
    scanf("\n");
    if(getchar() == 'y'){
       irad = 1;

       printf("Are fields on pressure, '0' or model '1' levels, this is only used for the vertical integrals of the energetics calculations.\n\n");
       scanf("%d", &lev_typ);
       if(lev_typ < 0 || lev_typ > 1){
          printf("****ERROR****, level type %d not a valid identifier.\n\n", lev_typ);
          exit(1);
       }

       if(lev_typ){
          printf("****WARNING****, order of integration will be changed so that vertical integration is performed first.\n\n");
       }

       printf("Do you want wind composites? 'y' or 'n'\n\n");
       scanf("\n");
       if(getchar() == 'y') iwind = 1;

       printf("What is the filename for the first radial field, U for winds?\n\n");
       scanf("%s", filrad);

       frad1 = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);

       if(irtrn != trnum){
          printf("****ERROR****, track file and radial field data not compatable: different numbers of tracks.\n\n");
          exit(1);
       }

       if(iwind){
          printf("What is the filename for the second radial field, V for winds?\n\n");
          scanf("%s", filrad);
          frad2 = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);

          if(idir) {
             if(!ndsmth){
                printf("Do you want to use the current point to determine system direction or average over several points\r\n"
                       "to improve direction smoothness, input '1' for single point or 'n' the number of points to use.  \r\n"
                       "This will depend on the track lifetime and minimum lifetime threshold.                           \n\n");
                scanf("%d", &ndsmth);
                if(ndsmth < 1){
                   printf("****ERROR****, number of points cannot be zero or negative for direction smoothing, exiting.\n\n");
                   exit(1);
                }

             }

             icomp_order = 1;
             printf("Do you want winds relative to direction velocity? 'y' or 'n'\n\n");
             scanf("\n");
             if(getchar() == 'y') {
                irel = 1;
                printf("What is the time step in seconds for speeds in m/s?\n\n");
                scanf("%lf", &tstp);
             }


             
          }
          else {

             printf("Do you want to composite U and V first then determine tangential and radial components, \r\n"
                    "or determine tangential and radial components then composite. Input '0' for composite   \r\n"
                    "first or '1' for components first.                                                      \n\n");
             scanf("%d", &icomp_order);
             if(icomp_order < 0 || icomp_order > 1){
                printf("****WARNING****, %d is wrong value for option, defaulting to '0'\n\n", icomp_order);
                icomp_order = 0;
             }

          }

          convert_track(tracks, trnum);

          if(icomp_order == 1){

/* check if tilts are considered */

             if(ifcnt){
                ifncnt = (int *)calloc(irnf, sizeof(int));
                mem_er((ifncnt == NULL) ? 0 : 1, irnf*sizeof(int));
                ifncntp = (int *)calloc(nfld, sizeof(int));
                mem_er((ifncntp == NULL) ? 0 : 1, nfld * sizeof(int));

                flev = fopen("level_ids.dat", "r");
                if(!flev){
                   printf("****ERROR****, unable to open level_ids.dat file\n\n");
                   exit(1);
                }

                for(i=0; i < irnf; i++) {
                    if(! fscanf(flev, "%d", ifncnt + i)){
                       printf("****ERROR****, inconsistent number of level ID's\n\n");
                       exit(1);
                    }

                    if(*(ifncnt + i) > nfld){
                       printf("****ERROR****, level Id. %d information does not exist for these tracks.\n\n", *(ifncnt + i));
                       exit(1);
                    }

                    if(*(ifncnt + i)){
                       if(*(nfwpos + *(ifncnt + i) - 1)){
                          *(ifncntp + i) = 0;
                          for(j=0; j < *(ifncnt + i) - 1; j++){
                              if(*(nfwpos + j)) *(ifncntp + i) += 3;
                              else *(ifncntp + i) += 1;
                          }
                        }
                        else{
                           printf("****WARNING****, selected field does not have positional information, setting to default.\n\n");
                           *(ifncnt + i) = 0;       
                        }
                    }
                }


                fclose(flev);

             }

/* Are fluxes required? */

             printf("Do you want fluxes with respect to another field, 'y' or 'n'?\n\n");
             scanf("\n");
             if(getchar() == 'y'){
                iflx = 1;
                printf("Use winds to compute KE or read from file, '0' for file and '1' for KE\n\n");
                scanf("%d", &iffl);
                if(iffl < 0 || iffl > 1){
                   printf("****ERROR****, incorrect option %d specified.\n\n", iffl);
                   exit(1);
                }
                if(!iffl){
                   printf("What is the filename for the additional field?\n\n");
                   scanf("%s", filrad);
                   fadd = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);

                   printf("Compute deviation from azimuthal average before computing flux, '0' for no or '1' for yes.\n\n");
                   scanf("%d", &idevazm);
                   if(idevazm){
                      alv = (double *)calloc(irnf, sizeof(double));
                      mem_er((alv == NULL) ? 0 : 1, irnf*sizeof(double));
                      nalv = (int *)calloc(irnf, sizeof(int));
                      mem_er((nalv == NULL) ? 0 : 1, irnf*sizeof(int));

                   }

                }
                else{
                   if(irel){
                      printf("Do you want to use the system relative winds to compute KE?, '0' for no, '1' for yes.\n\n");
                      scanf("%d", &irif);
                   }

                   printf("Do you want to combine the KE with another field? '0' for no and '1' for yes.\n\n");
                   scanf("%d", &ikecb);
                   if(ikecb < 0 || ikecb > 1){
                      printf("****ERROR****, incorrect option %d specified.\n\n", ikecb);
                      exit(1);
                   }
                   if(ikecb){
                      printf("Do you want the anomomaly over the region for this field, 'y' or 'n'\n\n");
                      scanf("\n");
                      if(getchar() == 'y') ike_cb_anom=1;
                      printf("What is the filename for the additional field?\n\n");
                      scanf("%s", filrad);
                      fadd = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);
                   }
           

                }

                if(!irif){
                   printf("Do you want the anomomaly over the region for the additional field, e.g. KE, 'y' or 'n'\n\n");
                   scanf("\n");
                   if(getchar() == 'y') ifanom=1;
                }

                printf("What is the radius to integrate over? Note, a radius larger than the maximum sampled radius \r\n"
                       "will be reset to the maximum sampled radius.                                                \n\n");
                scanf("%f", &nrad);
                asclr = FPI2 * (1.0 - cos(nrad * FP_PI)) * EARTH_RADIUS_M * EARTH_RADIUS_M;
                ascl = -asclr * RGAS / GGRAV;

                nrad = 90.0 - nrad;
                lscl = -EARTH_RADIUS_M * cos(nrad * FP_PI) * FPI2 / GGRAV;

                printf("ascl=%e, lscl=%e\n", ascl, lscl);
 
             }

             if(!iflx){
                printf("Do you want to compute advection with respect to another field , 'y' or 'n'.\r\n"
                       "Note this assumes the gradient components of the field are available.       \n\n");
                scanf("\n");
                if(getchar() == 'y'){
                   iadv = 1;
                   printf("What is the filename for the X(U)-component of the gradient of the additional field?\n\n");
                   scanf("%s", filrad);
                   fgrdx = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);
                   printf("What is the filename for the Y(V)-component of the gradient of the additional field?\n\n");
                   scanf("%s", filrad);
                   fgrdy = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);
                   printf("What is the radius to integrate over? Note, a radius larger than the maximum sampled radius \r\n"
                          "will be reset to the maximum sampled radius.                                                \n\n");
                   scanf("%f", &nrad);

/* scaling of advection term assumes we are working with gradients of geopotential in units of m^2 s^(-2) */

                   scladv = -FPI2 * (1.0 - cos(nrad * FP_PI)) * EARTH_RADIUS_M * EARTH_RADIUS_M / GGRAV;

                   nrad = 90.0 - nrad;
                }

             }

             printf("Specify level data for vertical integration, 'y' or 'n'.\n\n");
             scanf("\n");
             if(getchar() == 'y') {
               if(!lev_typ)levdat = leveldata(irnf);
               else lyes = 1;
             }

          }

       }

       else {
           printf("Do you want to combine with another field, 'y' or 'n'\n\n");
           scanf("\n");
           if(getchar() == 'y'){
              iaddf = 1;
              printf("****INFORMATION****, currently multiplication only.\n\n");
              printf("What is the filename for the additional field?\n\n");
              scanf("%s", filrad);
              fadd = open_radial_file(filrad, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);

              printf("Do you want the anomomaly over the region for the additional field, 'y' or 'n'\n\n");
              scanf("\n");
              if(getchar() == 'y') ifanom=1;

              printf("What is the radius to integrate over? Note, a radius larger than the maximum sampled radius \r\n"
                     "will be reset to the maximum sampled radius.                                                \n\n");
              scanf("%f", &nrad);
              asclr = FPI2 * (1.0 - cos(nrad * FP_PI)) * EARTH_RADIUS_M * EARTH_RADIUS_M;
              ascl = -asclr * RGAS / GGRAV;

              nrad = 90.0 - nrad;
              lscl = -EARTH_RADIUS_M * cos(nrad * FP_PI) * FPI2 / GGRAV;

              printf("Line scaling = %e area scaling = %e\n", lscl, ascl);

              printf("Do you want to normalize with pressure, for KE generation when fields combined are TEMP and OMEGA?\r\n"
                     "Input 'n' for no and 'y' for yes.                                                                 \n\n");
              scanf("\n");
              if(getchar() == 'y'){
                 ipnorm = 1;
                 if(!lev_typ)levdat = leveldata(irnf);
                 else lyes = 1;
              }
              else {
                 printf("Specify level data for vertical integration, 'y' or 'n'.\n\n");
                 scanf("\n");
                 if(getchar() == 'y') {
                    if(!lev_typ)levdat = leveldata(irnf);
                    else lyes = 1;
                 }
              }

           }
           else{
              printf("Perform integrations anyway, 'y' or 'n'\n\n");
              scanf("\n");
              if(getchar() == 'y'){
                 yinteg = 1;
                 printf("Do you want the anomomaly over the region for the additional field, 'y' or 'n'\n\n");
                 scanf("\n");
                 if(getchar() == 'y') ifanom=1;

                 printf("What is the radius to integrate over? Note, a radius larger than the maximum sampled radius \r\n"
                        "will be reset to the maximum sampled radius.                                                \n\n");
                 scanf("%f", &nrad);
                 asclr = FPI2 * (1.0 - cos(nrad * FP_PI)) * EARTH_RADIUS_M * EARTH_RADIUS_M;
                 ascl = -asclr * RGAS / GGRAV;

                 nrad = 90.0 - nrad;
                 lscl = -EARTH_RADIUS_M * cos(nrad * FP_PI) * FPI2 / GGRAV;

                 printf("Line scaling = %e area scaling = %e\n", lscl, ascl);

                 if(!lev_typ)levdat = leveldata(irnf);
                 else lyes = 1;
              } 
           }

       }

/* read in A and B coefficients for converting model levels to pressures */

       if(lyes){

          printf("What is the file that contains the A's and B's, this should have one more value than number of levels, \r\n"
                 "assuming levels are half levels and be those associated with the current model levels.                 \n\n");
          scanf("%s", flevnam);

          read_model_levs(flevnam, &aa, &bb, &nlab);

          if(nlab != irnf + 1){
             printf("****ERROR****, incorrect number of half levels.\n\n");
             exit(1);
          }

          printf("What is the file with sampled surface pressure?\n\n");
          scanf("%s", filsfp);

          irnft = irnf;
          irnf = 1;
          fsfp = open_radial_file(filsfp, &irtrn, &iptnum, &irnth, &irnr, &irnf, &idir, &ndsmth, &ifcnt, &itpadd, &iprojd);
          irnf = irnft;

       }

       if(levdat){
          vint = (double *)calloc(irnf-1, sizeof(double));
          mem_er((vint == NULL) ? 0 : 1, (irnf-1)*sizeof(double));
          for(i=0; i < irnf-1; i++) *(vint + i) = *(levdat + i) - *(levdat + i + 1);
       }


/* only do STD if not using variable pressure levels, i.e. model levels */
       if(!lyes){
          printf("Do you want STD as well as means for significance testing? 'y' or 'n'\n\n");
          scanf("\n");
          if(getchar() == 'y') istd = 1;
       }

       if(!irel){
          printf("Do you want anomalies with respect to the sampled region, 'y' or 'n'\n\n");
          scanf("\n");
          if(getchar() == 'y') ianom = 1;
       }

       if(ianom || ifanom || ike_cb_anom){
          printf("What radius is required for calculating the anomaly? Note, a radius larger than the maximum  \r\n"
                 "sampled radius will be reset to the maximum sampled radius.                                  \n\n");
          scanf("%f", &anom_rad);
          anom_rad = 90.0 - anom_rad;
       }

/*       if(iwind && ianom){
          printf("****WARNING****, no anomalies for winds, resetting option to no anomalies.\n\n");
          ianom = 0;
       }
*/

       irdim = irnr * irnth;

       sradf1 = (float *)calloc(irdim, sizeof(float));
       mem_er((sradf1 == NULL) ? 0 : 1, irdim * sizeof(float));
       slng = (float *)calloc(irnth, sizeof(float));
       mem_er((slng == NULL) ? 0 : 1, irnth * sizeof(float));
       slat = (float *)calloc(irnr, sizeof(float));
       mem_er((slat == NULL) ? 0 : 1, irnr * sizeof(float));
       coslat = (double *)calloc(irnr, sizeof(double));
       mem_er((coslat == NULL) ? 0 : 1, irnr * sizeof(double));
       fread(slng, irnth*sizeof(float), 1, frad1);
       fscanf(frad1, "%*c");
       fread(slat, irnr*sizeof(float), 1, frad1);
       fscanf(frad1, "%*c");

       arcll = (*(slng + 1) - *slng) * FP_PI;
       arclt = (*slat - *(slat + 1)) * FP_PI;

       if(iflx){
          tmp1 = (float *)calloc(irdim, sizeof(float));
          mem_er((tmp1 == NULL) ? 0 : 1, irdim * sizeof(float));
          tmp2 = (float *)calloc(irdim, sizeof(float));
          mem_er((tmp2 == NULL) ? 0 : 1, irdim * sizeof(float));
       }

       if(iadv){
          sgrdf1 = (float *)calloc(irdim, sizeof(float));
          mem_er((sgrdf1 == NULL) ? 0 : 1, irdim * sizeof(float));
          sgrdf2 = (float *)calloc(irdim, sizeof(float));
          mem_er((sgrdf2 == NULL) ? 0 : 1, irdim * sizeof(float));
       }

       if(lyes){
          sspres = (float *)calloc(irdim, sizeof(float));
          mem_er((sspres == NULL) ? 0 : 1, irdim * sizeof(float));
       }

       if(levdat || lyes){
          vsavg1 = (double *)calloc(irdim, sizeof(double));
          mem_er((vsavg1 == NULL) ? 0 : 1, irdim*sizeof(double));
          vsavg2 = (double *)calloc(irdim, sizeof(double));
          mem_er((vsavg2 == NULL) ? 0 : 1, irdim*sizeof(double));
       }

       if(ianom || ifanom){
          iranom = irnr;
          if(anom_rad <= *(slat + irnr - 1)) {anom_rad = *(slat + irnr - 1); iranom = irnr;}
          else{
             for(i=0; i < irnr; i++){
                 if(*(slat + i) < anom_rad) {anom_rad = *(slat + i - 1); iranom = i; break;}
             }
          }
       }

       ncount = (float **)calloc(irnf, sizeof(float *));
       mem_er((ncount == NULL) ? 0 : 1, irnf * sizeof(float));
       for(i=0; i < irnf; i++){
           *(ncount + i) = (float *)calloc(irdim, sizeof(float));
           mem_er((*(ncount + i) == NULL) ? 0 : 1, irdim * sizeof(float));
           for(j=0; j < irdim; j++) *(*(ncount + i) + j) = 0.0;
       }

       aavg = (double **)calloc(irnf, sizeof(double *));
       mem_er((aavg == NULL) ? 0 : 1, irnf*sizeof(double *));

       for(i=0; i < irnf; i++){
          *(aavg + i) = (double *)calloc(irdim, sizeof(double));
          mem_er((*(aavg + i) == NULL) ? 0 : 1, irdim*sizeof(double));

       }

       if(istd){
          astd = (double **)calloc(irnf, sizeof(double *));
          mem_er((astd == NULL) ? 0 : 1, irnf*sizeof(double *));

          for(i=0; i < irnf; i++){
             *(astd + i) = (double *)calloc(irdim, sizeof(double));
             mem_er((*(astd + i) == NULL) ? 0 : 1, irdim*sizeof(double));

          }

       }

       if((icomp_order && iflx) || iaddf || iadv || yinteg){

          if(!iadv){
             sadd = (float *)calloc(irdim, sizeof(float));
             mem_er((sadd == NULL) ? 0 : 1, irdim * sizeof(float));

             if(ikecb){
               sadd2 = (float *)calloc(irdim, sizeof(float));
               mem_er((sadd2 == NULL) ? 0 : 1, irdim * sizeof(float));
             }
          }

          linteg1 = (double **)calloc(irnf, sizeof(double *));
          mem_er((linteg1 == NULL) ? 0 : 1, irnf * sizeof(double *));

          rintg1 = (double *)calloc(irnf, sizeof(double));
          mem_er((rintg1 == NULL) ? 0 : 1, irnf * sizeof(double));

          lvtmp = (double *)calloc(irnf, sizeof(double));
          mem_er((lvtmp == NULL) ? 0 : 1, irnf * sizeof(double));

          for(i=0; i < irnf; i++){
              linteg1[i] = (double *)calloc(irnr, sizeof(double));
              mem_er((linteg1[i] == NULL) ? 0 : 1, irnr * sizeof(double));

          }

          if(iflx){

             linteg2 = (double **)calloc(irnf, sizeof(double *));
             mem_er((linteg2 == NULL) ? 0 : 1, irnf * sizeof(double *));

             linteg3 = (double **)calloc(irnf, sizeof(double *));
             mem_er((linteg3 == NULL) ? 0 : 1, irnf * sizeof(double *));

             for(i=0; i < irnf; i++){
                 linteg2[i] = (double *)calloc(irnr, sizeof(double));
                 mem_er((linteg2[i] == NULL) ? 0 : 1, irnr * sizeof(double));
                 linteg3[i] = (double *)calloc(irnr, sizeof(double));
                 mem_er((linteg3[i] == NULL) ? 0 : 1, irnr * sizeof(double));
             }

             rintg2 = (double *)calloc(irnf, sizeof(double));
             mem_er((rintg2 == NULL) ? 0 : 1, irnf * sizeof(double));

             rintg3 = (double *)calloc(irnf, sizeof(double));
             mem_er((rintg3 == NULL) ? 0 : 1, irnf * sizeof(double));

          }

          irr = irnr;
          if(nrad <= *(slat + irnr - 1)) {nrad = *(slat + irnr - 1); irr = irnr;}
          else{
             for(i=0; i < irnr; i++){
                 if(*(slat + i) < nrad) {nrad = *(slat + i - 1); irr = i; break;}
             }
          }

          printf("****INFORMATION****, current maximum radius is %f at grid point %d\n\n", nrad, irr);

       }

       if(iwind){

          slng2 = (float *)calloc(irnth, sizeof(float));
          mem_er((slng2 == NULL) ? 0 : 1, irnth * sizeof(float));
          slat2 = (float *)calloc(irnr, sizeof(float));
          mem_er((slat2 == NULL) ? 0 : 1, irnr * sizeof(float));

          fread(slng2, irnth*sizeof(float), 1, frad2);
          fscanf(frad2, "%*c");
          fread(slat2, irnr*sizeof(float), 1, frad2);
          fscanf(frad2, "%*c");

          for(i=0; i < irnr; i++){
              if(fabs(*(slat + i) - *(slat2 + i)) > GTOL){
                 printf("****ERROR****, grids do not match in radial direction\n\n");
                 exit(1);
              }
          }

          for(i=0; i < irnth; i++){
              if(fabs(*(slng + i) - *(slng2 + i)) > GTOL){
                 printf("****ERROR****, grids do not match in tangential direction\n\n");
                 exit(1);
              }
          }

          aavg1 = (double **)calloc(irnf, sizeof(double *));
          mem_er((aavg1 == NULL) ? 0 : 1, irnf*sizeof(double *));

          aavg2 = (double **)calloc(irnf, sizeof(double *));
          mem_er((aavg2 == NULL) ? 0 : 1, irnf*sizeof(double *));

          for(i=0; i < irnf; i++){
             *(aavg1 + i) = (double *)calloc(irdim, sizeof(double));
             mem_er((*(aavg1 + i) == NULL) ? 0 : 1, irdim*sizeof(double));
             *(aavg2 + i) = (double *)calloc(irdim, sizeof(double));
             mem_er((*(aavg2 + i) == NULL) ? 0 : 1, irdim*sizeof(double));
          }

          if(istd){
             astd1 = (double **)calloc(irnf, sizeof(double *));
             mem_er((astd1 == NULL) ? 0 : 1, irnf*sizeof(double *));

             astd2 = (double **)calloc(irnf, sizeof(double *));
             mem_er((astd2 == NULL) ? 0 : 1, irnf*sizeof(double *));

             for(i=0; i < irnf; i++){
                *(astd1 + i) = (double *)calloc(irdim, sizeof(double));
                mem_er((*(astd1 + i) == NULL) ? 0 : 1, irdim*sizeof(double));
                *(astd2 + i) = (double *)calloc(irdim, sizeof(double));
                mem_er((*(astd2 + i) == NULL) ? 0 : 1, irdim*sizeof(double));
             }
          }

          if(icomp_order){
             sradf2 = (float *)calloc(irdim, sizeof(float));
             mem_er((sradf2 == NULL) ? 0 : 1, irdim * sizeof(float));

             if(iflx){
                savg1 = (double **)calloc(irnf, sizeof(double *));
                mem_er((savg1 == NULL) ? 0 : 1, irnf*sizeof(double *));

                savg2 = (double **)calloc(irnf, sizeof(double *));
                mem_er((savg2 == NULL) ? 0 : 1, irnf*sizeof(double *));

                savg3 = (double **)calloc(irnf, sizeof(double *));
                mem_er((savg3 == NULL) ? 0 : 1, irnf*sizeof(double *));

                for(i=0; i < irnf; i++){
                   *(savg1 + i) = (double *)calloc(irdim, sizeof(double));
                   mem_er((*(savg1 + i) == NULL) ? 0 : 1, irdim*sizeof(double));
                   *(savg2 + i) = (double *)calloc(irdim, sizeof(double));
                   mem_er((*(savg2 + i) == NULL) ? 0 : 1, irdim*sizeof(double));
                   *(savg3 + i) = (double *)calloc(irdim, sizeof(double));
                   mem_er((*(savg3 + i) == NULL) ? 0 : 1, irdim*sizeof(double));
                }

                vke = (double *)calloc(irnf, sizeof(double));
                mem_er((vke == NULL) ? 0 : 1, irnf*sizeof(double));

             }

             else if (iadv){
                savg1 = (double **)calloc(irnf, sizeof(double *));
                mem_er((savg1 == NULL) ? 0 : 1, irnf*sizeof(double *));
                for(i=0; i < irnf; i++){
                   *(savg1 + i) = (double *)calloc(irdim, sizeof(double));
                   mem_er((*(savg1 + i) == NULL) ? 0 : 1, irdim*sizeof(double));
                }
             }

             vecg = (VEC *)calloc(irdim, sizeof(VEC));
             mem_er((vecg == NULL) ? 0 : 1, irdim*sizeof(VEC));

             newv = (VEC *)calloc(irdim, sizeof(VEC));
             mem_er((newv == NULL) ? 0 : 1, irdim*sizeof(VEC));

             etht = (VEC *)calloc(irdim, sizeof(VEC));
             mem_er((etht == NULL) ? 0 : 1, irdim*sizeof(VEC));

             ephi = (VEC *)calloc(irdim, sizeof(VEC));
             mem_er((ephi == NULL) ? 0 : 1, irdim*sizeof(VEC));

             erad = (VEC *)calloc(irdim, sizeof(VEC));
             mem_er((erad == NULL) ? 0 : 1, irdim*sizeof(VEC));

             etan = (VEC *)calloc(irdim, sizeof(VEC));
             mem_er((etan == NULL) ? 0 : 1, irdim*sizeof(VEC));

             for(i=0; i < irnr; i++){
                 rr = FP_PI2 - *(slat + i) * FP_PI;
                 if(rr < 0.0) rr = 0.0;
                 sincos(rr, &s1, &c1);
                 for(j=0; j < irnth; j++){
                     sincos(*(slng + j) * FP_PI, &s2, &c2);
                     vt = vecg + i * irnth + j;
                     vt->x = s1 * c2;
                     vt->y = s1 * s2;
                     vt->z = c1;
                 }
             }

             if(irel){
                erdd = (VEC *)calloc(irdim, sizeof(VEC));
                mem_er((erdd == NULL) ? 0 : 1, irdim*sizeof(VEC));

                etnn = (VEC *)calloc(irdim, sizeof(VEC));
                mem_er((etnn == NULL) ? 0 : 1, irdim*sizeof(VEC));

                for(i=0; i < irdim; i++){
                    etnn[i].x= -vecg[i].y;
                    etnn[i].y = vecg[i].x;
                    norm = sqrt(dotp(&etnn[i], &etnn[i]));
                    normv(&etnn[i], norm);

                    crosp(&etnn[i], &vecg[i], &erdd[i]);
                    norm = sqrt(dotp(&erdd[i], &erdd[i]));
                    normv(&erdd[i], norm);   
                }

             }

          }

       }

       for(i=0; i < irnr; i++) *(coslat + i) = cos(*(slat + i) * FP_PI);



/* calculate data block size */

       place1 = ftello(frad1);
       fread(sradf1, irdim*sizeof(float), 1, frad1);
       fscanf(frad1, "%*c");
       place2 = ftello(frad1);
       blklen = place2 - place1;
       fseeko(frad1, place1, SEEK_SET);

       printf("Do you want field extraction and averaging at:-        \r\n"
              "the maximum or minimum intensity, '0'                  \r\n"
              "the maximum or minimum tendency,  '1'                  \r\n"
              "the maximum vertical gradiant,    '2'                  \r\n"
              "or a particular time step,        '3'                  \n\n");

       scanf("%d", &iext);
       if(iext < 0 || iext > 3){
          printf("****ERROR****, incorect input %d\n", iext);
          exit(1);
       }
       if(iext == 0 || iext == 1){
          printf("Do you want minimum or maximum? Input '0' or '1'. \n\n");
          scanf("%d", &imxmn);
          if(imxmn < 0 || imxmn > 1){
             printf("****ERROR****, incorect input %d\n", imxmn);
             exit(1);
          }
       }
       else if(iext == 3){
          if(tracks->time){
             printf("What time do you want to composite at, YYYYMMDDHH?\n\n");
             scanf("%ld", &ctime);
          }
          else {
             printf("What time step do you want to composite at?\n\n");
             scanf("%ld", &ctime);
          }
       }

       printf("What offset time step is required relative to chosen time step?\n\n");
       scanf("%d", &toff);

       printf("Use all tracks or select tracks, '0' for all or '1' for select?\n\n");
       scanf("%d", &iselt);
       if(iselt < 0 || iselt > 1){
          printf("****ERROR****, wrong input for track selection, exiting.\n\n");
          exit(1);
       }
       if(iselt){
          printf("Number of tracks is %d which start and end track do you want?\n\n", trnum);
          scanf("%d %d", &itrst, &itren);
          --itrst;
          if(itrst < 0 || itren > trnum){
             printf("****ERROR****, chosen tracks are out of range.\n\n");
             exit(1);
          }
       }
       else{
         itrst = 0;
         itren = trnum;
       }

    }

    printf("Specify threshold for identification.\n\n");
    scanf("%f", &thresh);

    if(nff){
       printf("****INFORMATION****, there are additional fields  \r\n"
              "                     specify how to use them.     \n\n");
       printf("****INFORMATION****, at the moment only have the tests f > thresh and f0 - f1 > fth\n\n");
       printf("Which additional field is required for intensity thresholding, input '0' for default.\n\n");
       scanf("%d", &ift);
       if(ift < 0 || ift > nff){
          printf("****ERROR****, not a valid field Id.\n\n");
          exit(1);
       }
       ifft = 0;
       if(ift){
          for(i=0; i < ift; i++){
              if(*(nfwpos + i)) ifft += 3;
              else ifft += 1;
          }
          --ifft;
       }

       
       if(ift){
          if(*(nfwpos + ift - 1)){
             printf("Additional field has locational information, due you want to exclude missing values, \r\n"
                    "'0' for no and '1' for yes.                                                          \n\n");
             scanf("%d", &imiss);
          }
       }

       printf("Do you want to test vertical structure, 'y' or 'n'\n\n");
       scanf("\n");
       itv = getchar();

       if(ilm){
          printf("Do you want to check for land points, 'y' or 'n', if land then point is not used for TC identification.\n\n");
          scanf("\n");
          irog = getchar();
       }

       if(itv == 'y'){

          printf("What value do you want for 'fth'.\n\n");
          scanf("%f", &fth);
          printf("What value do you want for the number of points that satisfy fth.\n\n");
          scanf("%d", &nth);
          printf("Specify which pair of additonal fields to use, '0' implies use of original\r\n"
                 "vortcity center.                                                          \n\n");
          scanf("%d %d", &if1, &if2);
          if((if1 < 0 || if1 > nff) || (if2 < 0 || if2 > nff)){
             printf("****ERROR****, invalid field Id. chosen\n\n");
             exit(1);
          }

          iff1 = 0;
          if(if1){
             if(! *(nfwpos + if1 - 1)){
                printf("****ERROR****, not possible to check full vertical structure.\n\n");
                nvert = 1;
             }
             for(i=0; i < if1; i++){
                if(*(nfwpos + i)) {iff1 += 3; ++nl1;}
                else iff1 += 1;
             }
             --iff1;
          }

          iff2 = 0;
          if(if2){
             for(i=0; i < if2; i++){
                if(*(nfwpos + i)) {iff2 += 3; ++nl2;}
                else iff2 += 1;
                if(i+1 >= if1 && ! *(nfwpos + i)){
                   printf("****ERROR****, not possible to check full vertical structure.\n\n");
                   nvert = 1;
                }
             }
             --iff2;
          }

          nlp = nl2 - nl1 + 1;

          printf("Number of levels used is %d\n\n", nlp);

          printf("Do you want to test full vertical structure is present upto chosen upper level, 'y' or 'n' \r\n"
                 "This will use all additional fields assuming they are in vertical order, and that postional\r\n"
                 "information is present.                                                                    \n\n");
          scanf("\n");
          if(getchar() == 'y') {
            ifs = 1;
            if(nvert) exit(1);
          }

       }

       printf("If speed data is present which field is it, input '0' for none?\n\n");
       scanf("%d", &ispd);
       if(ispd){
          ifsp = 0;
          for(i=0; i < ispd; i++){
             if(*(nfwpos + i)) ifsp += 3;
             else ifsp += 1;
          }
          --ifsp;

       }

       printf("If precipitation data is available which field is it, input '0' for none?\n\n");
       scanf("%d", &isprec);
       if(isprec){
          ifsprec = 0;
          for(i=0; i < isprec; i++){
             if(*(nfwpos + i)) ifsprec += 3;
             else ifsprec += 1;
          }
          --ifsprec;

          fprc = fopen("precip.dat", "w");
          if(!fprc){
             printf("****ERROR****, can't open file precip.dat.\n\n");
             exit(1);
          }

       }

    }

    else{
       printf("****WARNING****, no additional fields are present for looking at vertical structure.\n\n");
       printf("                 only intensity thresholding available.                             \n\n");
    }

    ftim = fopen("timout.dat", "w");
    if(!ftim){
       printf("****WARNING****, can't open file timout.dat for write.\n\n");
    }

    fdout = fopen("fdif.dat", "w");
    if(!fdout){
       printf("****WARNING****, can't open file fdif.dat for write.\n\n");
    }

/* setup additional sampling of radial max, min and range */

    if(irad){

       tmpdat = (float *)calloc(irdim, sizeof(float));
       mem_er((tmpdat == NULL) ? 0 : 1, irdim * sizeof(float));

       rmposx = (off_t **)calloc(irnf, sizeof(off_t *));
       mem_er((rmposx == NULL) ? 0 : 1, irnf * sizeof(off_t *));

       rmposn = (off_t **)calloc(irnf, sizeof(off_t *));
       mem_er((rmposn == NULL) ? 0 : 1, irnf * sizeof(off_t *));

       rrpos = (off_t **)calloc(irnf, sizeof(off_t *));
       mem_er((rrpos == NULL) ? 0 : 1, irnf * sizeof(off_t *));

       frmax = fopen("frmax.dat", "w");
       if(!frmax){
          printf("****ERROR****, can't open file frmax.dat for write.\n\n");
          exit(1);
       }

       frmin = fopen("frmin.dat", "w");
       if(!frmax){
          printf("****ERROR****, can't open file frmin.dat for write.\n\n");
          exit(1);
       }

       frange = fopen("frange.dat", "w");
       if(!frange){
          printf("****ERROR****, can't open file frange.dat for write.\n\n");
          exit(1);
       }

       for(i=0; i < irnf; i++){
          *(rmposx + i) = (off_t *)calloc(irnth, sizeof(off_t));
          mem_er(((rmposx + i) == NULL) ? 0 : 1, irnth * sizeof(off_t));
          *(rmposn + i) = (off_t *)calloc(irnth, sizeof(off_t));
          mem_er(((rmposn + i) == NULL) ? 0 : 1, irnth * sizeof(off_t));
          *(rrpos + i) = (off_t *)calloc(irnth, sizeof(off_t));
          mem_er(((rrpos + i) == NULL) ? 0 : 1, irnth * sizeof(off_t));

          fprintf(frmax, "Level = %d\n", i+1);
          fprintf(frmin, "Level = %d\n", i+1);
          fprintf(frange, "Level = %d\n", i+1);
          for(j=0; j < irnth; j++){
              fprintf(frmax, "%8.4f & ", *(slng + j));
              fprintf(frmin, "%8.4f & ", *(slng + j));
              fprintf(frange, "%8.4f & ", *(slng + j));
              *(*(rmposx + i) + j) = ftello(frmax);
              *(*(rmposn + i) + j) = ftello(frmin);
              *(*(rrpos + i) + j) = ftello(frange);
              for(k=0; k < trnum; k++) {
                  fprintf(frmax, "%8.4f %e & ", 0.0, 0.0);
                  fprintf(frmin, "%8.4f %e & ", 0.0, 0.0);
                  fprintf(frange, "%8.4f %8.4f & ", 0.0, 0.0);
              }
              fprintf(frmax, "\n");
              fprintf(frmin, "\n");
              fprintf(frange, "\n");
          }
       }

/* read threshold value for radial range */
       fthresh = fopen("thresh.dat", "r");
       if(!fthresh){
          printf("****WARNING****, can't open file thresh.dat for read.\n\n");
       }
       else {
          fscanf(fthresh, "%f", &rthresh);
          fclose(fthresh);
       }

    }

/* open file for writing rotation angle and speeds */

    fdir = fopen("dir_speed.dat", "w");
    if(!fdir){
       printf("****WARNING****, can't open file dir_speed.dat for write.\n\n");
       exit(1);
    }

    for(i=0; i < itrst; i++) nwtr += (tracks + i)->num;

    printf("Start %d\n", nwtr);

    for(i=itrst; i < itren; i++){
       atr = tracks + i;
       atr2 = tracks2 + i;
       fp = atr->trpt;

       fdifmax = 0.0;

       nctr = nwtr;
       ifld = -1;

       fmxmn = (imxmn) ? SMALL : LARGE;

       frst = last = -1;

       if(ilm){
          fmin = 100000.0;
          for(j=0; j < lglng; j++){
              fff = fabs(*(glng + j) - fp->xf);
              if(fff < fmin) {fmin = fff; ilng = j;}
          }

          fmin = 100000.0;
          for(j=0; j < lglat; j++){
              fff = fabs(*(glat + j) - fp->yf);
              if(fff < fmin) {fmin = fff; ilat = j;}
          }

          igen = 1;
          if(*(ilms + ilat * lglng + ilng)){
            if(fp->yf < -20.0 || fp->yf > 20.0) igen = 0;
          }
          else{
            if(fp->yf < -30.0 || fp->yf > 30.0) igen = 0;
          }

          if(!igen){
             for(j=0; j < atr->num; j++){
                free((atr->trpt + j)->add_fld);
                free((atr2->trpt + j)->add_fld);
             }
             free(atr->trpt);
             free(atr2->trpt);
             atr->num = 0;
             atr2->num = 0;
             continue;
          }

       }
       else if(igg){
          if(fp->yf < glat1 || fp->yf > glat2) {
             for(j=0; j < atr->num; j++){
                 free((atr->trpt + j)->add_fld);
                 free((atr2->trpt + j)->add_fld);
             }
             free(atr->trpt);
             free(atr2->trpt);
             atr->num = 0;
             atr2->num = 0;
             continue;
          }
       }

       if(ireg){
          fmax = 0.0;
          inreg = 0;
          for(j=0; j < atr->num; j++){
             fp = atr->trpt + j;
             if(fp->zf > fmax){
                fmax = fp->zf;
                lngm = fp->xf;
                latm = fp->yf;
                if(atr->time) tim = fp->time;
                else tim = fp->fr_id;
             }

          }

          if(lng1 > lng2){
             if((lat2 - latm) * (latm - lat1) >= 0.0 && 
                ((lng2 - lngm) * lngm >= 0.0 || (360.0 - lngm) * (lngm - lng1) >= 0.0)) inreg = 1;
          }

          else if((lng2 - lngm) * (lngm - lng1) >= 0.0 && (lat2 - latm) * (latm - lat1) >= 0.0) inreg = 1;

          if(inreg && tim >= tim1 && tim <= tim2) inreg = 1;
          else inreg = 0; 

          if(!inreg){
             for(j=0; j < atr->num; j++){
                 free((atr->trpt + j)->add_fld);
                 free((atr2->trpt + j)->add_fld);
             }
             free(atr->trpt);
             free(atr2->trpt);
             atr->num = 0;
             atr2->num = 0;
             continue;
          }

       }

       istc = 0;
       nvr = 0;

       for(j=0; j < atr->num; j++){

           fp = atr->trpt + j;

           if(irog == 'y'){
              if(orog_test(glng, glat, lglng, lglat, ilms, fp->xf, fp->yf)) {istc = 0; continue;}
           }

           if(!ift) str = fp->zf;
           else {
              str = fp->add_fld[ifft];
              if(str > ADD_CHECK) continue;
              if(imiss && fp->add_fld[ifft - 2] > ADD_CHECK) continue;
           }

           if(ifs){
              if(!if1){
                 iv = 1;
                 nn = 0;
                 for(k=0; k < if2; k++){
                     if(fp->add_fld[nn] < ADD_CHECK) ++iv; 
                     nn += 3;
                 }
              }
              else {
                 iv = 0;
                 nn = iff1 - 2;
                 for(k=if1-1; k < if2; k++){
                     if(fp->add_fld[nn] < ADD_CHECK) ++iv; 
                     nn += 3;
                 }
              }
           }
           else iv = nlp;

           if(itv == 'y') {
              if(!if1)
                 fdif = fp->zf - fp->add_fld[iff2];
              else
                 fdif = fp->add_fld[iff1] - fp->add_fld[iff2];
              if(fdif > fdifmax) fdifmax = fdif;
           }
           else fdif = LARGE;

           if(str >= thresh && fdif >= fth && iv == nlp){

              ++istc;
              if(frst < 0) frst = j;
              last = j;
              if(istc > nvr) nvr = istc;
              if(irad){
                 if(!iext){
                    if(imxmn){
                       if(str > fmxmn) {fmxmn = str; ifld = j;}
                    }
                    else {
                       if(str < fmxmn) {fmxmn = str; ifld = j;}
                    }
                 }
                 else if(iext == 1){
                    if(j < atr->num - 1){

                       if(!ift) str2 = (fp + 1)->zf;
                       else {
                          str2 = (fp + 1)->add_fld[ifft];
                          if(str2 > ADD_CHECK) continue;
                          if(imiss && (fp + 1)->add_fld[ifft - 2] > ADD_CHECK) continue;
                       }

                       tend = str2 - str;
                       if(imxmn){
                          if(tend > fmxmn) {fmxmn = tend; ifld = j;}
                       }
                       else {
                          if(tend < fmxmn) {fmxmn = tend; ifld = j;}
                       } 
                          
                    }
                 }
                 else if(iext == 2){
                   if(fdif > fmxmn) {fmxmn = fdif; ifld = j;}
                 }
                 else if(iext == 3){
                   if(atr->time){
                     if(fp->time == ctime) ifld = j;
                   }
                   else{
                     if((long int)(fp->fr_id) == ctime) ifld = j;
                   }
                 }
              }
           }
           else istc = 0;
       }

/*      printf("time point = %d, %d\n", ifld, (atr->trpt + ifld)->fr_id); */

       nwtr += atr->num;

       if(nvr < nth || ifld < 0) {
          for(j=0; j < atr->num; j++){
              free((atr->trpt + j)->add_fld);
              free((atr2->trpt + j)->add_fld);
          }
          free(atr->trpt);
          free(atr2->trpt);
          atr->num = 0;
          atr2->num = 0;
       }
       else {

          ifld += toff;

/*          fprintf(ftim, "%5d\n", (atr->trpt + ifld)->fr_id); */
          fprintf(ftim, "%5d\n", ifld + 1);
          fprintf(fdout, "%f\n", fdifmax);

          if(irad && (ifld >= 0 && ifld < atr->num)){

             if(lyes){
                fseeko(fsfp, (nctr + ifld) * blklen + place1, SEEK_SET);
                fread(sspres, irdim*sizeof(float), 1, fsfp);

             }

             fseeko(frad1, (nctr + ifld) * irnf * blklen + place1, SEEK_SET);

             if(iaddf) fseeko(fadd, (nctr + ifld) * irnf * blklen + place1, SEEK_SET);

             if(!icomp_order){
                for(j=0; j < irnf; j++){

                   if(lyes) model_level_diff(&ad, &bd, aa, bb, j, irnf);

                   fread(sradf1, irdim*sizeof(float), 1, frad1);
                   fscanf(frad1, "%*c");
                   if(ianom) anomaly(sradf1, coslat, irnr, irnth, iranom);

                   if(iaddf){
                      fread(sadd, irdim*sizeof(float), 1, fadd);
                      fscanf(fadd, "%*c");
                      if(ifanom) anomaly(sadd, coslat, irnr, irnth, iranom);

/* combine primary and secondary fields, could add more options here */

                      for(k=0; k < irdim; k++) {
                         if(*(sradf1 + k) > ADD_CHECK || *(sadd + k) > ADD_CHECK) continue;
                         *(sradf1 + k) = *(sradf1 + k) * *(sadd + k);
                         if(ipnorm) {
                            if(levdat) *(sradf1 + k) /= *(levdat + j);
                            else  *(sradf1 + k) /= *(aa + j) + *(bb + j) * *(sspres + k);
                         }
                      }
                   }

/* write radial distributions */

                   radial_extreme(frmax, frmin, frange, sradf1, slat, rthresh, irnth, irnr, *(rmposx + j), *(rmposn + j), *(rrpos + j));

                   for(k=0; k < irdim; k++) {
                       if(*(sradf1 + k) > ADD_CHECK) continue;
                       if(lyes) *(*(aavg + j) + k) += *(sradf1 + k) * (ad + bd * *(sspres + k));
                       else *(*(aavg + j) + k) += *(sradf1 + k);
                       if(istd) *(*(astd + j) + k) += *(sradf1 + k) * *(sradf1 + k);
                       *(*(ncount + j) + k) += 1.0;
                   }

                }

             }

             if(iwind){
                fseeko(frad2, (nctr + ifld) * irnf * blklen + place1, SEEK_SET);

                if(iflx && (!iffl || ikecb)) fseeko(fadd, (nctr + ifld) * irnf * blklen + place1, SEEK_SET);

                if(iadv){
                   fseeko(fgrdx, (nctr + ifld) * irnf * blklen + place1, SEEK_SET);
                   fseeko(fgrdy, (nctr + ifld) * irnf * blklen + place1, SEEK_SET);
                }

                if(!icomp_order){
                   for(j=0; j < irnf; j++){
                      if(lyes) model_level_diff(&ad, &bd, aa, bb, j, irnf);
                      fread(sradf1, irdim*sizeof(float), 1, frad2);
                      fscanf(frad2, "%*c");

                      if(ianom) anomaly(sradf1, coslat, irnr, irnth, iranom);

                      for(k=0; k < irdim; k++) {
                         if(*(sradf1 + k) > ADD_CHECK) continue;
                         if(lyes) *(*(aavg2 + j) + k) += *(sradf1 + k) * (ad + bd * *(sspres + k));
                         else *(*(aavg2 + j) + k) += *(sradf1 + k);
                         if(istd) *(*(astd2 + j) + k) += *(sradf1 + k) * *(sradf1 + k);
                      }

                   }

                }

             }

             if(iwind && icomp_order){

                fpc = atr->trpt + ifld;
                xx = fpc->xf * FP_PI;
                yy = FP_PI2 - fpc->yf * FP_PI;
                if(yy < 0.) yy = 0.0;
                sincos(xx, &s1, &c1);
                sincos(yy, &s2, &c2);

                pt1.x = fpc->pp[0];
                pt1.y = fpc->pp[1];
                pt1.z = fpc->pp[2];

                setup_rotmatrix(s1, c1, s2, c2, rot);

/* has data has been sampled on a grid that has been rotated to direction of storm */

                if(idir){

                   vsdir.x = -c1 * c2;
                   vsdir.y = -s1 * c2;
                   vsdir.z = s2;

/* compute system velocity for relative velocities */

                   if(irel) pdirl = prop_speed(atr, ifld) * EARTH_RADIUS_M / tstp;

/* calculate directional vector */
                   pt1.x = fpc->pp[0];
                   pt1.y = fpc->pp[1];
                   pt1.z = fpc->pp[2];

                   arot =   dirangle(atr, &cc1, &ss1, ifld, ndsmth);

                   fprintf(fdir, "%e %e\n", arot/FP_PI, pdirl);

                   norm = dotp(&vsdir, &pt1);
                   crosp(&vsdir, &pt1, &vt1);
                   mulv(&vt1, &vt1, ss1);
                   mulv(&pt1, &vt2, ((1.0 - cc1) * norm));
                   mulv(&vsdir, &vt3, cc1);
                   addv(&vt1, &vt2, &vn);
                   addv(&vn, &vt3, &tvec);

                }
               
                if(!ifcnt)
                   setup_radvec(rot, vecg, &pt1, newv, etht, ephi, erad, etan, ss1, cc1, irdim, idir);

                for(j=0; j < irnf; j++){

                   if(lyes) model_level_diff(&ad, &bd, aa, bb, j, irnf);

                   if(ifcnt){
                      find_last(fpc, &fpt, ifncnt, ifncntp, j, itpadd);

                      xx = (fpt.x).xy * FP_PI;
                      yy = FP_PI2 - (fpt.y).xy * FP_PI;

/* determine sampling point at this level */

                      if(iprojd){

                         sincos(xx, &s1, &c1);
                         sincos(yy, &s2, &c2);

                         ptt.x = c1 * s2;
                         ptt.y = s1 * s2;
                         ptt.z = c2;

                         if(idir) fvec(&pt1, &tvec, &ptt, &vn);
                         else {
                           if(ifld < atr->num - 1){
                              pt1.x = fpc->pp[0];
                              pt1.y = fpc->pp[1];
                              pt1.z = fpc->pp[2];
                              pt2.x = (fpc+1)->pp[0];
                              pt2.y = (fpc+1)->pp[1];
                              pt2.z = (fpc+1)->pp[2];
                           }
                           else {
                              pt1.x = (fpc-1)->pp[0];
                              pt1.y = (fpc-1)->pp[1];
                              pt1.z = (fpc-1)->pp[2];
                              pt2.x = fpc->pp[0];
                              pt2.y = fpc->pp[1];
                              pt2.z = fpc->pp[2];
                           }

                           fvec(&pt1, &pt2, &ptt, &vn);
                         }

                         yy = acos(vn.z);
                         xx = atan2(vn.y, vn.x);
                         if(xx < 0.0) xx = FPI2 + xx;

                      }

                      sincos(xx, &s1, &c1);
                      sincos(yy, &s2, &c2);

                      setup_rotmatrix(s1, c1, s2, c2, rot);

                      pt1.x = c1 * s2;
                      pt1.y = s1 * s2;
                      pt1.z = c2;

                      setup_radvec(rot, vecg, &pt1, newv, etht, ephi, erad, etan, ss1, cc1, irdim, idir);

                   }


                   fread(sradf1, irdim*sizeof(float), 1, frad1);
                   fscanf(frad1, "%*c");
                   fread(sradf2, irdim*sizeof(float), 1, frad2);
                   fscanf(frad2, "%*c");

                   if(iflx){
                      if(!iffl){
                         fread(sadd, irdim*sizeof(float), 1, fadd);
                         fscanf(fadd, "%*c");
                         if(ifanom) anomaly(sadd, coslat, irnr, irnth, iranom);
                      }
                      else if(!irif){
                         memcpy(tmp1, sradf1, irdim * sizeof(float));
                         memcpy(tmp2, sradf2, irdim * sizeof(float));
                         if(ifanom) {
                            anomaly(tmp1, coslat, irnr, irnth, iranom);
                            anomaly(tmp2, coslat, irnr, irnth, iranom);
                         }
                         for(k=0; k < irdim; k++) {
                            uu = *(tmp1 + k);
                            vv = *(tmp2 + k);
                            if(uu < ADD_CHECK && vv < ADD_CHECK)
                               if(lyes) *(sadd + k) = 0.5 * (uu * uu + vv * vv) * (ad + bd * *(sspres + k));
                               else *(sadd + k) = 0.5 * (uu * uu + vv * vv);
                            else
                               *(sadd + k) = ADD_UNDEF;

                         }

                      }
                      
                   }

                   else if(iadv){
                      fread(sgrdf1, irdim*sizeof(float), 1, fgrdx);
                      fscanf(fgrdx, "%*c");
                      fread(sgrdf2, irdim*sizeof(float), 1, fgrdy);
                      fscanf(fgrdy, "%*c");
                   }

                   if(ikecb){
                      fread(sadd2, irdim*sizeof(float), 1, fadd);
                      fscanf(fadd, "%*c");
                      if(ike_cb_anom) anomaly(sadd2, coslat, irnr, irnth, iranom);
                   }

                   if(ianom){
                      anomaly(sradf1, coslat, irnr, irnth, iranom);
                      anomaly(sradf2, coslat, irnr, irnth, iranom);
                   }

                   for(k=0; k < irdim; k++){
                       uu = *(sradf1 + k);
                       vv = *(sradf2 + k);

                       if(uu > ADD_CHECK || vv > ADD_CHECK) continue;
                       *(*(ncount + j) + k) += 1.0;

                       mulv(&ephi[k], &vx, uu);
                       mulv(&etht[k], &vy, (-vv));
                       addv(&vx, &vy, &vel);

                       uu = dotp(&etan[k], &vel);
                       vv = dotp(&erad[k], &vel);

                       if(iadv){
                          gx = *(sgrdf1 + k);
                          gy = *(sgrdf2 + k);

                          mulv(&ephi[k], &vx, gx);
                          mulv(&etht[k], &vy, gy);             
                          addv(&vx, &vy, &vel);

                          gx = dotp(&etan[k], &vel);
                          gy = dotp(&erad[k], &vel);

                       }

                       if(irel){
                          ang = sqrt(1.0 - vecg[k].x * vecg[k].x);
                          ang = acos(ang);
                          if(vecg[k].x < 0.0) ang = -ang;
                          sincos(ang, &s1, &c1);
                          vtmp.x = -c1;
                          vtmp.y = 0.0;
                          vtmp.z = s1;
                          mulv(&vtmp, &vtmp, pdirl);
                          uu -= dotp(&vtmp, &etnn[k]);
                          vv -= dotp(&vtmp, &erdd[k]);

                       }

                       if(iflx){
                          if(idevazm){
                             *(tmp1 + k) = vv;
                             *(tmp2 + k) = *(sadd + k);
                          }
                          else {
                             if(!irif){
                               if(lyes){
                                  *(*(savg1 + j) + k) += uu * *(sadd + k) * (ad + bd * *(sspres + k));
                                  *(*(savg2 + j) + k) += vv * *(sadd + k) * (ad + bd * *(sspres + k));
/*				  if(j == irnf - 1){
				    *(*(savg1 + j) + k) += uu * *(sadd + k) * ((1.0 - *(bb + j)) * *(sspres + k) - *(aa + j));
				    *(*(savg2 + j) + k) += vv * *(sadd + k) * ((1.0 - *(bb + j)) * *(sspres + k) - *(aa + j));
				  }
*/
                                  if(ikecb) *(*(savg3 + j) + k) += *(sadd + k) * *(sadd2 + k) * (ad + bd * *(sspres + k));
                                  else *(*(savg3 + j) + k) += *(sadd + k) * (ad + bd * *(sspres + k));
                               }
                               else {
                                  *(*(savg1 + j) + k) += uu * *(sadd + k);
                                  *(*(savg2 + j) + k) += vv * *(sadd + k);
                                  if(ikecb) *(*(savg3 + j) + k) += *(sadd + k) * *(sadd2 + k);
                                  else *(*(savg3 + j) + k) += *(sadd + k);
                               }
                             }
                             else {
                               ke = 0.5 * (uu * uu + vv * vv);
                               if(lyes){
                                  *(*(savg1 + j) + k) += uu * ke * (ad + bd * *(sspres + k));
                                  *(*(savg2 + j) + k) += vv * ke * (ad + bd * *(sspres + k));
                                  if(ikecb) *(*(savg3 + j) + k) += ke * *(sadd2 + k) * (ad + bd * *(sspres + k));
                                  else *(*(savg3 + j) + k) += ke * (ad + bd * *(sspres + k));

                               }
                               else {
                                  *(*(savg1 + j) + k) += uu * ke;
                                  *(*(savg2 + j) + k) += vv * ke;
                                  if(ikecb) *(*(savg3 + j) + k) += ke * *(sadd2 + k);
                                  else *(*(savg3 + j) + k) += ke;
                               }
                             }
                          }
                       }

                       else if(iadv){
                          if(lyes)
                             *(*(savg1 + j) + k) += uu * gx + vv * gy * (ad + bd * *(sspres + k));
                          else
                             *(*(savg1 + j) + k) += uu * gx + vv * gy;
                       }

                       *(tmpdat + k) = sqrt(uu * uu + vv * vv);
                       if(lyes){
                          *(*(aavg1 + j) + k) += uu * (ad + bd * *(sspres + k));
                          *(*(aavg2 + j) + k) += vv * (ad + bd * *(sspres + k));
                          *(*(aavg + j) + k) += *(tmpdat + k) * (ad + bd * *(sspres + k));
                       }
                       else {                       
                          *(*(aavg1 + j) + k) += uu;
                          *(*(aavg2 + j) + k) += vv;
                          *(*(aavg + j) + k) += *(tmpdat + k);
                       }
                       if(istd){
                          *(*(astd1 + j) + k) += uu * uu;
                          *(*(astd2 + j) + k) += vv * vv;
                          *(*(astd + j) + k) += uu * uu + vv * vv;
                       }

                   }

                   if(idevazm){
                      lavg1 = line_anomaly(tmp1, irnth, irr);
                      lavg2 = line_anomaly(tmp2, irnth, irr);
                      aravg = anomaly(sadd, coslat, irnr, irnth, irr);
                      if(lavg1 < ADD_CHECK && lavg2 <  ADD_CHECK && aravg < ADD_CHECK){
                         *(alv + j) += (lavg2 - aravg) * lavg1;
                         ++(*(nalv + j));
                      }
                      for(k=0; k < irdim; k++){
                          if(*(tmp1 + k) > ADD_CHECK || *(tmp2 + k) > ADD_CHECK) continue;
                          *(*(savg2 + j) + k) += *(tmp1 + k) * *(tmp2 + k);
                      }
                   }

                   radial_extreme(frmax, frmin, frange, tmpdat, slat, rthresh, irnth, irnr, *(rmposx + j), *(rmposn + j), *(rrpos + j));

                }

             }

          }

          atr2->trpt = atr2->trpt + frst;
          atr2->num = last - frst + 1;
   
          if(isprec){
            sumprc = 0.0;
            for(j=0; j < atr2->num; j++){
                 fp = atr2->trpt + j;
                 sumprc += fp->add_fld[ifsprec];
            }
            fprintf(fprc, "%e\n", 6.0 * sumprc);
            sumpp += 6.0 * sumprc;
          }

          if(ispd){

             for(j=0; j < atr2->num; j++){
                 fp = atr2->trpt + j;
                 pdi += pow(fp->add_fld[ifsp], 3.0);
             }
          }
       }

    }

    fclose(fdir);
    fclose(ftim);
    fclose(fdout);
    if(irad) {
       fclose(frmax); fclose(frmin); fclose(frange); free(tmpdat);
       for(i=0; i < irnf; i++) {free(*(rmposx + i)); free(*(rmposn + i)); free(*(rrpos + i));}
       free(rmposx); free(rmposn); free(rrpos);
    }

    if(ifcnt) {free(ifncnt); free(ifncntp);}

    if(fprc) fclose(fprc);

    strncpy(filout, filin, MAXCHR);
    strcat(filout, ".tcident");

    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    meantrd(fout, tracks, trnum, 's', gpr, ipr, alat, alng);

    fclose(fout);

    fout = fopen(hurr, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", hurr);
       exit(1);
    }

    meantrd(fout, tracks2, trnum, 's', gpr, ipr, alat, alng);

    fclose(fout);

    for(i=0; i < trnum; i++){
        atr = tracks + i;
        if(atr->num) ++ntrack;
        for(j=0; j < atr->num; j++)free((atr->trpt + j)->add_fld);
        free(atr->trpt);
    }
    free(tracks);

    if(fsfp) fclose(fsfp);

    if(irad){
       fclose(frad1);
       if(iwind) fclose(frad2);
       if(fadd) fclose(fadd);
       if(iadv) {fclose(fgrdx); fclose(fgrdy);}



       for(i=0; i < irnf; i++){
           for(j=0; j < irdim; j++) {
              if(*(*(ncount + i) + j)) {
                 *(*(aavg + i) + j) /= *(*(ncount + i) + j);
                 if(istd) *(*(astd + i) + j) /= *(*(ncount + i) + j);
              }
              else {
                 *(*(aavg + i) + j) = ADD_UNDEF;
                 if(istd) *(*(astd + i) + j) = ADD_UNDEF;
              }
           }

       }


       nll = irnf;

       if(lyes){
          printf("****WARNING****, all levels will be used in vertical integration.\n\n");
          nll1 = 1;
          nll2 = irnf;
       }
       else {
          printf("There are %d levels, specify the levels between which you want to restrict the vertical integration, l1, l2?\r\n"
                 "i.e. first NL levels will be used.                                                                           \n\n", irnf);
          scanf("%d %d", &nll1, &nll2);
          if((nll1 < 1 || nll1 > irnf) || (nll2 < 1 || nll2 > irnf)){
             printf("****ERROR****, incorrect levels identifiers specified, exiting.\n\n");
             exit(1);
          }
       }
       nll = nll2;

       if(iaddf){

/* perform azimuthal line integral of fluxes */

          sumarc = FPI2 * (1.0 - sin(*(slat + irr - 1) * FP_PI));

          for(i=0; i < irnf; i++){

             line_area_integ(linteg1[i] , rintg1 + i, aavg[i], arcll, arclt, coslat, sumarc, irnth, irr);

          }

/* vertical integration */

          if(levdat || lyes) {

             for(i=0; i < irnf; i++) *(lvtmp + i) = *(linteg1[i] + irr - 1);
             vlintg1 = vert_integ(lvtmp, vint, nll1, nll2, lyes);
             vaintg1 = vert_integ(rintg1, vint, nll1, nll2, lyes);

             printf("Vertical integrals; Line = %e; Area = %e\n\n", vlintg1, vaintg1);
             printf("Scaled vertical integrals; Line = %e; Area = %e\n\n", vlintg1 * lscl, vaintg1 * ascl);

          }

          if(*(rintg1 + nll - 1) < ADD_CHECK)
             printf("Scaled area averaged value at the last pressure level = %e\n", *(rintg1 + nll - 1) * asclr / GGRAV);
          else
             printf("****ERROR****, missing value for area averaged last pressure level.\n\n");

       }
       else if(yinteg){
/* perform azimuthal line integral of fluxes */

          sumarc = FPI2 * (1.0 - sin(*(slat + irr - 1) * FP_PI));

          for(i=0; i < irnf; i++){

             line_area_integ(linteg1[i] , rintg1 + i, aavg[i], arcll, arclt, coslat, sumarc, irnth, irr);

          }
          if(levdat || lyes) {

             for(i=0; i < irnf; i++) *(lvtmp + i) = *(linteg1[i] + irr - 1);
             vlintg1 = vert_integ(lvtmp, vint, nll1, nll2, lyes);
             vaintg1 = vert_integ(rintg1, vint, nll1, nll2, lyes);

             printf("Vertical integrals; Line = %e; Area = %e\n\n", vlintg1, vaintg1);
             printf("Scaled vertical integrals; Line = %e; Area = %e\n\n", vlintg1 * lscl, vaintg1 * asclr / GGRAV);
          }


       }

       if(iwind){
          for(i=0; i < irnf; i++){
             for(j=0; j < irdim; j++) {
                if(*(*(ncount + i) + j))  {
                   *(*(aavg1 + i) + j) /= *(*(ncount + i) + j);
                   *(*(aavg2 + i) + j) /= *(*(ncount + i) + j);
                   if(istd){
                      *(*(astd1 + i) + j) /= *(*(ncount + i) + j);
                      *(*(astd2 + i) + j) /= *(*(ncount + i) + j);
                   }
                   if(iflx){
                      *(*(savg1 + i) + j) /= *(*(ncount + i) + j);
                      *(*(savg2 + i) + j) /= *(*(ncount + i) + j);
                      *(*(savg3 + i) + j) /= *(*(ncount + i) + j);
                   }
                   else if(iadv)
                      *(*(savg1 + i) + j) /= *(*(ncount + i) + j);
                }
                else {
                   *(*(aavg1 + i) + j) = ADD_UNDEF;
                   *(*(aavg2 + i) + j) = ADD_UNDEF;
                   if(istd){
                      *(*(astd1 + i) + j) = ADD_UNDEF;
                      *(*(astd2 + i) + j) = ADD_UNDEF;
                   }
                   if(iflx){
                      *(*(savg1 + i) + j) = ADD_UNDEF;
                      *(*(savg2 + i) + j) = ADD_UNDEF;
                      *(*(savg3 + i) + j) = ADD_UNDEF;
                   }
                   else if(iadv)
                      *(*(savg1 + i) + j) = ADD_UNDEF;
                }
             }
          }

          if(iflx) {

/* perform azimuthal line and area integrals of fluxes */

             sumarc = FPI2 * (1.0 - sin(*(slat + irr - 1) * FP_PI));

             for(i=0; i < irnf; i++){

                line_area_integ(linteg1[i], rintg1 + i, savg1[i], arcll, arclt, coslat, sumarc, irnth, irr);
                line_area_integ(linteg2[i], rintg2 + i, savg2[i], arcll, arclt, coslat, sumarc, irnth, irr); 

/* area averaged KE */

                line_area_integ(linteg3[i], rintg3 + i, savg3[i], arcll, arclt, coslat, sumarc, irnth, irr);

/*                for(j=0; j < irdim; j++) *(sradf1 + j) = (float)(*(savg3[i] + j));
                *(vke + i) = anomaly(sradf1, coslat, irnr, irnth, irr);               */   

             }

             if(ikecb){
               if(*(rintg3 + nll - 1) < ADD_CHECK)
                  printf("Scaled area averaged value at the last pressure level = %e\n", *(rintg3 + nll - 1) * asclr / GGRAV);
               else
                  printf("****ERROR****, missing value for area averaged last pressure level.\n\n");
             }

             if(levdat || lyes) {

                for(i=0; i < irnf; i++) *(lvtmp + i) = *(linteg1[i] + irr - 1);
                vlintg1 = vert_integ(lvtmp, vint, nll1, nll2, lyes);
                vaintg1 = vert_integ(rintg1, vint, nll1, nll2, lyes);
                for(i=0; i < irnf; i++) *(lvtmp + i) = *(linteg2[i] + irr - 1);
                vlintg2 = vert_integ(lvtmp, vint, nll1, nll2, lyes);
                vaintg2 = vert_integ(rintg2, vint, nll1, nll2, lyes);


                for(i=0; i < irdim; i++){
                    for(j=0; j < irnf; j++){
                        *(lvtmp + j) = *(savg1[j] + i);
                    }
                    *(vsavg1 + i) = vert_integ(lvtmp, vint, nll1, nll2, lyes);
                    for(j=0; j < irnf; j++){
                        *(lvtmp + j) = *(savg2[j] + i);
                    }
                    *(vsavg2 + i) = vert_integ(lvtmp, vint, nll1, nll2, lyes);
                }

/*                vkeintg = vert_integ(vke, vint, nll); */

                vkeintg = vert_integ(rintg3, vint, nll1, nll2, lyes);

                printf("Area Integrated KE: ");
                for(i=0; i < irnf; i++){
                    printf("%6.4e ", *(rintg3 + i));
                }
                printf("\n\n");


                printf("Vertical integrals of tangential; Line = %e; Area = %e\n\n", vlintg1, vaintg1);
                printf("Scaled vertical integrals of tangential; Line = %e; Area = %e\n\n", vlintg1 * lscl, vaintg1 * ascl);
                printf("Vertical integrals of radial; Line = %e; Area = %e\n\n", vlintg2, vaintg2);
                printf("Scaled vertical integrals of radial; Line = %e; Area = %e\n\n", vlintg2 * lscl, vaintg2 * ascl);
                printf("Vertical integral of area averaged kinetic energy = %e\n\n", vkeintg * asclr / GGRAV);

                if(idevazm){
                   for(j=0; j < irnf; j++){
                       if(*(nalv + j) > 0) *(alv + j) /= *(nalv + j);
                   }
                   vphi = vert_integ(alv, vint, nll1, nll2, lyes);
                   printf("Vertical integral of (phil - phia)*vnl = %e\n\n", vphi*lscl);
                }

/* find radius of maximum flux */

                fmxmn = SMALL;
                for(i=0; i < irr; i++){
                   vlintg2 = 0.0;
                   for(j=0; j < irnf; j++) *(lvtmp + j) = *(linteg1[j] + i);
                   vlintg2 = vert_integ(lvtmp, vint, nll1, nll2, lyes);

                   if(vlintg2 > fmxmn){fmxmn = vlintg2; iflxrd = i;}
                }
                printf("Max flux = %e at radius = %f\n", fmxmn, *(slat + iflxrd));
             }

          }

          else if (iadv){

/* perform azimuthal line and area integrals of advection */

             sumarc = FPI2 * (1.0 - sin(*(slat + irr - 1) * FP_PI));

             for(i=0; i < irnf; i++)
                line_area_integ(linteg1[i], rintg1 + i, savg1[i], arcll, arclt, coslat, sumarc, irnth, irr);

             if(levdat || lyes) {


                for(i=0; i < irnf; i++) *(lvtmp + i) = *(linteg1[i] + irr - 1);
                vlintg1 = vert_integ(lvtmp, vint, nll1, nll2, lyes);
                vaintg1 = vert_integ(rintg1, vint, nll1, nll2, lyes);

                printf("Vertical integrals of Line = %e; Area = %e\n\n", vlintg1, vaintg1);
                printf("Scaled vertical integral of Area = %e\n\n", vaintg1 * scladv);

             }


          }

/* compute radial and tangential wind speeds from averaged U and V */

          if(!icomp_order && !idir){

             if(istd){
                printf("****WARNING****, no standard deviations for vector components calculated from composite first approach.\n\n");
             }

             for(i=0; i < irnf; i++) memcpy(*(aavg1 + i), *(aavg + i), irdim * sizeof(float));

             for(i=0; i < irnf; i++){

/* compute wind speed */

                 for(j=0; j < irdim; j++) {
                     uu = *(*(aavg1 + i) + j);
                     vv = *(*(aavg2 + i) + j);
                     *(*(aavg + i) + j) = sqrt(uu * uu + vv * vv);
                 }

             }

             for(i=0; i < irnr; i++){

                rad = FP_PI2 - *(slat + i) * FP_PI;
                if(rad < 0.0) rad = 0.0;
                sincos(rad, &s1, &c1);

                for(j=0; j < irnth; j++){

                    thet = *(slng + j) * FP_PI;
                    sincos(thet, &s2, &c2);

                    ptt.x = s1 * c2;
                    ptt.y = s1 * s2;
                    ptt.z = c1;

                    sq = sqrt(ptt.y * ptt.y + ptt.z * ptt.z);

                    erd.x = c1 * c2;
                    erd.y = c1 * s2;
                    erd.z = -s1;

                    etn.x = -s2;
                    etn.y = c2;
                    etn.z = 0;

                    for(k=0; k < irnf; k++){

                       uu = *(*(aavg1 + k) + i * irnth + j);
                       vv = *(*(aavg2 + k) + i * irnth + j);   

                       ww.x = -sq * vv;
                       ww.y = (uu * ptt.z + vv * ptt.x * ptt.y) / sq;
                       ww.z = (vv * ptt.x * ptt.z - uu * ptt.y) / sq;

                       *(*(aavg1 + k) + i * irnth + j) = dotp(&etn, &ww);
                       *(*(aavg2 + k) + i * irnth + j) = dotp(&erd, &ww);

                    }
 
                 }

             }

          }

       }


       if(istd){
          for(i=0; i < irnf; i++){
             for(j=0; j < irdim; j++) {
                 if(*(*(aavg + i) + j) < ADD_CHECK){
                   *(*(astd + i) + j) -= *(*(aavg + i) + j) * *(*(aavg + i) + j);
                   *(*(astd + i) + j) = (*(*(astd + i) + j) > 0.0) ? sqrt(*(*(astd + i) + j)) : 0.0; 
                 }
                 else *(*(astd + i) + j) = ADD_UNDEF;

                 if(iwind){
                    if(*(*(aavg1 + i) + j) < ADD_CHECK) {
                       *(*(astd1 + i) + j) -= *(*(aavg1 + i) + j) * *(*(aavg1 + i) + j);
                       *(*(astd1 + i) + j) = (*(*(astd1 + i) + j) > 0.0) ? sqrt(*(*(astd1 + i) + j)) : 0.0;
                    }
                    else *(*(astd1 + i) + j) = ADD_UNDEF;

                    if(*(*(aavg2 + i) + j) < ADD_CHECK) {
                       *(*(astd2 + i) + j) -= *(*(aavg2 + i) + j) * *(*(aavg2 + i) + j);
                       *(*(astd2 + i) + j) = (*(*(astd2 + i) + j) > 0.0) ? sqrt(*(*(astd2 + i) + j)) : 0.0;
                    }
                    else *(*(astd2 + i) + j) = ADD_UNDEF;
                 }
             }
          }          

       }

/* write to file */

       write_fields("reg_avg.dat", aavg, slng, slat, irnth, irnr, irnf);
       if(istd) write_fields("reg_std.dat", astd, slng, slat, irnth, irnr, irnf);

       if(iwind){
          write_fields("reg_avg_tan.dat", aavg1, slng, slat, irnth, irnr, irnf);
          write_fields("reg_avg_rad.dat", aavg2, slng, slat, irnth, irnr, irnf);
          if(istd && icomp_order){
             write_fields("reg_std_tan.dat", astd1, slng, slat, irnth, irnr, irnf);
             write_fields("reg_std_rad.dat", astd2, slng, slat, irnth, irnr, irnf);
          }
          if(iflx){
             write_fields("reg_avg_flux_tan.dat", savg1, slng, slat, irnth, irnr, irnf);
             write_fields("reg_avg_flux_rad.dat", savg2, slng, slat, irnth, irnr, irnf);
          }
          else if(iadv)
             write_fields("reg_avg_adv.dat", savg1, slng, slat, irnth, irnr, irnf);
       }

       printf("Write fields as NETCDF files, 'y' or 'n'\n\n");
       scanf("\n");
       if(getchar() == 'y'){
          levdat = write_fields_netcdf("reg_avg.nc", aavg, slng, slat, levdat, irnth, irnr, irnf, iext, imxmn, iwind, 0, ctime);
          if(istd) write_fields_netcdf("reg_std.nc", astd, slng, slat, levdat, irnth, irnr, irnf, iext, imxmn, iwind, -1, ctime);

          if(iaddf){
             write_intflux_netcdf("int_comb.nc", linteg1, slat, levdat, irr, irnf, imxmn, iext, 0, ctime);
             write_profile_netcdf("pro_comb.nc", rintg1, levdat, irnf, iext, imxmn, ctime);
          }

          if(iwind){
             write_fields_netcdf("reg_avg_tan.nc", aavg1, slng, slat, levdat, irnth, irnr, irnf, iext, imxmn, iwind, 1, ctime);
             write_fields_netcdf("reg_avg_rad.nc", aavg2, slng, slat, levdat, irnth, irnr, irnf, iext, imxmn, iwind, 2, ctime);
             if(istd && icomp_order){
                write_fields_netcdf("reg_std_tan.nc", astd1, slng, slat, levdat, irnth, irnr, irnf, iext, imxmn, iwind, 1, ctime);
                write_fields_netcdf("reg_std_rad.nc", astd2, slng, slat, levdat, irnth, irnr, irnf, iext, imxmn, iwind, 2, ctime);
             }
             if(iflx){

                write_fields_netcdf("reg_avg_flux_tan.nc", savg1, slng, slat, levdat, irnth, irr, irnf, iext, imxmn, iwind, 3, ctime);
                write_fields_netcdf("reg_avg_flux_rad.nc", savg2, slng, slat, levdat, irnth, irr, irnf, iext, imxmn, iwind, 4, ctime);
                write_intflux_netcdf("int_flux_tan.nc", linteg1, slat, levdat, irr, irnf, imxmn, iext, 1, ctime);
                write_intflux_netcdf("int_flux_rad.nc", linteg2, slat, levdat, irr, irnf, imxmn, iext, 2, ctime);

                for(i=0; i < irnf; i++) {
                   *(rintg1 + i) = *(linteg1[i] + irr - 1);
                   *(rintg2 + i) = *(linteg2[i] + irr - 1);
                }

vaintg2 = vert_integ(rintg2, vint, nll1, nll2, lyes);
printf("%e %e\n", lscl, vaintg2 * lscl);

                write_profile_netcdf("pro_flux_tan.nc", rintg1, levdat, irnf, iext, imxmn, ctime);
                write_profile_netcdf("pro_flux_rad.nc", rintg2, levdat, irnf, iext, imxmn, ctime);

                if(levdat || lyes){
                   write_fields_netcdf("reg_vrtavg_flux_tan.nc", &vsavg1, slng, slat, levdat, irnth, irr, 1, iext, imxmn, iwind, 3, ctime);
                   write_fields_netcdf("reg_vrtavg_flux_rad.nc", &vsavg2, slng, slat, levdat, irnth, irr, 1, iext, imxmn, iwind, 4, ctime);
                } 

             }
             else if(iadv){
                write_fields_netcdf("reg_avg_adv.nc", savg1, slng, slat, levdat, irnth, irr, irnf, iext, imxmn, iwind, 5, ctime);
                write_intflux_netcdf("int_adv.nc", linteg1, slat, levdat, irr, irnf, imxmn, iext, 3, ctime);
                write_profile_netcdf("pro_adv.nc", rintg1, levdat, irnf, iext, imxmn, ctime);

             }
          }

       }

       free(aa); free(bb);

       free(vsavg1); free(vsavg2);
       free(sradf1); free(sradf2); 
       free(tmp1); free(tmp2);
       free(sgrdf1); free(sgrdf2);
       free(slng);
       free(slat);
       free(coslat);
       if(lyes) free(sspres);
       if(iwind) {
          free(slng2); free(slat2);
          for(i=0; i < irnf; i++) free(*(aavg1 + i));
          free(aavg1);
          for(i=0; i < irnf; i++) free(*(aavg2 + i));
          free(aavg2);

          if(istd){
             for(i=0; i < irnf; i++) free(*(astd1 + i));
             free(astd1);
             for(i=0; i < irnf; i++) free(*(astd2 + i));
             free(astd2);
          }

          if(icomp_order){

             if(iflx){
                for(i=0; i < irnf; i++) {
                    free(*(savg1 + i)); free(*(savg2 + i)); free(*(savg3 + i)); 
                    free(linteg1[i]); free(linteg2[i]); free(linteg3[3]);
                }
                free(savg1);
                free(savg2);
                free(savg3);
                free(linteg1); free(linteg2); free(linteg3);
                free(rintg1); free(rintg2); free(rintg3);
                free(vke);
                free(alv);
                free(nalv);
                free(tmp1);
                free(tmp2);
                free(sadd2);
             }
             else if(iadv){
                for(i=0; i < irnf; i++) free(linteg1[i]);
                free(linteg1);
                free(rintg1);
                free(sgrdf1);
                free(sgrdf2);
             }
             free(lvtmp);
          }

          free(vecg); free(newv);
          free(etht); free(ephi);
          free(erad); free(etan);
          free(erdd); free(etnn);
       }
       free(sadd);

       if(levdat){
          free(levdat);
          free(vint);
       }

       for(i=0; i < irnf; i++) {free(*(aavg + i)); free(*(ncount + i));}
       free(aavg);
       if(istd){
          for(i=0; i < irnf; i++) free(*(astd + i));
          free(astd);
       }
       free(ncount);
    }


    printf("Number of tracks is %d\n\n", ntrack);
    printf("PDI is %e\n\n", pdi);
    printf("Precip (avg) is %e\n\n", sumpp / ntrack);

    return 0;

}


/* write fields in ascii/binary format */

void write_fields(char *name, double **fld, float *slng, float *slat, int irnth, int irnr, int irnf)
{
    int i, j;
    int ii, irdim;

    float *fldt=NULL;

    FILE *fout=NULL;

    irdim = irnth * irnr;

    fldt = (float *)calloc(irdim, sizeof(float));
    mem_er((fldt == NULL) ? 0 : 1, irdim*sizeof(float));

    fout = fopen(name, "w");

    fprintf(fout, "%d %d %d\n", irnth, irnr, irnf);
    ii = 0;
    for(i=0; i < irnth; i++){
       fprintf(fout, "%8.4f ", *(slng + i));
       ++ii;
       if(ii == 10){fprintf(fout, "\n"); ii = 0;}
    }
    if(ii) fprintf(fout, "\n");

    ii = 0;
    for(i=0; i < irnr; i++){
       fprintf(fout, "%8.4f ", *(slat + i));
       ++ii;
       if(ii == 10){fprintf(fout, "\n"); ii = 0;}
    }
    if(ii) fprintf(fout, "\n");

    for(i=0; i < irnf; i++){
       for(j=0; j < irdim; j++) *(fldt + j) = *(*(fld + i) + j);
       fprintf(fout, "FIELD %d\n", i+1);
       fwrite(fldt, irdim * sizeof(float), 1, fout);
       fprintf(fout, "\n");
    }

    fclose(fout);

    free(fldt);

    return;

}


/* determine the deviation from the areal mean */

double anomaly(float *sradf, double *coslat, int irnr, int irnth, int iran)
{
     int i, j;
     int irdim=0;

     double sum=0.0, sumcl=0.0; 

     irdim = irnr * irnth;

     for(i=0; i < iran; i++) {
        for(j=0; j < irnth; j++){
           if(*(sradf + i * irnth + j) > ADD_CHECK) continue;
           sum += *(coslat + i) * *(sradf + i * irnth + j);
           sumcl += *(coslat + i);
        }
     }

     if(sumcl > 0.0)
        sum /= sumcl;
     else
        sum = ADD_UNDEF;

     for(i=0; i < irdim; i++) {
        if(*(sradf + i) > ADD_CHECK) continue;
        *(sradf + i) -= sum;
     }

    return sum;
}

/* calculate deviation from the mean along the azimuthal direction */

double line_anomaly(float *tmp, int irnth, int irr)
{
    int i, j;

    int nl=0;

    float ff;
    double linsum=0.0;

    for(i=0; i < irr; i++){
       linsum = 0.0;
       nl = 0;
       for(j=0; j < irnth; j++){
           ff = *(tmp + i * irnth + j);
           if(ff > ADD_CHECK) continue;
           linsum += ff;
           ++nl;
       }

       if(nl > 0){
         linsum /= nl;
         for(j=0; j < irnth; j++) {
             ff = *(tmp + i * irnth + j);
             if(ff > ADD_CHECK) continue;
             *(tmp + i * irnth + j) -= linsum;
         }
       }
       else{
         linsum = ADD_UNDEF;
         for(j=0; j < irnth; j++) *(tmp + i * irnth + j) = ADD_UNDEF;
       }
    }

    return linsum;

}

/* determine the extreme values of the composite along radii */

void radial_extreme(FILE *frmax, FILE *frmin, FILE * frange, float *sradf, float *slat, float rthresh, int irnth, int irnr, off_t *rmposx, off_t *rmposn, off_t *rrpos)
{

    int i, j;
    int imn=0, imx=0;

    float rmax=0.0, rmin=0.0;
    float rmps1=0.0, rmps2=0.0;
    float rv1=0.0, rv2=0.0;
    float rmx=0.0, rmn=0.0;

    for(i=0; i < irnth; i++){
       rmax = SMALL;
       rmin = LARGE;
       rmps1 = *slat;
       rmps2 = *slat;
       imn = imx = 0;
       for(j=0; j < irnr; j++){
           rv1 = *(sradf + j * irnth + i);
           rv2 = *(sradf + (irnr - j - 1) * irnth + i);
           if(rv1 > rmax){rmax = rv1; rmps1 = *(slat + j);}
           if(rv1 < rmin){rmin = rv1; rmps2 = *(slat + j);}
           if(rv1 >= rthresh && !imn) {rmn = *(slat + j); imn = 1;}
           if(rv2 >= rthresh && !imx) {rmx = *(slat + irnr - j - 1); imx = 1;}
       }

       fseeko(frmax, *(rmposx + i), SEEK_SET);
       fprintf(frmax, "%8.4f %e & ", rmps1, rmax);
       *(rmposx + i) = ftell(frmax);

       fseeko(frmin, *(rmposn + i), SEEK_SET);
       fprintf(frmin, "%8.4f %e & ", rmps2, rmin);
       *(rmposn + i) = ftell(frmin);

       fseeko(frange, *(rrpos + i), SEEK_SET);
       fprintf(frange, "%8.4f %8.4f & ", rmn, rmx);
       *(rrpos + i) = ftell(frange);
    }

    return;
}

/* read level data from stdin */

float *leveldata(int nlev)
{
    int i;
    float *levdat=NULL;

    levdat = (float *)calloc(nlev, sizeof(float));
    mem_er((levdat == NULL) ? 0 : 1, nlev * sizeof(float));
    printf("Specify the level values in the correct order.\n\n");
    for(i=0; i < nlev; i++){
        printf("Value for level %d: ", i+1);
        scanf("%f", levdat + i);
        printf("\n");
    }

    return levdat;
}

/* determine rotation angle for alinging grid with system direction */

double dirangle(struct tot_tr *altr, double *cn, double *sn, int pt_id, int ndsmth)
{
   int i=0;
   int n2, st=0, en=0;
   int narot=0;

   double aa=0.0, arot=0.0, norm=0.0;
   double xx, yy;
   double s1, s2, c1, c2;

   struct fet_pt_tr *atr=NULL;

   VEC vsdir={0.0, 0.0, 0.0};
   VEC pt1, pt2, tvec, vtmp;
   VEC vt1, vt2, vt3, vn;

   n2 = ndsmth / 2;

   atr = altr->trpt + pt_id;

   if(altr->num <= 1){
      printf("****ERROR****, insufficient number of points in this track to determine direction.\n\n");
      exit(1);
   }

   if(ndsmth % 2) {st = pt_id - n2; en = pt_id + n2;}
   else {st = pt_id - n2 - 1; en = pt_id + n2;}

   if(st < 0) st = 0;
   if(en > altr->num - 2) en = altr->num - 2;
   if(st > en) st = en;

   for(i = st; i <= en; i++){
       atr = altr->trpt + i;
       xx = atr->xf * FP_PI;
       yy = FP_PI2 - atr->yf * FP_PI;
       if(yy < 0.) yy = 0.0;
       sincos(xx, &s1, &c1);
       sincos(yy, &s2, &c2);
       vsdir.x = -c1 * c2;
       vsdir.y = -s1 * c2;
       vsdir.z = s2;

       pt1.x = atr->pp[0];
       pt1.y = atr->pp[1];
       pt1.z = atr->pp[2];
       pt2.x = (atr+1)->pp[0];
       pt2.y = (atr+1)->pp[1];
       pt2.z = (atr+1)->pp[2];
       crosp(&pt1, &pt2, &vtmp);
       crosp(&vtmp, &pt1, &tvec);
       norm = sqrt(dotp(&tvec, &tvec));
       normv(&tvec, norm);

/* compute rotation angle */

       aa = acos(dotp(&vsdir, &tvec));

/* compute which side of tvec that vsdir is for anticlockwise rotation */

       norm = dotp(&vtmp, &vsdir);
       if(norm < 0.0) aa *= -1.0;
       sincos(aa, &s1, &c1);
/* rotate standard direction for check */
       norm = dotp(&vsdir, &pt1);
       crosp(&vsdir, &pt1, &vt1);
       mulv(&vt1, &vt1, s1);
       mulv(&pt1, &vt2, ((1.0 - c1) * norm));
       mulv(&vsdir, &vt3, c1);
       addv(&vt1, &vt2, &vn);
       addv(&vn, &vt3, &vn);

       if(fabs(dotp(&vn, &tvec) - 1.0) > 1.0e-6){
         printf("****ERROR****, orientating radial grid to storm direction incorrect, dotp=%e > 1.0e-6.\n\n", fabs(dotp(&vn, &tvec) - 1.0));
         exit(1);
       }

       arot += aa;
       ++narot;

   }

   arot /= narot;

   sincos(arot, sn, cn);

   return arot;
}

/* determine system propogation speed over three time steps */

double prop_speed(struct tot_tr *altr, int ifld)
{

    double spp=0.0;

    VEC ptt, pt0, pt1;

    struct fet_pt_tr *fpt=NULL, *fp0=NULL, *fp1=NULL;
    
    fpt = altr->trpt + ifld;
    ptt.x = fpt->pp[0];
    ptt.y = fpt->pp[1];
    ptt.z = fpt->pp[2];

    if(ifld+1 < altr->num && ifld - 1 >= 0){
      fp0 = fpt - 1;
      fp1 = fpt + 1;
      pt0.x = fp0->pp[0];
      pt0.y = fp0->pp[1];
      pt0.z = fp0->pp[2];
      pt1.x = fp1->pp[0];
      pt1.y = fp1->pp[1];
      pt1.z = fp1->pp[2];
      spp += acos(dotp(&ptt, &pt0));
      spp += acos(dotp(&ptt, &pt1));
      spp *= 0.5;
    }
    else if(ifld+1 < altr->num){
      fp1 = fpt + 1;
      pt1.x = fp1->pp[0];
      pt1.y = fp1->pp[1];
      pt1.z = fp1->pp[2];
      spp = acos(dotp(&ptt, &pt1));
    }
    else if(ifld - 1 >= 0){
      fp0 = fpt - 1;
      pt0.x = fp0->pp[0];
      pt0.y = fp0->pp[1];
      pt0.z = fp0->pp[2];
      spp = acos(dotp(&ptt, &pt0));
    }
    else {
      printf("****ERROR****, track point number inconsistency, exiting.\n\n");
      exit(1);
    }

    return spp;
}

/* perform line and area integrals using trapazoidal rule */

void line_area_integ(double *linteg, double *rintg, double *savg, double arcll, double arclt, double *coslat, double sumarc, int irnth, int irr)
{
    int j, k;
    int ilg=0, imis=0;
    
    double arcl=0.0;
    double val=0.0;
    double sumarea=0.0;
    double *sumarcl=NULL;

    sumarcl = (double *)calloc(irr, sizeof(double));
    mem_er((sumarcl == NULL) ? 0 : 1, irr*sizeof(double));

    *rintg = 0.0;
    for(j=0; j < irr; j++){

        arcl = *(coslat + j) * arcll;

        *(linteg + j) = 0.0;
        *(sumarcl + j)=0.0;
        ilg = 0;

        for(k=0; k < irnth; k++) {
           val = *(savg + j * irnth + k);
           if(val > ADD_CHECK) {imis = 1; continue;}
           else {
              ilg = 1;
              *(linteg + j) += val;
              *(sumarcl + j) += arcl;
           }
        }

        if(ilg)
           *(linteg + j) *= arcl;
        else
           *(linteg + j) = ADD_UNDEF;

    }

/* radial integral of fluxes */

    if(imis){

       for(j=1; j < irr-1; j++) {
           if(*(linteg + j) < ADD_CHECK){
              sumarea += *(sumarcl + j) * arclt;
              *rintg += *(linteg + j);
           }

       }

       *rintg *= arclt;

       if(*linteg < ADD_CHECK && *(linteg + irr - 1) < ADD_CHECK){
          *rintg += 0.5 * arclt * (*linteg + *(linteg + irr - 1));
          sumarea += 0.5 * arclt * (*sumarcl + *(sumarcl + irr - 1));
       }
       else if(*linteg < ADD_CHECK){
          *rintg += 0.5 * arclt * *linteg;
          sumarea += 0.5 * arclt * *sumarcl;
       }
       else if(*(linteg + irr - 1) < ADD_CHECK){
          *rintg += 0.5 * arclt * *(linteg + irr - 1);
          sumarea += 0.5 * arclt * *(sumarcl + irr - 1);
       }

       if(sumarea > 0.0)
          *rintg /= sumarea;
       else
          *rintg = ADD_UNDEF;

       for(j=0; j < irr; j++){

          if(*(linteg + j) < ADD_CHECK) 
             *(linteg + j) /= *(sumarcl + j);
          else
             *(linteg + j) = ADD_UNDEF;
       }

    }

    else{

       for(j=1; j < irr-1; j++)
           *rintg += *(linteg + j);

       *rintg *= arclt;
       *rintg += 0.5 * arclt * (*linteg + *(linteg + irr - 1));

       *rintg /= sumarc;

       for(j=0; j < irr; j++){
           arcl = *(coslat + j) * FPI2;
          *(linteg + j) /= arcl;
       }

    }

    free(sumarcl);

    return;
}


/* open radial data file for read */

FILE *open_radial_file(char *filrad, int *irtrn, long int *iptnum, int *irnth, int *irnr, int *irnf, int *idir, int *ndsmth, int *ifcnt, int *itpadd, int *iprojd)
{

    int irtrn2=0;
    int irnth2=0, irnr2=0, irnf2=0;
    int idir2=0, ndsmth2=0;
    int ifcnt2=0, itpadd2=0, iprojd2=0;

    long int iptnum2=0;

    static int ifrst=0;

    FILE *frad=NULL;
    char line[MAXCHR];

    frad = fopen(filrad, "r");
    if(!frad){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filrad);
       exit(1);
    }
    fgets(line, MAXCHR, frad);

    if(!ifrst){
       sscanf(line, "%d %ld %d %d %d %d %d %d %d %d", irtrn, iptnum, irnth, irnr, irnf, idir, ndsmth, ifcnt, itpadd, iprojd);
       ifrst = 1;
    }
    else{
       sscanf(line, "%d %ld %d %d %d %d %d %d %d %d", &irtrn2, &iptnum2, &irnth2, &irnr2, &irnf2, &idir2, &ndsmth2, &ifcnt2, &itpadd2, &iprojd2);

       if(*irtrn != irtrn2 || *iptnum != iptnum2 || *irnth != irnth2 ||
          *irnr != irnr2   || *irnf != irnf2 || *idir != idir2  || *ndsmth != ndsmth2 ||
          *ifcnt != ifcnt2 || *itpadd != itpadd2 || *iprojd != iprojd2                   ){
           printf("****ERROR****, incompatability between radial field files.\n\n");
           exit(1);
       }
    }

    return frad;
}

/* calculate vertical integral */

double vert_integ(double *dat, double *vlev, int nl1, int nl2, int lyes)
{
    int i=0;

    double vint=0.0;

    for(i=nl1-1; i < nl2; i++){
        if(*(dat + i) > ADD_CHECK){
           printf("****ERROR****, a level has a missing value.\n\n");
           vint = ADD_UNDEF;
           return vint;
        }
    }

    if(!lyes){
       for(i=nl1-1; i < nl2-1; i++){
           vint += 0.5 * (*(dat + i) + *(dat + i + 1)) * *(vlev + i); 
       }
    }
    else {
       for(i=nl1-1; i < nl2; i++) vint += *(dat + i);
    }

    return vint;
}

/* function sets up the rotation matrix */

void setup_rotmatrix(double s1, double c1, double s2, double c2, VEC *rot)
{

    rot[0].x = c1 * c2;
    rot[0].y = -s1;
    rot[0].z = c1 * s2;
    rot[1].x = s1 * c2;
    rot[1].y = c1;
    rot[1].z = s1 * s2;
    rot[2].x = -s2;
    rot[2].z = c2;

    return;
}

/* setup local radial vectors */

void setup_radvec(VEC *rot, VEC *vecg, VEC *pt1, VEC *newv, VEC *etht, VEC *ephi, VEC *erad, VEC *etan, double ss1, double cc1, int irdim, int idir)
{

    int i;

    VEC vt1, vt2, vt3, vn;

    double norm=0.0;


    for(i=0; i < irdim; i++){

        newv[i].x = dotp(&rot[0], &vecg[i]);
        newv[i].y = dotp(&rot[1], &vecg[i]);
        newv[i].z = dotp(&rot[2], &vecg[i]);

        if(idir){
           norm = dotp(&newv[i], pt1);
           crosp(&newv[i], pt1, &vt1);
           mulv(&vt1, &vt1, ss1);
           mulv(pt1, &vt2, ((1.0 - cc1) * norm));
           mulv(&newv[i], &vt3, cc1);
           addv(&vt1, &vt2, &vn);
           addv(&vn, &vt3, &newv[i]);
        }

        ephi[i].x = -(newv[i].y);
        ephi[i].y = newv[i].x;
        ephi[i].z = 0.0;
        norm = sqrt(dotp(&ephi[i], &ephi[i]));
        normv(&ephi[i], norm);
                   
        crosp(&ephi[i], &newv[i], &etht[i]);
        norm = sqrt(dotp(&etht[i], &etht[i]));
        normv(&etht[i], norm);

        crosp(pt1, &newv[i], &etan[i]);
        norm = sqrt(dotp(&etan[i], &etan[i]));
        normv(&etan[i], norm);

        crosp(&etan[i], &newv[i], &erad[i]);
        norm = sqrt(dotp(&erad[i], &erad[i]));
        normv(&erad[i], norm);
    }

    return;
} 

/* if a missing point is found find equivelent point at previous levels */

void find_last(struct fet_pt_tr *atr, struct feature_pts *fpt, int *ifncnt, int *ifncntp, int ic, int itpadd)
{

   int i, j, ii;
   int nmnf=0;
   int nxf=0;

   if(*(ifncnt + ic)){
     (fpt->x).xy = *(atr->add_fld + *(ifncntp + ic));
     (fpt->y).xy = *(atr->add_fld + *(ifncntp + ic) + 1);

     if((fpt->x).xy > ADD_CHECK){

       if(itpadd && ic > 0){
          nmnf = *(ifncntp + ic);
          ii = *(ifncnt + ic) - 2;

          nxf = 0;
          for(i=ii; i >= 0; i--){
              if(*(nfwpos + i)) {
                 nmnf -= 3;
                 nxf = 0;
                 for(j=ic-1; j>=0; j--) {
                    if(i+1 == *(ifncnt + j)) nxf = 1;
                 }

                 if(*(atr->add_fld + nmnf) > ADD_CHECK || !nxf) continue;

                 break;
              }
              else --nmnf;            
          }

          if(nxf){
            (fpt->x).xy = *(atr->add_fld + nmnf);
            (fpt->y).xy = *(atr->add_fld + nmnf + 1);
          }
          else {
            (fpt->x).xy = atr->xf;
            (fpt->y).xy = atr->yf;
          }

        }
        else {
          (fpt->x).xy = atr->xf;
          (fpt->y).xy = atr->yf;
        }
     }
   }
   else {
     (fpt->x).xy = atr->xf;
     (fpt->y).xy = atr->yf;
   }

   return;
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

/* function to read model level A's and B's */

void read_model_levs(char *filnm, float **aa, float **bb, int *nl)
{

    int i=0;

    FILE *flev=NULL;

    flev = fopen(filnm, "r");
    if(!flev){
      printf("****ERROR****, unable to open file %s\n\n", filnm);
      exit(1);
    }

    fscanf(flev, "%d", nl);

    *aa = (float *)calloc(*nl, sizeof(float));
    mem_er((*aa == NULL) ? 0 : 1, *nl * sizeof(float));

    *bb = (float *)calloc(*nl, sizeof(float));
    mem_er((*bb == NULL) ? 0 : 1, *nl * sizeof(float));

    for(i=0; i < *nl; i++){
        fscanf(flev, "%f %f", *aa + i, *bb + i);
    }

/* compute full level coefficients */

    for(i=0; i < *nl - 1; i++){
        *(*aa + i) = 0.5 * (*(*aa + i) + *(*aa + i + 1));
        *(*bb + i) = 0.5 * (*(*bb + i) + *(*bb + i + 1));
    } 


    fclose(flev);


}

/* function to return pressure level difference on model levels in terms of A's and B's */

void model_level_diff(float *a, float *b, float *aa, float *bb, int nl, int irnf)
{
    if(!nl){
       *a = 0.5 * (*(aa + 1) - *aa);
       *b = 0.5 * (*(bb + 1) - *bb);
    }
    else if(nl == irnf - 1){
       *a = 0.5 * (*(aa + irnf - 1) - *(aa + irnf - 2));
       *b = 0.5 * (*(bb + irnf - 1) - *(bb + irnf - 2));
    }
    else {
       *a = 0.5 * (*(aa + nl + 1) - *(aa + nl - 1));
       *b = 0.5 * (*(bb + nl + 1) - *(bb + nl - 1));
    }

    return;
}





