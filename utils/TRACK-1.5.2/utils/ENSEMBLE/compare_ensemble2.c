#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"

#define TOLMATCH    2.0
#define DISTLARGE  180.0
#define TOLNUM     0.6
#define MAXCHR     50
#define FILENAME    "match_"
#define ID1        -10.0
#define ID2        10.0
#define NID        40
#define IDD1       0.
#define IDD2       2.
#define NIDD       10
#define TPBW       5.0
#define PERC       0.6       /* overlap with region */
#define TREG       0.5      /* overlap with time   */

#define MATCH_STR      0
#define TRACK_MAX_OUT  0
#define MISS       ADD_CHECK
#define IDIFF          1


/* program to compare two track ensembles */

typedef struct _reg{
    float rlng1;
    float rlng2;
    float rlat1;
    float rlat2;
    long int rtim1;
    long int rtim2;
    int ig;
    int reg_type;
} REG;                      /* region structure */

typedef struct _tim{
    long int stim;     /* start time */
    int tstp;          /* current time step */
    int tsmp;          /* new time sampling */
    int iosmp;         /* output sampled file switch */
    int itsmp;         /* use time sampling switch */
} TSAMP;               /* time sampling structure  */


struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
int toverlap(struct tot_tr * , struct tot_tr * , long int * , long int * );
float trdist(struct tot_tr * , struct tot_tr * , float , long int , long int , int * , int * , int * , int , int , int , int , int , int , int );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
float phpstr(int , struct fet_pt_tr * );
int inregion(struct tot_tr * , REG * , int , int , int , int );
float intensity(struct tot_tr * , REG * , int , int , int , int , int , int * , int , int , int );
void orog_nearest(struct fet_pt_tr * , float * , float * , int , int , int * , int * );
int region(struct fet_pt_tr * , REG * , long int , int , int);
struct tot_tr *time_sample(struct tot_tr * , int trnum, TSAMP * );
long int new_time(long int , int );

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int main(int argc, char **argv)
{

    int i,j;
    int nump=0, numpm=0;
    int totnum=0;
    int jtrmin=0;
    int nmat=0;

    int *nm1=NULL;

    int trnum1, trnum2;
    int gpr1=0, ipr1=0, gpr2=0, ipr2=0;
    int nbin=0;
    int iper_num_t;
    int str_type=0, sep_type=0, sep_ws=0;
    int nmbin=0;
    int nid=NID, nidd=NIDD;
    int in1=0, in2=0, it1=0, it2=0;
    int ior=0, iorr=0;
    int ftype='s';
    int vtype=0;
    int irord=0;
    int nmtch=0;
    int numst=1;

    int inreg1=0, inreg2=0, inreg=0;
    int ireg=0;
    int idisty=0;
    int imax=0, imax1=0, imax2=0;
    int nfl1=0, nfl2=0;
    int infl1=0, infl2=0;
    int iff1=0, iff2=0;
    int ifpos1=-1, ifpos2=-1;
    int upos1=0, upos2=0;
    int iprint=0, nint1=0, nint2=0;

    int ilngo=0, ilato=0;
    int ilng=0, ilat=0;
    int istf=0, instr=0, instr1=0, instr2=0;
    int loc1=0, loc2=0, lloc1=0, lloc2=0;
    
    int ioff=0;

    int nff1=0, nff2=0;
    int nfld1=0, nfld2=0;
    int *nfwpos1=NULL, *nfwpos2=NULL;

    long int is1=0, is2=0;
    
    TSAMP tsamp1={0, 0, 0, 0, 0}, tsamp2={0, 0, 0, 0, 0};

    char filout[MAXCHR], ctype[2];
    char buff[MAXCHR];

/* Histograms:- 1st match, 1st no match, 2nd match, 2nd no match */

    float *hist[4]={NULL, NULL, NULL, NULL};
    float *nptm=NULL;
    float *npid=NULL;
    float *npdd=NULL;
    float **ntpdt=NULL;

    float id1=ID1, id2=ID2;
    float idd1=IDD1, idd2=IDD2;

    float alat1=0.0, alng1=0.0, alat2=0.0, alng2=0.0;
    float dist=0.0;
    float distmin, distm;
    float tolnum=TOLNUM, tolmatch=TOLMATCH;
    float trpm1=0.0, trpm2=0.0;

    float frat=0.0, ffrat=0.0;
    float ditty=0.0;
    float t1=0.0, t2=0.0, dt=0.0, did=0.0, ddd=0.0;
    float tpbw=TPBW;

    float *temp=NULL;
    float *temp2=NULL;
    float *temp3=NULL;
    float *temp4=NULL;

    float str1=0.0, str2=0.0, str, strm;
    float ss1, ss2;
    float stf1=0.0, stf2=0.0;
    float ntrm=0.0;
    float nmatch=0.0;
    float nomch1=0.0, nomch2=0.0;

    float xpos1=0.0, xpos2=0.0, ypos1=0.0, ypos2=0.0;

    float sum_per_t;
    float st1=0.0, st2=0.0;

    float forog=0.0;
    float *orogm=NULL;
    float *olng=NULL, *olat=NULL;

    REG reg={0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0};

    FILE *fin1=NULL, *fin2=NULL;
    FILE *fout=NULL;

    FILE *fcmp=NULL;

    FILE *ftpdt=NULL;

    FILE *fstr=NULL, *fdiff=NULL;

    struct tot_tr *tr1=NULL, *tr2=NULL;
    struct tot_tr *tm1=NULL, *tm2=NULL;
    struct tot_tr *atr1=NULL, *atr2=NULL;
    
    struct tot_tr *trsmp1=NULL, *trsmp2=NULL;
    struct tot_tr *tmsmp1=NULL, *tmsmp2=NULL;
    struct tot_tr *trtmp=NULL;
    
    struct tot_tr *ntm1=NULL, *ntm2=NULL;

    struct fet_pt_tr *fp1=NULL, *fp2=NULL, *fptmp=NULL;

    if(argc < 13){
       printf("Usage: ensemble [track file1] [track file2] [temp1] [temp2] [nbin]                              \r\n"
              "                [mean (0), max (1) or min (2) intensity : offset]                               \r\n"
	      "                [Mean (0)  or Min (1) seperation : use whole track (0) or start of track (1)]   \r\n"
              "                [format of main distributions #/month (0), frequency (1), pdf (2)]              \r\n"
              "                [Add. field Id. 1 : use location (0 or 1) : use time sampling]                  \r\n"
	      "                [Add. field Id. 2 : use location : use time sampling]                           \r\n"
	      "                [Ftype, 's' or 'vn', n=0 (speed), n=1 (tendency), n=2 (grwth)]                  \r\n"
	      "                [Region matching type, '-1' (no region matching), '0' (track in region),        \r\n"
	      "                '1' (extrema in region), '2' (60%% of points in region),'3' (extrema of region)]\r\n"
	      "                [Re-order tracks, '1' for yes or '0' for no] [Dist. Thresh. (optional)]         \r\n"
	      "                [Time. Thresh. (optional)][Intensity Range, I1, I2 (optional)]                  \n\n");
       exit(1);
    }

    sscanf(argv[3], "%f", &t1);
    sscanf(argv[4], "%f", &t2);
    sscanf(argv[5], "%d", &nbin);
    sscanf(argv[6], "%d:%d", &str_type, &ioff);
    sscanf(argv[7], "%d:%d", &sep_type, &sep_ws);
    sscanf(argv[8], "%d", &idisty);
    sscanf(argv[9], "%d:%d:%d", &infl1, &upos1, &(tsamp1.itsmp));
    sscanf(argv[10], "%d:%d:%d", &infl2, &upos2, &(tsamp2.itsmp));
    sscanf(argv[11], "%s", ctype); 
    ftype = ctype[0];
    sscanf(argv[12], "%d", &(reg.reg_type));

    if(argc > 13){
       if(argc >= 14){
          sscanf(argv[13], "%d", &irord);
       }
        
       if(argc > 14 && (argc != 16 && argc != 18)){
          printf("****ERROR****, in optional command line arguments.\n\n");
	  exit(1);
       }
       if(argc >= 16){
          sscanf(argv[14], "%f", &tolmatch);
          sscanf(argv[15], "%f", &tolnum);
       }
       if(argc == 18){
          sscanf(argv[16], "%f", &st1);
          sscanf(argv[17], "%f", &st2);
          if(st2 < st1){
             printf("****ERROR****, 2nd intensity value smaller than 1st, existing\n\n");
             exit(1);
          }
       }  
    }
    
    printf("%d %d\n", upos1, upos2);
    
    if((upos1 < 0 || upos1 > 1) || (upos2 < 0 || upos2 > 1)) {
       printf("****ERROR****, incorrect specifier for whether to use location of additional fields.\n\n");
       exit(1);
    }
    
    if(str_type < 0 || str_type > 2){
       printf("****ERROR****, wrong intensity ID\n\n");
       exit(1);
    }
    
    if(sep_type < 0 || sep_type > 1){
       printf("****ERROR****, wrong seperation ID\n\n");
       exit(1);
    }
    
    if(sep_ws < 0 || sep_ws > 1){
       printf("****ERROR****, wrong track length ID\n\n");
       exit(1);
    }
    
    if(tsamp1.itsmp < 0 || tsamp1.itsmp > 1){
       printf("****ERROR****, wrong track sampling ID for ensemble 1\n\n");
       exit(1);
    } 
    
    if(tsamp2.itsmp < 0 || tsamp2.itsmp > 1){
       printf("****ERROR****, wrong track sampling ID for ensemble 2\n\n");
       exit(1);
    } 
    
    if(sep_ws){
       printf("what is the number of points to match over from the start of track overlap?\n\n");
       scanf("%d", &numst);
       if(numst <= 0){
          printf("****ERROR****, can't have negative or zero number of points.\n\n");
	  exit(1);
       }
    }


    printf("**********************************************************\r\n"
           "   WARNINIG WARNINIG WARNINIG WARNINIG WARNINIG WARNINIG   \n");
    if(ftype == 's'){
       printf("****INFORMATION****, the chosen file type is '%c' which requires a scaler file\r\n"
              "                     to be used, i.e. tr_trs of ff_trs files.                   \n", ftype);
    }
    else if(ftype == 'v'){
       printf("****INFORMATION****, the chosen file type is '%c' which requires a vector file\r\n"
              "                     to be used, i.e. phase_trs files.                          \n", ftype);
       sscanf(&ctype[1], "%d", &vtype);
    }
    else{
       printf("****ERROR****, file type '%c' not known.\n\n", ftype);
       exit(1);
    }
    printf("**********************************************************\r\n");

    printf("%c %d\n", ftype, vtype);

    fcmp = fopen("cmp.dat", "r");
    if(fcmp){
       printf("****WARNING****, file exists with subsidary bining data for \r\n"
              "                 intensity and seperation distributions     \r\n"
              "                 this will be used instead of defaults.     \n\n");
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f", &tpbw);
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f %f %d", &id1, &id2, &nid);
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f %f %d", &idd1, &idd2, &nidd);
       fclose(fcmp); fcmp = NULL;

    }

    fcmp = fopen("region.dat", "r");
    if(fcmp){
       printf("****WARNING****, file exists with region data for sub-setting tracks,   \r\n"
              "                 this information may be used to restrict the analysis  \r\n"
              "                 to the specified region.                               \n\n");
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f %f", &(reg.rlng1), &(reg.rlng2));
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f %f", &(reg.rlat1), &(reg.rlat2));
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%ld %ld", &(reg.rtim1), &(reg.rtim2));
       if(fgets(buff, MAXCHR, fcmp)) {
          sscanf(buff, "%f %f", &stf1, &stf2);
          if(stf1 > stf2){
             printf("****ERROR****, intensity thresholds specified in region.dat file are inconsistent, stf1 > stf2.\n\n");
             exit(1);
          }
          if(fgets(buff, MAXCHR, fcmp)) {
             sscanf(buff, "%d %d", &loc1, &loc2);
          }
          else {
             printf("****ERROR****, need additional field location information when applying intensity thresholds. \n\n");
             exit(1);
          }
          istf = 1;
          printf("****WARNING****, using intensity thresholding on chosen system extreme type.\n\n");
       }

       ireg = 1;
       if(reg.rlng2 < reg.rlng1){
          printf("****WARNING****, stradeling the Greenwich meridian\n\n");
          reg.ig = 1;
       }
       fclose(fcmp); fcmp = NULL;

       if(reg.reg_type < -1 || reg.reg_type > 3){
          printf("****ERROR****, region type %d not valid.\n\n", reg.reg_type);
          exit(1);
       }

       if(reg.reg_type == 1 && !str_type){
          printf("****ERROR****, illegal combination of region type and intensity type.\n\n");
          exit(1);
       }

    }
    else reg.reg_type = -1;

    if(reg.reg_type < 0){
       ireg = 0;
       printf("****INFORMATION****, not using region sub-setting.\n\n");
    }
    else printf("****INFORMATION****, using region type %d\n\n", reg.reg_type);

    fcmp = fopen("orog.dat", "r");
    if(fcmp){
       printf("****WARNING****, file exists with orographic data for masking,          \r\n"
              "                 this information will be used to restrict the analysis \r\n"
              "                 to regions below a specified orographic level. This    \r\n"
              "                 assumes data is surface geopotential in units of metres.\n\n");
       printf("Specify an level in metres.\n\n");
       scanf("%f", &forog);

       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%d %d", &ilngo, &ilato);
       olng = (float *)calloc(ilngo, sizeof(float));
       olat = (float *)calloc(ilato, sizeof(float));

       orogm = (float *)calloc(ilngo*ilato, sizeof(float));
       for(i=0; i < ilngo; i++) fscanf(fcmp, "%f", olng + i);
       for(i=0; i < ilato; i++) fscanf(fcmp, "%f", olat + i);
       fgets(buff, MAXCHR, fcmp);
       fgets(buff, MAXCHR, fcmp);
       fread(orogm, ilngo*ilato*sizeof(float), 1, fcmp);

       fclose(fcmp); fcmp = NULL;

    }


    fin1 = fopen(argv[1], "r");
    if(!fin1){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[1]);
       exit(1);
    }

    tr1 = read_tracks(fin1, &trnum1, &gpr1, &ipr1, ftype, &alat1, &alng1, NULL, NULL);
    nff1 = nff;
    nfld1 = nfld;
    nfwpos1 = nfwpos;
    nfwpos = NULL;
    fseek(fin1, 0L, SEEK_SET);
    tm1 = read_tracks(fin1, &trnum1, &gpr1, &ipr1, ftype, &alat1, &alng1, NULL, NULL);

    nfl1= nff;

    iper_num_t = iper_num;
    sum_per_t = sum_per;

    fclose(fin1);
    
    if(tsamp1.itsmp){
       if(tr1->time){
          printf("What is the starting time (YYYYMMDDHH) and current timestep (hrs) for first track data set?\n\n");
          scanf("%ld %d", &(tsamp1.stim), &(tsamp1.tstp));
	  printf("What is the time sampling, e.g. 1, 2, 3...?\n\n");
	  scanf("%d", &(tsamp1.tsmp));
	  
	  trsmp1 = time_sample(tr1, trnum1, &tsamp1);
	  tmsmp1 = time_sample(tm1, trnum1, &tsamp1);
	  
	  printf("Output matches at original sampling or new sampling, '0' or '1'\n\n");
	  scanf("%d", &(tsamp1.iosmp));
	  
	  if(tsamp1.iosmp < 0 || tsamp1.iosmp > 1){
	     printf("****ERROR****, output sampling Id. not known.\n\n");
	     exit(1);
	  }

/* swap pointers */
	  
	  trtmp = tr1;
	  tr1 = trsmp1;
	  trsmp1 = trtmp;
	  
	  trtmp = tm1;
	  tm1 = tmsmp1;
	  tmsmp1 = trtmp;
       }
       else{
          printf("Time sampling only possible using real time, convert to real time.\n\n");
	  exit(1);
       }
    }

    fin2 = fopen(argv[2], "r");
    if(!fin2){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[2]);
       exit(1);
    }

    tr2 = read_tracks(fin2, &trnum2, &gpr2, &ipr2, ftype, &alat2, &alng2, NULL, NULL);
    nff2 = nff;
    nfld2 = nfld;
    nfwpos2 = nfwpos;
    nfwpos = NULL;
    fseek(fin2, 0L, SEEK_SET);
    tm2 = read_tracks(fin2, &trnum2, &gpr2, &ipr2, ftype, &alat2, &alng2, NULL, NULL);

    nfl2 = nff;

    fclose(fin2);
    
    if(tsamp2.itsmp){
       if(tr2->time){
          printf("What is the starting time (YYYYMMDDHH) and current timestep (hrs) for first track data set?\n\n");
          scanf("%ld %d\n", &(tsamp2.stim), &(tsamp2.tstp)); 
	  printf("What is the new timestep (hrs)?\n\n");
	  scanf("%d", &(tsamp2.tsmp)); 

	  trsmp2 = time_sample(tr2, trnum2, &tsamp2);
	  tmsmp2 = time_sample(tm2, trnum2, &tsamp2);
	  
	  printf("Output matches at original sampling or new sampling, '0' or '1'\n\n");
	  scanf("%d", &(tsamp1.iosmp));
	  
	  if(tsamp1.iosmp < 0 || tsamp1.iosmp > 1){
	     printf("****ERROR****, output sampling Id. not known.\n\n");
	     exit(1);
	  }
	  
	  trtmp = tr2;
	  tr2 = trsmp2;
	  trsmp2 = trtmp;
	  
	  trtmp = tm2;
	  tm2 = tmsmp1;
	  tmsmp2 = trtmp;
	  
       }  
       else {
          printf("Time sampling only possible using real time, convert to real time.\n\n");
	  exit(1);       
       }
    }
    
    
 /* allocate track arrays for re-ordering tracks */
 
    if(irord){
    
       if(tsamp1.itsmp || tsamp2.itsmp){
          printf("****ERROR****, can't currently use track re-ordering with track time sampling.\n\n");
	  exit(1);
       }
 
       ntm1 = (struct tot_tr * )calloc(trnum1, sizeof(struct tot_tr));
       mem_er((ntm1 == NULL) ? 0 : 1, trnum1 * sizeof(struct tot_tr));
    
       ntm2 = (struct tot_tr * )calloc(trnum2, sizeof(struct tot_tr));
       mem_er((ntm2 == NULL) ? 0 : 1, trnum2 * sizeof(struct tot_tr));  
    
    } 

    if(nfl1 != nfl2){
       printf("****WARNING****, there are different numbers of additional fields in the track files.\n\n");

    }

    if(infl1 < 0 || infl1 > nfl1) {
       printf("****ERROR***, requested default or additional field does not exist for first track file.\n\n");
       exit(1);
    }

    if(infl2 < 0 || infl2 > nfl2) {
       printf("****ERROR***, requested default or additional field does not exist for second track file.\n\n");
       exit(1);
    }

    if(loc1 < 0 || loc1 > nfl1){
       printf("****ERROR***, LOC1, requested default or additional field does not exist for first track file.\n\n");
       exit(1);
    }

    if(loc2 < 0 || loc2 > nfl2){
       printf("****ERROR***, LOC2, requested default or additional field does not exist for first track file.\n\n");
       exit(1);
    }

    if(istf){

       if(loc1){
          if(! *(nfwpos1 + loc1 - 1)){
             printf("****ERROR****, LOC1, additional fields do not have location information for output to str.dat.\n\n");
             exit(1);
          }
          lloc1=0;
          for(i=0; i < loc1; i++){
             if(*(nfwpos1 + i)) lloc1 += 3;
             else lloc1 += 1;
          }
          --lloc1;
       }

       if(loc2){
          if(! *(nfwpos2 + loc2 - 1)){
             printf("****ERROR****, LOC2, additional fields do not have location information for output to str.dat.\n\n");
             exit(1);
          }
          lloc2=0;
          for(i=0; i < loc2; i++){
             if(*(nfwpos2 + i)) lloc2 += 3;
             else lloc2 += 1;
          }
          --lloc2;
       }


    }

    if(infl1) {
      iff1 = 0;
      for(i=0; i < infl1; i++){
         if(*(nfwpos1 + i)) iff1 += 3;
         else iff1 += 1;
      }
      --iff1;
      if(*(nfwpos1 + infl1 - 1)) ifpos1 = iff1 - 2;
      else upos1 = 0;

    }
    else upos1 = 0;

    if(infl2) {
      iff2 = 0;
      for(i=0; i < infl2; i++){
         if(*(nfwpos2 + i)) iff2 += 3;
         else iff2 += 1;
      }
      --iff2;
      if(*(nfwpos2 + infl2 - 1)) ifpos2 = iff2 - 2;
      else upos2 = 0;

    }
    else upos2 = 0;

    temp = (float *)calloc(nbin, sizeof(float));
    mem_er((temp == NULL) ? 0 : 1, nbin * sizeof(float));


    nmbin = (int)floor((100.0 * (1.0 - tolnum) / tpbw) + 0.5) + 1;

    for(i=0; i < 4; i++){
        hist[i] = (float *)calloc(nbin + 1, sizeof(float));
        mem_er((hist[i] == NULL) ? 0 : 1, (nbin + 1) * sizeof(float));
    }

    temp2 = (float *)calloc(nmbin, sizeof(float));
    mem_er((temp2 == NULL) ? 0 : 1, nmbin * sizeof(float));

    temp3 = (float *)calloc(nid, sizeof(float));
    mem_er((temp3 == NULL) ? 0 : 1, nid * sizeof(float));

    temp4 = (float *)calloc(nidd, sizeof(float));
    mem_er((temp4 == NULL) ? 0 : 1, nidd * sizeof(float));

    nptm = (float *)calloc(nmbin, sizeof(float));
    mem_er((nptm == NULL) ? 0 : 1, nmbin * sizeof(float));

    npid = (float *)calloc(nid, sizeof(float));
    mem_er((npid == NULL) ? 0 : 1, nid * sizeof(float));

    npdd = (float *)calloc(nidd, sizeof(float));
    mem_er((npdd == NULL) ? 0 : 1, nidd * sizeof(float));

    ntpdt = (float **)calloc(nidd, sizeof(float *));
    mem_er((ntpdt == NULL) ? 0 : 1, nidd * sizeof(float *));

    for(i=0; i < nidd; i++){
        *(ntpdt + i) = (float *)calloc(nmbin, sizeof(float));
        mem_er((*(ntpdt + i) == NULL) ? 0 : 1, nmbin * sizeof(float));        
    }
    
/* assign temperature bins */

    dt = (t2 - t1) / nbin;
    did = (id2 - id1) / nid;
    ddd = (idd2 - idd1) / nidd;

    for(i=0; i < nbin; i++) *(temp + i) = t1 + 0.5 * dt + i * dt; 


    nm1 = (int *)calloc(trnum1, sizeof(int));
    mem_er((nm1 == NULL) ? 0 : 1, trnum1 * sizeof(int));

/* assign bins for the percentage of points that match */

    tpbw /= 100.0;

    for(i=0; i < nmbin; i++) *(temp2 + i) = tolnum + tpbw * ((float)i + 0.5);

    for(i=0; i < nid; i++) *(temp3 + i) = id1 + 0.5 * did + i * did;
    
    for(i=0; i < nidd; i++) *(temp4 + i) = idd1 + 0.5 * ddd + i * ddd;

    if((gpr1 != gpr2 || ipr1 != ipr2) ||
        (fabs(alat1 - alat2) > 1.0e-4 || fabs(alng1 - alng2) > 1.0e-4)){
        printf("***ERROR***, projection parameters do not match\n\n");
        exit(1);
    }

    if((tr1->time && !tr2->time) || (!tr1->time && tr2->time)){
       printf("***ERROR***, temporal identities do not match\n\n");
       exit(1);
    }

    if(iper_num != iper_num_t || fabs(sum_per - sum_per_t) > 1.0e-6){
       printf("****ERROR****, number of time periods do not match\n");
       printf("               iper_num= %d and %d, sum_per= %f and %f\n\n", iper_num, iper_num_t, sum_per, sum_per_t);
       exit(1);

    }


    fstr = fopen("str.dat", "w");
    if(!fstr){
       printf("****ERROR****, cannot open file for write.\n\n");
       exit(1); 
    }

    if(IDIFF){
       fdiff = fopen("diff.dat", "w");
       if(!fdiff){
          printf("****ERROR****, cannot open file for write.\n\n");
          exit(1); 
       }
    }

    for(i=0; i < trnum1; i++){

        atr1 = tr1 + i;

        if(atr1->num > 0){

/* determine max or mean intensity */

           str1 = intensity(atr1, &reg, str_type, infl1, iff1, ftype, vtype, &imax1, loc1, lloc1, ioff);

           if(str1 > MISS) continue;

/*           if(str_type){
              if(MATCH_STR){
                 if(ftype == 's') strm = (infl2) ? (atr1->trpt + imax1)->add_fld[iff2] : (atr1->trpt + imax1)->zf;
                 else strm = phpstr(vtype, atr1->trpt + imax1);
              }
           }
*/

/* check if track goes through region */

           if(ireg) inreg1 = inregion(atr1, &reg , imax1, str_type, loc1, lloc1);
           else inreg1 = 1;

           if(!inreg1) continue;

/* check if intensity extrema are within thresholds */

           if(istf){
              if((stf2 - str1) * (str1 - stf1) >= 0.0) instr1 = 1;
              else instr1 = 0;
           }
           else instr1 = 1;

           if(!instr1) continue;

/* check if extrema should be masked by orography*/

           if(orogm && str_type && inreg1){
              fp1 = atr1->trpt + imax1;
              orog_nearest(fp1, olng, olat, ilngo, ilato, &ilng, &ilat);
              if(*(orogm + ilat * ilngo + ilng) > forog) inreg1 = 0;
           }

           if(inreg1 && (str1 < t1 || str1 > t2)){
             printf("***ERROR***, strength outside range %f, for track %d.\n\n", str1, i+1);
             exit(1);
           }

           distmin = DISTLARGE;
           distm = DISTLARGE;
           jtrmin = -1;
           numpm = 0;
           frat = 0.0;

           nmat = 0;

           for(j=0 ; j < trnum2; j++){

              atr2 = tr2 + j;

              if(atr2->num > 0){

                 str = intensity(atr2, &reg, str_type, infl2, iff2, ftype, vtype, &imax, loc2, lloc2, ioff);

                 if(str > MISS) continue;
/* check if track goes through region */

                 if(ireg) inreg = inregion(atr2, &reg, imax, str_type, loc2, lloc2);
                 else inreg = 1;

                 if(!inreg) continue;

/* check if intensity extrema are within thresholds */

                 if(istf){
                    if((stf2 - str) * (str - stf1) >= 0.0) instr = 1;
                    else instr = 0;
                 }
                 else instr = 1;

                 if(!instr) continue;

                 if(orogm && str_type && inreg){
                    fp2 = atr2->trpt + imax;
                    orog_nearest(fp2, olng, olat, ilngo, ilato, &ilng, &ilat);
                    if(*(orogm + ilat * ilngo + ilng) > forog) inreg = 0;
                 }

/* Check for overlap in time */	 

                 if(toverlap(atr1, atr2, &is1, &is2)){

                    dist = trdist(atr1, atr2, tolmatch, is1, is2, &nump, &it1, &it2, sep_type, sep_ws, numst, upos1, upos2, ifpos1, ifpos2);

                    ffrat = 2.0 * nump / (float)(atr1->num + atr2->num); 
		    
                    if(ffrat >= tolnum && dist <= tolmatch){

                       ++nmat;

                       if(jtrmin < 0 ){

                          jtrmin = j;
                          distmin = (sep_type == 0) ? dist : trdist(atr1, atr2, tolmatch, is1, is2, &nump, &it1, &it2, 0, sep_ws, numst, upos1, upos2, ifpos1, ifpos2);
                          distm = dist;
                          frat = ffrat;
                          numpm = nump;
                          in1 = it1;
                          in2 = it2;
                          inreg2 = inreg;
                          instr2 = instr;		  
                          imax2 = imax;
                          str2 = str;

                       }

                       else if(jtrmin >= 0 && dist < distm) {

                          jtrmin = j;
                          distmin = (sep_type == 0) ? dist : trdist(atr1, atr2, tolmatch, is1, is2, &nump, &it1, &it2, 0, sep_ws, numst, upos1, upos2, ifpos1, ifpos2);
                          distm = dist;
                          frat = ffrat;
                          numpm = nump;
                          in1 = it1;
                          in2 = it2;
                          inreg2 = inreg;
                          instr2 = instr;
                          imax2 = imax;
                          str2 = str;

                       }


                    }


                 }


              }

           }


           if(jtrmin >= 0 && inreg1 && inreg2 && instr1 && instr2 && imax1 >= 0 && imax2 >= 0){ 

              if(str2 > MISS) continue;

              *(nm1 + i) = nmat;
              nmatch += 1.0;

              if(str2 < t1 || str2 > t2){
                printf("***ERROR***, strength outside range %f.\n\n", str2);
                exit(1);
              }

              ++(*(hist[0] + (int)((str1 - t1) / dt)));
              ++(*(hist[2] + (int)((str2 - t1) / dt)));

              if(!MATCH_STR){
                 iprint = 0;
                 if(ifpos1 >= 0 && ifpos2 >= 0){
                    if((atr1->trpt + imax1)->add_fld[ifpos1] < MISS && ((tr2 + jtrmin)->trpt + imax2)->add_fld[ifpos2] < MISS) iprint = 1;
                    else {
                       if(TRACK_MAX_OUT) {
		         (tm1 + i)->num = 0; (tm2 + jtrmin)->num = 0;
			 if(tsamp1.itsmp)(tmsmp1 + i)->num = 0; 
			 if(tsamp2.itsmp)(tmsmp2 + jtrmin)->num = 0;
		       }
                    }
                 }
                 else if(ifpos1 >= 0){
                    if((atr1->trpt + imax1)->add_fld[ifpos1] < MISS) iprint = 1;
                    else{
                       if(TRACK_MAX_OUT) {
		         (tm1 + i)->num = 0;
			 if(tsamp1.itsmp) (tmsmp1 + i)->num = 0;
		       }
                    }   
                 }
                 else if(ifpos2 >= 0){
                    if(((tr2 + jtrmin)->trpt + imax2)->add_fld[ifpos2] < MISS) iprint = 1;
                    else {
                       if(TRACK_MAX_OUT) {
		         (tm2 + jtrmin)->num = 0;
			 if(tsamp2.itsmp) (tmsmp2 + jtrmin)->num = 0;
		       }
                    }
                 }
                 else iprint = 1;

                 if(iprint) {
                    nint1 = nint2 = 0;
                    for(j=0; j < atr1->num; j++){
                       fptmp = atr1->trpt + j;
                       if(ftype == 's') str = (infl1) ? fptmp->add_fld[iff1] : fptmp->zf;
                       else str = phpstr(vtype, fptmp);
                       if((st2 - str) * (str - st1) >= 0.0) ++nint1;

                    }
                    for(j=0; j < (tr2 + jtrmin)->num; j++){
                       fptmp = (tr2 + jtrmin)->trpt + j;
                       if(ftype == 's') str = (infl2) ? fptmp->add_fld[iff2] : fptmp->zf;
                       else str = phpstr(vtype, fptmp);
                       if((st2 - str) * (str - st1) >= 0.0) ++nint2;
                    }

                    if(loc1){
                       xpos1 = (atr1->trpt + imax1)->add_fld[lloc1-2];
                       ypos1 = (atr1->trpt + imax1)->add_fld[lloc1-1];
                       if(xpos1 > MISS){
                          xpos1 = (atr1->trpt + imax1)->xf;
                          ypos1 = (atr1->trpt + imax1)->yf;
                       }
                    }
                    else {
                       xpos1 = (atr1->trpt + imax1)->xf;
                       ypos1 = (atr1->trpt + imax1)->yf;
                    }

                    if(loc2){
                       xpos2 = ((tr2 + jtrmin)->trpt + imax2)->add_fld[lloc2-2];
                       ypos2 = ((tr2 + jtrmin)->trpt + imax2)->add_fld[lloc2-1];
                       if(xpos2 > MISS){
                          xpos2 = ((tr2 + jtrmin)->trpt + imax2)->xf;
                          ypos2 = ((tr2 + jtrmin)->trpt + imax2)->yf;
                       }                      
                    }
                    else {
                       xpos2 = ((tr2 + jtrmin)->trpt + imax2)->xf;
                       ypos2 = ((tr2 + jtrmin)->trpt + imax2)->yf;
                    }

                    if(atr1->time)
                       fprintf(fstr, "%d %d %d %d %10ld %10ld %f %f %f %f %f %f %d %d %d %d\n", i, jtrmin, imax1, imax2, (atr1->trpt + imax1)->time, ((tr2 + jtrmin)->trpt + imax2)->time, xpos1, ypos1, str1, xpos2, ypos2, str2, atr1->num, (tr2 + jtrmin)->num, nint1, nint2);
                    else
                       fprintf(fstr, "%d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d\n", i, jtrmin, imax1, imax2, (atr1->trpt + imax1)->fr_id, ((tr2 + jtrmin)->trpt + imax2)->fr_id, xpos1, ypos1, str1, xpos2, ypos2, str2, atr1->num, (tr2 + jtrmin)->num, nint1, nint2);
                    }
              }
              else {
                 iprint = 0;
		 if(tr1->time  && tr2->time){
		    for(j=0; j < (tr2 + jtrmin)->num; j++) {
		        if((atr1->trpt + imax1)->time == ((tr2 + jtrmin)->trpt + j)->time) imax2 = j;
		    }
		 }
		 else {
		    for(j=0; j < (tr2 + jtrmin)->num; j++) {
		        if((atr1->trpt + imax1)->fr_id == ((tr2 + jtrmin)->trpt + j)->fr_id) imax2 = j;
		    }
		 }
		 
                 if(ftype == 's') strm = (infl2) ? ((tr2 + jtrmin)->trpt + imax2)->add_fld[iff2] : ((tr2 + jtrmin)->trpt + imax2)->zf;
                 else strm = phpstr(vtype, (tr2 + jtrmin)->trpt + imax2);		 
		 
                 if(ifpos1 >= 0 && ifpos2 >= 0){
                    if((atr1->trpt + imax1)->add_fld[ifpos1] < MISS && ((tr2 + jtrmin)->trpt + imax2)->add_fld[ifpos2] < MISS) iprint = 1;
                    else {
                       if(TRACK_MAX_OUT) {
		         (tm1 + i)->num = 0; (tm2 + jtrmin)->num = 0;
			 if(tsamp1.itsmp)(tmsmp1 + i)->num = 0; 
			 if(tsamp2.itsmp)(tmsmp2 + jtrmin)->num = 0;
		       }
                    }
                 }
                 else if(ifpos1 >= 0){
                    if((atr1->trpt + imax1)->add_fld[ifpos1] < MISS) iprint = 1;
                    else{
                       if(TRACK_MAX_OUT) {
		         (tm1 + i)->num = 0;
			 if(tsamp1.itsmp) (tmsmp1 + i)->num = 0;
		       }
                    }
                 }
                 else if(ifpos2 >= 0){
                    if(((tr2 + jtrmin)->trpt + imax2)->add_fld[ifpos2] < MISS) iprint = 1;
                    else {
                       if(TRACK_MAX_OUT) {
		         (tm2 + jtrmin)->num = 0;
			 if(tsamp2.itsmp) (tmsmp2 + jtrmin)->num = 0;
		       }
                    }
                 }
                 else iprint = 1;

                 if(iprint){
                    nint1 = nint2 = 0;
                    for(j=0; j < atr1->num; j++){
                       fptmp = atr1->trpt + j;
                       if(ftype == 's') str = (infl1) ? fptmp->add_fld[iff1] : fptmp->zf;
                       else str = phpstr(vtype, fptmp);
                       if((st2 - str) * (str - st1) >= 0.0) ++nint1;
                    }
                    for(j=0; j < (tr2 + jtrmin)->num; j++){
                       fptmp = (tr2 + jtrmin)->trpt + j;
                       if(ftype == 's') str = (infl2) ? fptmp->add_fld[iff2] : fptmp->zf;
                       else str = phpstr(vtype, fptmp);
                       if((st2 - str) * (str - st1) >= 0.0) ++nint2;
                    }

                    if(loc1){
                       xpos1 = (atr1->trpt + imax1)->add_fld[lloc1-2];
                       ypos1 = (atr1->trpt + imax1)->add_fld[lloc1-1];
                       if(xpos1 > MISS){
                          xpos1 = (atr1->trpt + imax1)->xf;
                          ypos1 = (atr1->trpt + imax1)->yf;
                       }
                    }
                    else {
                       xpos1 = (atr1->trpt + imax1)->xf;
                       ypos1 = (atr1->trpt + imax1)->yf;
                    }

                    if(loc2){
                       xpos2 = ((tr2 + jtrmin)->trpt + imax2)->add_fld[lloc2-2];
                       ypos2 = ((tr2 + jtrmin)->trpt + imax2)->add_fld[lloc2-1];
                       if(xpos2 > MISS){
                          xpos2 = ((tr2 + jtrmin)->trpt + imax2)->xf;
                          ypos2 = ((tr2 + jtrmin)->trpt + imax2)->yf;
                       }
                    }
                    else {
                       xpos2 = ((tr2 + jtrmin)->trpt + imax2)->xf;
                       ypos2 = ((tr2 + jtrmin)->trpt + imax2)->yf;
                    }

                    if(atr1->time)
                       fprintf(fstr, "%d %d %d %d %10ld %10ld %f %f %f %f %f %f %d %d %d %d\n", i, jtrmin, imax1, imax2, (atr1->trpt + imax1)->time, ((tr2 + jtrmin)->trpt + imax2)->time, xpos1, ypos1, str1, xpos2, ypos2, strm, atr1->num, (tr2 + jtrmin)->num, nint1, nint2);
                    else
                       fprintf(fstr, "%d %d %d %d %d %d %f %f %f %f %f %f %d %d %d %d\n", i, jtrmin, imax1, imax2, (atr1->trpt + imax1)->fr_id, ((tr2 + jtrmin)->trpt + imax2)->fr_id, xpos1, ypos1, str1, xpos2, ypos2, strm, atr1->num, (tr2 + jtrmin)->num, nint1, nint2);

                 }
              }

              ++(*(nptm + (int)((frat - tolnum) / tpbw)));
              atr1->num = 0;
              (tr2 + jtrmin)->num = 0;
	      if(tsamp1.itsmp) (trsmp1 + i)->num = 0;
	      if(tsamp2.itsmp) (trsmp2 + jtrmin)->num = 0;
	      
	      if(irord){
	         *(ntm1 + nmtch) = *(tm1 + i);
		 *(ntm2 + nmtch) = *(tm2 + jtrmin);
	         ++nmtch;
	      }

              totnum += numpm;

              for(j=0; j < numpm; j++){

                  fp1 = atr1->trpt + in1 + j;
                  fp2 = (tr2 + jtrmin)->trpt + in2 + j;
/*                  ditty = fp2->zf - fp1->zf; */
                  ss1 = (infl1) ? fp1->add_fld[iff1] : fp1->zf;
                  ss2 = (infl2) ? fp2->add_fld[iff2] : fp2->zf;
                  if(ss1 > MISS || ss2 > MISS) continue;
                  ditty = ss2 - ss1;

                  if(IDIFF) fprintf(fdiff, "%e %f %f %f %f %ld\n", ditty, fp1->xf, fp1->yf, fp2->xf, fp2->yf,
		  fp1->time); 
		  
                  if(ditty < id1 || ditty > id2){
                     if(!ior){
                        printf("****WARNING****, track point intensity difference out of range for histogram.\n\n");

                        ior = 1;
                     }
                     printf("Value = %f\n", ditty);
                  }

                  else  ++(*(npid + (int)((ditty - id1) / did)));

              }

              if(distmin > idd2){

                 if(!iorr){
                    printf("****WARNING****, track point seperation distance out of range for histogram.\n\n");

                    iorr = 1;
                 }
                 printf("Value = %f\n", distmin);
              }

              else  {
                   ++(*(npdd + (int)((distmin - idd1) / ddd)));
                   ++(*(*(ntpdt + (int)((distmin - idd1) / ddd)) + (int)((frat - tolnum) / tpbw)));
              }


              
           }

           else if(inreg1){
              ++(*(hist[1] + (int)((str1 - t1) / dt)));

              nomch1 += 1.0;
           }

           if(!inreg1) atr1->num = -1;
           if(jtrmin >=0 && !inreg2) (tr2 + jtrmin)->num = -1;

/*           if(!inreg1) atr1->num = 0;
           if(!inreg2) (tr2 + jtrmin)->num = 0; */


        }

    }

    if(IDIFF) fclose(fdiff);
    
    fclose(fstr);

    for(i=0; i < trnum2; i++){

        atr2 = tr2 + i;

        if(atr2->num > 0){

            str2 = intensity(atr2, &reg, str_type, infl2, iff2, ftype, vtype, &imax, loc2, lloc2, ioff);
            if(str2 > MISS) continue; 

            if(ireg) inreg = inregion(atr2, &reg, imax, str_type, loc2, lloc2);
            else inreg = 1;

            if(inreg){
	    
               if(str2 < t1 || str2 > t2){
                  printf("***ERROR***, mean strength outside range %f.\n\n", str2);
                  exit(1);
               }

               ++(*(hist[3] + (int)((str2 - t1) / dt)));
               nomch2 += 1.0;
            }
            else {
               atr2->num = -1;
            }

        }

    }
    
    printf("%d %d\n", (int) nmatch, nmtch);

    for(i=0; i < trnum1; i++){
        atr1 = tr1 + i;
        atr2 = tm1 + i;

        if(atr1->num > 0) {
	   atr2->num = 0;
	   if(tsamp1.itsmp) (tmsmp1 + i)->num = 0;
	}
        else if(atr1->num < 0) {
	   atr1->num = 0; atr2->num = 0;
	   if(tsamp1.itsmp) {(trsmp1 + i)->num = 0; (tmsmp1 + i)->num = 0;}
	}

    }

    for(i=0; i < trnum2; i++){
        atr1 = tr2 + i;
        atr2 = tm2 + i;

        if(atr1->num > 0) {
	   atr2->num = 0;
	   if(tsamp2.itsmp) (tmsmp2 + i)->num = 0;
	}
        else if(atr1->num < 0) {
	   atr1->num = 0; atr2->num = 0;
	   if(tsamp2.itsmp) {(trsmp2 + i)->num = 0; (tmsmp2 + i)->num = 0;}	   
	}


    }


/* write out sub-sets of track ensembles */

    nff = nff1;
    nfld = nfld1;
    nfwpos = nfwpos1;
    
    if(tsamp1.itsmp && !tsamp1.iosmp){
       trtmp = tr1;
       tr1 = trsmp1;
       trsmp1 = trtmp;
	  
       trtmp = tm1;
       tm1 = tmsmp1;
       tmsmp1 = trtmp;    
    }

    strcpy(filout, FILENAME);
    strcat(filout, "ens1_yes.dat");

    fout = fopen(filout, "w");

    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    if(irord){
        meantrd(fout, ntm1, nmtch, ftype, gpr1, ipr1, alat1, alng1);   
    }
    else {
       meantrd(fout, tm1, trnum1, ftype, gpr1, ipr1, alat1, alng1);
    }

    fclose(fout); fout = NULL;

    strcpy(filout, FILENAME);
    strcat(filout, "ens1_no.dat");
    
    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    meantrd(fout, tr1, trnum1, ftype, gpr1, ipr1, alat1, alng1);

    fclose(fout); fout = NULL;

    nff = nff2;
    nfld = nfld2;
    nfwpos = nfwpos2;
    
    if(tsamp2.itsmp && !tsamp2.iosmp){
       trtmp = tr2;
       tr2 = trsmp2;
       trsmp2 = trtmp;
	  
       trtmp = tm2;
       tm2 = tmsmp2;
       tmsmp2 = trtmp;    
    }

    strcpy(filout, FILENAME);
    strcat(filout, "ens2_yes.dat");
    
    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    if(irord){
       meantrd(fout, ntm2, nmtch, ftype, gpr2, ipr2, alat2, alng2);    
    }
    else {
       meantrd(fout, tm2, trnum2, ftype, gpr2, ipr2, alat2, alng2);
    }

    fclose(fout); fout = NULL;

    strcpy(filout, FILENAME);
    strcat(filout, "ens2_no.dat");
    
    fout = fopen(filout, "w");
    if(!fout){
       printf("***ERROR***, unable to open file %s for 'w'\n\n", filout);
       exit(1);
    }

    meantrd(fout, tr2, trnum2, ftype, gpr2, ipr2, alat2, alng2);

    fclose(fout); fout = NULL;



    for(i=0; i < trnum1; i++){
        if(*(nm1 + i) > 1) ntrm += 1.0;
    }


    if(iper_num && !idisty){

       ntrm /= sum_per;


       for(i=0; i < nbin; i++){
           for(j=0; j < 4; j++)
               *(hist[j] + i) /= sum_per;
       }

/*       for(i=0; i < nmbin; i++)
           *(nptm + i) /= sum_per; */

/*       for(i=0; i < nid; i++)
           *(npid + i) /= nmatch; */

   

    }

    if(idisty == 1){
       for(i=0; i < nbin; i++) {
           *(hist[0] + i) /= nmatch;
           *(hist[1] + i) /= nomch1;
           *(hist[2] + i) /= nmatch;
           *(hist[3] + i) /= nomch2;
       }

    }

    else if(idisty == 2){
       for(i=0; i < nbin; i++) {
           *(hist[0] + i) /= (nmatch * dt);
           *(hist[1] + i) /= (nomch1 * dt);
           *(hist[2] + i) /= (nmatch * dt);
           *(hist[3] + i) /= (nomch2 * dt);
       }
    }

    for(i=0; i < nbin; i++){
        for(j=0; j < 4; j++)
            *(hist[j] + nbin) += *(hist[j] + i);
    }


    trpm1 = (float)trnum1 / sum_per;
    trpm2 = (float)trnum2 / sum_per;

    printf("# Tracks that match is %f\n\n", nmatch);

    for(i=0; i < nbin; i++)printf("%9.6f\n", *(hist[0] + i));

    printf("************ SUMMARY STATISTICS ********************\n\n");

    printf("First file %s, second file %s\n\n", argv[1], argv[2]);
    printf("Distance threshold is %f deg.\n", tolmatch);
    printf("Relative No. of points threshold is %f\n", tolnum);
    printf("No. of periods is %f \n", sum_per);
    printf("No. of tracks in ensemble 1 is %d\n", trnum1);
    printf("No. of tracks per month in ensemble 1 is %f\n", trpm1); 
    printf("No. of tracks in ensemble 2 is %d\n", trnum2);
    printf("No. of tracks per month in ensemble 1 is %f\n", trpm2);
    printf("No of tracks that match for distance and No. of points thresholds is %f\n", nmatch / sum_per);
    printf("No. of tracks in ensemble 1 that match more than one track \r\n"
           "in ensemble 2 for distance and No. of points thresholds is %f\n\n", ntrm);

    printf("Percentage of tracks that match in ensemble 1 is %f\n", *(hist[0] + nbin) / trpm1);
    printf("Percentage of tracks that don't match in ensemble 1 is %f\n", *(hist[1] + nbin) / trpm1);

    printf("Percentage of tracks that match in ensemble 2 is %f\n", *(hist[2] + nbin) / trpm2);
    printf("Percentage of tracks that don't match in ensemble 2 is %f\n", *(hist[3] + nbin) / trpm2);

    printf("*****************Histograms**************************\n\n");
    printf("          1st Match   1st No Match   2nd Match   2nd No Match\n\n");
    for(i=0; i < nbin; i++){
        printf("%9.6f    %9.6f          %9.6f        %9.6f           %9.6f\n", *(temp + i),  *(hist[0] + i), *(hist[1] + i), *(hist[2] + i), *(hist[3] + i));
    }
    printf("Sum          %9.6f          %9.6f        %9.6f           %9.6f\n",   *(hist[0] + nbin), *(hist[1] + nbin), *(hist[2] + nbin), *(hist[3] + nbin));

    printf("-----------------------------------------------------\n");
    printf("      Number of Points That Match Distribution       \n");
    printf("Fract.  ");
    for(i=0; i<nmbin; i++)printf("%6.3f ", *(temp2 + i));
    printf("\n");
    printf("No.     ");
    for(i=0; i<nmbin; i++)printf("%6.2f ", *(nptm + i) / (nmatch * tpbw));
    printf("\n");

    printf("\n\n");

    printf("-----------------------------------------------------\n");
    printf("      Distribution for intensity differences         \n");
    printf("Intens.  ");
    for(i=0; i<nid; i++)printf("%8.4f ", *(temp3 + i));
    printf("\n");
    printf("No.      ");
    for(i=0; i<nid; i++)printf("%8.5f ", *(npid + i) / (totnum * did));

    printf("\n\n");

    printf("\n\n");

    printf("-----------------------------------------------------\n");
    printf("      Distribution for seperation distances          \n");
    printf("Seperat..  ");
    for(i=0; i<nidd; i++)printf("%8.4f ", *(temp4 + i));
    printf("\n");
    printf("No.      ");
    for(i=0; i<nidd; i++)printf("%8.5f ", *(npdd + i) / (nmatch * ddd));

    printf("\n\n");

/* write 2D statistics */

    ftpdt = fopen("temp-dist.stat", "w");

/*     fprintf(ftpdt, "%d %d\n", nidd, nmbin); */
    for(i=0; i < nidd; i++){

/*        for(j=0; j < nmbin; j++) fprintf(ftpdt, "%8.4f ", *(*(ntpdt + i) + j) / (nmatch * nmbin * ddd * tpbw));
        fprintf(ftpdt, "\n"); */
          for(j=0; j < nmbin; j++) 
             fprintf(ftpdt, "%8.4f %8.4f %8.4f \n", *(temp2 + j), *(temp4 + i), *(*(ntpdt + i) + j) / (nmatch * nmbin * ddd * tpbw));
        

    }

    fclose(ftpdt);
    
/* clean up */

    for(i=0; i < trnum1; i++) {
       atr1 = tr1 + i;
       atr2 = tm1 + i;
       if(nff1){
         for(j=0; j < atr1->num; j++) free((atr1->trpt + j)->add_fld);
	 for(j=0; j < atr2->num; j++) free((atr2->trpt + j)->add_fld);
       }
       free(atr1->trpt);
       free(atr2->trpt);
    }
    
    free(tr1);
    free(tm1);
    if(irord) free(ntm1);
    
    for(i=0; i < trnum2; i++) {
       atr1 = tr2 + i;
       atr2 = tm2 + i;
       if(nff2){
         for(j=0; j < atr1->num; j++) free((atr1->trpt + j)->add_fld);
	 for(j=0; j < atr2->num; j++) free((atr2->trpt + j)->add_fld);
       }
       free(atr1->trpt);
       free(atr2->trpt);
    }
    
    free(tr2);
    free(tm2);
    if(irord) free(ntm2);    

    return 0;

}

float phpstr(int vtype, struct fet_pt_tr *fp)
{

   float str;
   switch(vtype){
   case 0:
      str = fp->zf; 
      break;
   case 1:
      str = fp->tend;
      break;
   case 2:
      str = fp->gwthr;
      break;
   case 3:
      str = fp->vec[0];
      break;
   case 4:
      str = fp->vec[1];
      break;
   default:
      printf("****ERROR****, str type not known for this file, exiting.\n\n");
      exit(1);
   }
   return str;
}


void orog_nearest(struct fet_pt_tr *fp, float *olng, float *olat, int ilngo, int ilato, int *ilng, int *ilat)
{

     int il, ir, ix, iy;

     float fa, fb;

     il = 0;
     ir = ilngo - 1;

     if((*(olng + ir) - fp->xf) < 0.0) *ilng = ir;
     else if((fp->xf - *olng) < 0.0) *ilng = il;
     else {
        while(ir - il > 1){
           ix = (il + ir) / 2;
           if((*(olng + ir) - fp->xf) * (fp->xf - *(olng + ix)) > 0.0) il = ix; 
           else ir = ix;
        }

        fa = fabs(*(olng + il) - fp->xf);
        fb = fabs(*(olng + ir) - fp->xf);
        *ilng = (fa < fb) ? il : ir;
     }   

     il = 0;
     ir = ilato - 1;

     if((*(olat + ir) - fp->yf) < 0.0) *ilat = ir;
     else if((fp->yf - *olat) < 0.0) *ilat = il;
     else {
        while(ir - il > 1){
           iy = (il + ir) / 2;
           if((*(olat + ir) - fp->yf) * (fp->yf - *(olat + iy)) > 0.0) il = iy;
           else ir = iy;
        }

        fa = fabs(*(olat + il) - fp->yf);
        fb = fabs(*(olat + ir) - fp->yf);
        *ilat = (fa < fb) ? il : ir;

     }

     return;

}

int inregion(struct tot_tr *atr, REG *reg, int imax, int str_type, int loc, int lloc)
{

    int i;
    int ntim=0;
    int nreg=0;
    int inreg=0;

    long int tim;

    struct fet_pt_tr *fp;

    for(i=0; i < atr->num; i++){
        fp = atr->trpt + i;
        if(atr->time) tim = fp->time;
        else tim = fp->fr_id;
        if(tim >= reg->rtim1 && tim <= reg->rtim2) ++ntim; 
    }
    if((float) ntim / (float) (atr->num) >= TREG) {
        if(! reg->reg_type){
           for(i=0; i < atr->num; i++){
              fp = atr->trpt + i;
              if(region(fp, reg, atr->time, loc, lloc)){inreg = 1; break;}
           }
         }
/* check is maximum is in region */

         else if(str_type && (reg->reg_type == 1 || reg->reg_type == 3)){

            fp = atr->trpt + imax;
            if(region(fp, reg, atr->time, loc, lloc))inreg = 1;
         }
/* check if n% of points are in region */
         else if(reg->reg_type == 2){
           for(i=0; i < atr->num; i++){
              fp = atr->trpt + i;
              nreg += region(fp, reg, atr->time, loc, lloc);
            }
            if((float)nreg / (float)(atr->num) >= PERC) inreg = 1;
            if(str_type){
               fp = atr->trpt + imax;
               if(!region(fp, reg, atr->time, loc, lloc))inreg = 0;
            }

         }

    }

    return inreg;

}

float intensity(struct tot_tr *atr, REG * reg, int str_type, int infl, int iff, int ftype, int vtype, int *imax, int loc, int lloc, int ioff)
{

    int i;
    int npt=0;

    float str=0.0, stmp=0.0;

    struct fet_pt_tr *fp;

    *imax = 0;

    if(!str_type){

       for(i=0; i < atr->num; i++) {
          fp = atr->trpt + i;
          if(ftype == 's') stmp = (infl) ? fp->add_fld[iff] : fp->zf;
          else stmp = phpstr(vtype, fp);
          if(stmp < MISS){
             str += stmp;
             ++npt;
          }

       }

       if(npt) str /= npt;
       else str = ADD_UNDEF;

    }

    else if(str_type && reg->reg_type == 3){
       *imax = 0;
       str = ADD_UNDEF;
       for(i=0; i < atr->num; i++) {
           fp = atr->trpt + i;
           if(region(fp, reg, atr->time, loc, lloc)){
              if(ftype == 's') str = (infl) ? fp->add_fld[iff]: fp->zf;
              else str = phpstr(vtype, fp);
              *imax = i;
           }
           if(str < MISS) break;
       }
       if(str > MISS) return str;
       for(i=0; i < atr->num; i++) {
           fp = atr->trpt + i;
           if(region(fp, reg, atr->time, loc, lloc)){
              if(ftype == 's') stmp = (infl) ? fp->add_fld[iff] : fp->zf;
              else stmp = phpstr(vtype, fp);
              if(stmp < MISS){
                 if(str_type == 1){
                    if(stmp > str) {str = stmp; *imax = i;}
                 }
                 else if(str_type == 2){
                    if(stmp < str) {str = stmp; *imax = i;}
                 }

              }
           }
       }
       
       if(ioff){
          *imax += ioff;
	  if(*imax < 0 || *imax > atr->num - 1) *imax = -1;
       }

    }

    else{
       *imax = 0;
       for(i=0; i < atr->num; i++) {
           fp = atr->trpt + i;
           if(ftype == 's') str = (infl) ? fp->add_fld[iff]: fp->zf;
           else str = phpstr(vtype, fp);
           *imax = i;
           if(str < MISS) break;
       }
       if(str > MISS) return str;
       for(i=0; i < atr->num; i++) {
           fp = atr->trpt + i;
           if(ftype == 's') stmp = (infl) ? fp->add_fld[iff] : fp->zf;
           else stmp = phpstr(vtype, fp);
           if(stmp < MISS){
              if(str_type == 1){
                 if(stmp > str) {str = stmp; *imax = i;}
              }
              else if(str_type == 2){
                 if(stmp < str) {str = stmp; *imax = i;}
              }
           }
       }
       
       if(ioff){
          *imax += ioff;
	  if(*imax < 0 || *imax > atr->num - 1) *imax = -1;
       }

    }

    return str;
}

int region(struct fet_pt_tr *fp, REG *reg, long int ttim, int loc, int lloc)
{

    int inreg=0;
    long int tim;
    
    float xpos=0.0, ypos=0.0;

    if(reg->reg_type < 0) {inreg = 1; return inreg;}

    if(ttim) tim = fp->time;
    else tim = fp->fr_id;
    
    if(!loc){
       xpos = fp->xf;
       ypos = fp->yf;
    }
    else{
       xpos = fp->add_fld[lloc-2];
       ypos = fp->add_fld[lloc-1];
    }
    
    if(xpos > MISS) return inreg;
    
    if(! reg->ig){
      if((reg->rlng2 - xpos) * (xpos - reg->rlng1) >= 0.0 &&
         (reg->rlat2 - ypos) * (ypos - reg->rlat1) >= 0.0 &&
         (tim >= reg->rtim1 && tim <= reg->rtim2)                   ){
         inreg = 1;
      }
    }
    else {
      if((reg->rlat2 - ypos) * (ypos - reg->rlat1) >= 0.0 &&
         (tim >= reg->rtim1 && tim <= reg->rtim2)               && 
         ((reg->rlng2 - xpos) * (xpos - 0.0) >= 0.0 || 
         (360.0 - xpos) * (xpos - reg->rlng1) >= 0.0)){
         inreg = 1;
      }

    }

    return inreg;
}

struct tot_tr *time_sample(struct tot_tr *trr, int trnum, TSAMP *tsamp)
{

    int i=0, j=0;
    int ifr=0;
    int ipos=0, idf=0;
    int noff = 0;
    int itrid=0;
    int ntmp=0, npt=0;
    int num=0;
    int id1=0;
    
    long int ta;
    
    struct tot_tr *trsmp=NULL;
    struct tot_tr *atr=NULL, *antr=NULL;
    struct fet_pt_tr *at=NULL, *anr=NULL;
    
    trsmp = (struct tot_tr * )calloc(trnum, sizeof(struct tot_tr));
    mem_er((trsmp == NULL) ? 0 : 1, trnum * sizeof(struct tot_tr));
    
    for(i=0; i < trnum; i++){
        atr = trr + i;
        antr = trsmp + i;
	memcpy(antr, atr, sizeof(struct tot_tr));
        ntmp = (atr->num / tsamp->tsmp) + tsamp->tsmp;
        antr->trpt = (struct fet_pt_tr * )calloc(ntmp, sizeof(struct fet_pt_tr));
        mem_er((antr->trpt == NULL) ? 0 : 1, ntmp * sizeof(struct fet_pt_tr));
    }
    
    for(i=0; i < trnum; i++){
    
        atr = trr + i;
	antr = trsmp + i;
	
	if(atr->time){
	   ta = tsamp->stim;
	   ifr = 1;
	
	   if(ta > atr->time){
	      printf("Starting time for sampling occurs later than start of track.\n\n");
	      exit(1);
	   }
	
	   while(ta < atr->time){
               ta = new_time(ta, tsamp->tstp);
               ++ifr;
           }
	   
	   num = 0;
           for(j=0; j < atr->num; j++) {
              at = atr->trpt + j;
              at->fr_id = ifr + num;
              num += 1 + at->nfm;
	      
	      if(at->nfm){
	         printf("****ERROR****, missing time steps not currently supported.\\n\n");
		 exit(1);
	      }
           }
       
       }


       npt = 0;
       ipos = atr->num + 1;

       for(j=0; j < atr->num; j++) {
           at = atr->trpt + j;
           idf = at->fr_id - 1 - noff;
           if(!(idf % tsamp->tsmp) && idf >= 0) {ipos = j; itrid = (idf / tsamp->tsmp) + 1; break;}
       }
       
       if(ipos >= atr->num) continue;
       
       memcpy(antr->trpt + npt, at, sizeof(struct fet_pt_tr));
       id1 = npt; ++npt; antr->num = npt;
       anr = antr->trpt + id1;
       anr->fr_id = itrid;
       ++itrid;
       
       ipos += tsamp->tsmp;
       at += tsamp->tsmp;
       
       while(ipos < atr->num){
           memcpy(antr->trpt + npt, at, sizeof(struct fet_pt_tr)); 
           anr = antr->trpt + npt;
           ++npt; antr->num = npt;
           anr->fr_id = itrid;
	   ++itrid;
	   ipos += tsamp->tsmp;
           if(ipos >= atr->num) break;
           at += tsamp->tsmp;
       }
       
       antr->num = npt;

    }

    return trsmp;
}
