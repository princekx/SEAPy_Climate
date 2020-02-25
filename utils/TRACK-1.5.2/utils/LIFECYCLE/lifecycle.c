#include <Stdio.h>
#include <stdlib.h>
#include "splice.h"

#define  PLARGE   1.0e+10
#define  NLARGE  -1.0e+10

#define  NBIN     201

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

int tom='g';

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void centre(struct tot_tr * , REG * , int , int , int );
int region(struct fet_pt_tr * , REG * , long int );

int main(void)
{

    int i, j;
    int trnum=0;
    int gpr, ipr;
    int itr=0, ift=0, ifr=0, ifft=0, iffr=0;
    int iref=0;
    int itype=0;
    int icent=NBIN/2;
    int nbin[NBIN];
    int ib=0;
    int ireg=0;
    
    char buff[MAXCHR];

    float alat, alng;
    float str;
    float astr[NBIN];
    float sum=0.0;
    float sum10=0.0;

    FILE *ftr=NULL, *flif=NULL, *ftnd=NULL;
    FILE *fcmp=NULL;

    char filin[200];

    struct tot_tr *tracks=NULL, *trr=NULL;
    struct fet_pt_tr *fp=NULL;
    
    REG reg={0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0};
    
    fcmp = fopen("region.dat", "r");
    if(fcmp){
       printf("****WARNING****, file exists with region data for sub-setting tracks, this    \r\n"
              "                 information may be used to restrict the analysis to the      \r\n"
              "                 specified region.                                            \r\n"
	      "                                                                              \r\n"
	      "                 Note region matching currently only wrt to tracked locations.\n\n");
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f %f", &(reg.rlng1), &(reg.rlng2));
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%f %f", &(reg.rlat1), &(reg.rlat2));
       fgets(buff, MAXCHR, fcmp);
       sscanf(buff, "%ld %ld", &(reg.rtim1), &(reg.rtim2));
       ireg = 1;
       if(reg.rlng2 < reg.rlng1){
          printf("****WARNING****, stradeling the Greenwich meridian\n\n");
          reg.ig = 1;
       }
       fclose(fcmp); fcmp = NULL; 
       
       printf("What type of region matching is required for system centering, \r\n"
              "Input  '0' extrema in region                                   \r\n"
	      "       '1' extrema of region                                   \n\n");
       scanf("%d", &(reg.reg_type));
       if(reg.reg_type < 0 || reg.reg_type > 1){
          printf("****ERROR****, incorrect region matching specified.\n\n");
	  exit(1);  
       }  
    }    
    else reg.reg_type = -1;    

    printf("What is the track file to read?\n\n");
    scanf("%s", filin);

    ftr = fopen(filin, "r");
    if(!ftr){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
       exit(1);
    }

    tracks = read_tracks(ftr, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(ftr);

    printf("Do you want an individual track or the average over all tracks, '0' for single or '1' for average.\n\n");
    scanf("%d", &itype);

    if(!itype){
       printf("Which track is required?\n\n");
       scanf("%d", &itr);
       if(itr < 1 || itr > trnum){
          printf("****ERROR****, track number does not exist.\n\n");
	  exit(1);
       }
    }
    else{
       for(i=0; i < NBIN; i++) {astr[i] = 0.0; nbin[i] = 0;}
    }

    printf("Center tracks wrt '0' max. intensity     \r\n"
           "                  '1' min. intensity     \n\n");
    scanf("%d", &iref);

    if(nff){
       printf("****INFORMATION****, there are %d additional fields associated with these tracks \r\n", nff);

       printf("Which field is required for history, input '0' for default?\n\n");
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

       printf("Which field is required for reference, e.g. max intensity?\n\n");
       scanf("%d", &ifr);
       iffr = 0;
       if(ifr){
          for(i=0; i < ifr; i++){
              if(*(nfwpos + i)) iffr += 3;
              else iffr += 1;
          }
          --iffr;
       }


    }

    flif = fopen("history.dat", "w");
    if(!flif){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", "history.dat");
       exit(1);
    }

    if(!itype){

       trr = tracks + itr - 1;

       centre(trr, &reg, iref, ifr, iffr);

       for(i=0; i < trr->num; i++){
          fp = trr->trpt + i;
          str = (ift) ? fp->add_fld[ifft] : fp->zf;
          fprintf(flif, "%d %e\n", fp->fr_id, str);
       }

    }
    else {
       for(i=0; i < trnum; i++){
           trr = tracks + i;
           centre(trr, &reg, iref, ifr, iffr);

           for(j=0; j < trr->num; j++){
               fp = trr->trpt + j;

               str = (ift) ? fp->add_fld[ifft] : fp->zf;
	       if(str > ADD_CHECK) continue;
               
	       ib = icent + fp->fr_id;
               if(ib < 0 || ib >= NBIN) continue;
               astr[ib] += str;
               nbin[ib] += 1;

           }
       }

       for(i=0; i < NBIN; i++) {
           if(nbin[i]) astr[i] /= (float)nbin[i];
           sum += astr[i];
           fprintf(flif, "%d %e\n", -icent + i, astr[i]);
       }

       for(i=90; i< 111; i++) sum10 += astr[i];

       ftnd = fopen("tendency.dat", "w");
       if(!ftnd){
          printf("***ERROR***, unable to open file %s for 'r'\n\n", "tendency.dat");
          exit(1);
       }

       for(i=0; i < NBIN-1; i++)
           fprintf(ftnd, "%d %e\n", -icent + i, astr[i+1] - astr[i]);

       fclose(ftnd);
       
    }

    fclose(flif);

    printf("sum=%f sum10=%f\n", sum, sum10);

    return 0;
}


void centre( struct tot_tr *trr, REG *reg, int iref, int ifr, int iffr)
{
    int i;
    int ipref=-1;
    int frid=0;

    float ref, str;

    struct fet_pt_tr *fp=NULL;

    if(!iref) ref = NLARGE;
    else ref = PLARGE;

    for(i=0; i < trr->num; i++){
       fp = trr->trpt + i;
       fp->fr_id = i+1;
    }

/* find reference point */

/* extrema in region or no region */
    
    if(reg->reg_type <= 0){
       for(i=0; i < trr->num; i++){
          fp = trr->trpt + i;
          str = (ifr) ? fp->add_fld[iffr] : fp->zf;
          if(str > ADD_CHECK) continue;
          if(!iref){
             if(str > ref) {ref = str; ipref = i;}
          }
          else{
             if(str < ref) {ref = str; ipref = i;}
          }
       }

       if(!region(trr->trpt + ipref, reg, trr->time)) ipref = -1;
       
    }
    
/* extrema of region */

    else {

       for(i=0; i < trr->num; i++){
          fp = trr->trpt + i;
          if(!region(fp, reg, trr->time)) continue;
          str = (ifr) ? fp->add_fld[iffr] : fp->zf;
          if(str > ADD_CHECK) continue;
          if(!iref){
             if(str > ref) {ref = str; ipref = i;}
          }
          else{
             if(str < ref) {ref = str; ipref = i;}
          }
       }    
    
    }
    
    if(ipref >= 0){frid = (trr->trpt + ipref)->fr_id;}
    else {
      trr->num = 0;
      return;
    }

    for(i=0; i < trr->num; i++){
       fp = trr->trpt + i;
       fp->fr_id -= frid;
    }


    return;
}

int region(struct fet_pt_tr *fp, REG *reg, long int ttim)
{

    int inreg=0;
    long int tim;
    
    float xpos=0.0, ypos=0.0;

    if(reg->reg_type < 0) {inreg = 1; return inreg;}

    if(ttim) tim = fp->time;
    else tim = fp->fr_id;
    
    xpos = fp->xf;
    ypos = fp->yf;
    
    if(! reg->ig){
      if((reg->rlng2 - xpos) * (xpos - reg->rlng1) >= 0.0 &&
         (reg->rlat2 - ypos) * (ypos - reg->rlat1) >= 0.0 &&
         (tim >= reg->rtim1 && tim <= reg->rtim2)            ){
         inreg = 1;
      }
    }
    else {
      if((reg->rlat2 - ypos) * (ypos - reg->rlat1) >= 0.0 &&
         (tim >= reg->rtim1 && tim <= reg->rtim2)         && 
         ((reg->rlng2 - xpos) * (xpos - 0.0) >= 0.0 || 
         (360.0 - xpos) * (xpos - reg->rlng1) >= 0.0)        ){
         inreg = 1;
      }

    }

    return inreg;
}
