#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"
#include "vec.h"
#include "geo_values.h"
#include "st_fo.h"

#define  NFSTT   64
#define  TSTEP   6
#define  ILARGE  100000000

#define  TIMTOL  6
#define  TMSTOL  1
#define  DISTOL  5.0

/* program to compute landfall and properties */

typedef struct lndf {
    int tstep;
    int trid;
    int lfdir;
    long int time;
    float lng;
    float lat;
    float str;
    VEC vp;
} LNDF;

typedef struct lsmask {
    int ilm;
    int lglng;
    int lglat;
    int *ilms;
    float *glng;
    float *glat;
} LSM;


struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
LNDF *find_lndf(struct tot_tr * , VEC * , int , int * , int , int , int , int , LSM * );
void free_tracks(struct tot_tr * , int );
long int new_time(long int , int );
long int timesep(long int , long int );
void convert_track_add(struct tot_tr * , int , int , int );
int intersect_t(VEC * , VEC * , VEC * );
double geodist(VEC * , VEC *);

int orog_test(float * , float * , int , int , int * , float , float );

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int tf=4;

int main(int argc, char **argv)
{
    int i=0, j=0;
    int trnum=0;
    int gpr=0, ipr=0;    
    int ncnt=0;
    int ilndf=0;
    int ilng=0, ilat=0, iint=0;
    int upos=0, ipos=0;
    int lsmd=0;
    
    int nlnd=0;
    
    char buff[MAXCHR];
    char filin[MAXCHR];
    char filcnt[MAXCHR];

    FILE *fin=NULL;
    FILE *fcnt=NULL;
    FILE *fout=NULL;
    FILE *flm=NULL;
    
    
    float alat=0.0, alng=0.0;

    float *lngcnt=NULL, *latcnt=NULL;
    
    float *lms=NULL;
    
    double phi=0.0, thet=0.0;
    double s1=0.0, s2=0.0, c1=0.0, c2=0.0;

    struct tot_tr *tracks=NULL;
    
    VEC *vcnt=NULL;
    
    LNDF *lndf=NULL;
    
    LSM lsm={0, 0, 0, NULL, NULL, NULL};
    
    if(argc != 6) {
      printf("Usage: landfall [trackfile] [boundary file] [longitude id] [latitude id] [intensity id]\n\n");
      exit(1);
    }
    
    sscanf(argv[1], "%s", filin);
    sscanf(argv[2], "%s", filcnt);
    
    printf("Do you want just first landfall or all land-ocean transitions, '0' for all and '1' for first landfall only.\n\n");
    scanf("%d", &ilndf);
    
    fcnt = fopen(filcnt, "r");
    if(!fcnt){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filcnt);
       exit(1);
    } 
    
    sscanf(argv[3], "%d", &ilng);    
    sscanf(argv[4], "%d", &ilat);
    sscanf(argv[5], "%d", &iint);
    
    if(!ilng || !ilat){
       if(ilng != ilat){
          printf("****ERROR****, if default location is required both longitude and latitude must be the defaults.\n\n");
	  exit(1);
       }
       
       upos = 0;
    }
    else {upos = 1; ipos = ilng - 1;}

    fgets(buff, MAXCHR, fcnt);
    sscanf(buff, "%d", &ncnt);
    
    lngcnt = (float *)calloc(ncnt, sizeof(float));
    mem_er((lngcnt == NULL) ? 0 : 1, ncnt*sizeof(float));
    
    latcnt = (float *)calloc(ncnt, sizeof(float));
    mem_er((latcnt == NULL) ? 0 : 1, ncnt*sizeof(float));    
      
    vcnt = (VEC *)calloc(ncnt, sizeof(VEC));
    mem_er((vcnt == NULL) ? 0 : 1, ncnt*sizeof(VEC));
    
    for(i=0; i < ncnt; i++){
        fgets(buff, MAXCHR, fcnt);
	sscanf(buff, "%f %f", lngcnt+i, latcnt+i);
	
	phi = *(lngcnt + i) * FP_PI; 
        thet = FP_PI2 - *(latcnt + i) * FP_PI;
        if(thet < 0.) thet = 0.;
	          
	sincos(phi, &s1, &c1);
        sincos(thet, &s2, &c2);

        (vcnt + i)->x = s2 * c1;
        (vcnt + i)->y = s2 * s1;
        (vcnt + i)->z = c2;
	
    }
        
    fclose(fcnt); 
    
    flm = fopen("lmask.dat", "r");
    if(flm){
       printf("****INFORMATION****, using land mask file lmask.dat\n\n");
       lsm.ilm = 1;
       fgets(buff, MAXCHR, flm);
       sscanf(buff, "%d %d", &(lsm.lglng), &(lsm.lglat));
       lsmd = lsm.lglng * lsm.lglat;
       lsm.glng = (float *)calloc(lsm.lglng, sizeof(float));
       lsm.glat = (float *)calloc(lsm.lglat, sizeof(float));
       lms = (float *)calloc(lsmd, sizeof(float));
       lsm.ilms = (int *)calloc(lsmd, sizeof(int));
       for(i=0; i < lsm.lglng; i++) fscanf(flm, "%f", lsm.glng + i); 
       for(i=0; i < lsm.lglat; i++) fscanf(flm, "%f", lsm.glat + i);
       fgets(buff, MAXCHR, flm);
       fgets(buff, MAXCHR, flm);
       fread(lms, lsm.lglng * lsm.lglat*sizeof(float), 1, flm);
       fclose(flm);
       for(i=0; i < lsmd; i++) *(lsm.ilms + i) = (int)(*(lms + i)); 

    }

    fout = fopen("landfall.dat", "w");
	  
		
    printf("%s\n", filin);

    fin = fopen(filin, "r");
    if(!fin){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
       exit(1);
    }

    tracks = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
        
    fclose(fin); fin = NULL;
       
    convert_track_add(tracks, trnum, upos, ipos);
		
/* landfall */	

    for(i=0; i < trnum; i++){	

        lndf = find_lndf(tracks+i, vcnt, ncnt, &nlnd, ilndf, ilng, ilat, iint, &lsm);
		

	if(tracks->trpt->time){
           for(j=0; j < nlnd; j++) fprintf(fout, "%d %d %ld %f %f %e %d\n", i+1, (lndf + j)->trid, (lndf + j)->time, (lndf + j)->lng, (lndf + j)->lat, (lndf + j)->str, (lndf + j)->lfdir);
        }
        else {
	   for(j=0; j < nlnd; j++) fprintf(fout, "%d %d %d %f %f %e %d\n", i+1, (lndf + j)->trid, (lndf + j)->tstep, (lndf + j)->lng, (lndf + j)->lat, (lndf + j)->str, (lndf + j)->lfdir);
	}
	
				
        free(lndf);
     
     }
		
     free_tracks(tracks, trnum);
        
     fclose(fout);
    
     free(vcnt);

     return 0;
}


LNDF *find_lndf(struct tot_tr *track, VEC *vcnt, int ncnt, int *nlnd, int ilndf, int ilng, int ilat, int iint, LSM *lsm)
{

    int i=0, j=0, k=0;
    int ic=0, iout=0;
    int iflnd=0;
    int ivf=0;
    int inr=0;
    int nlp1=0, nlp2=0;
    int ioff=0;
    
    float xx=0.0, yy=0.0;
    
    double norm=0.0;
    
    struct fet_pt_tr *fp=NULL, *fp1=NULL;
    
    VEC *v1=NULL, *v2=NULL;
    VEC a, b, ab, vf;
    VEC vv;
    
    LNDF *lndf=NULL;

    ic = 0;
    
    *nlnd = 0;
	
    for(i=1; i < track->num; i++){
	fp = track->trpt + i - 1; 
	fp1 = track->trpt + i;
	
        if(fp->pp[0] > ADD_CHECK || fp1->pp[0] > ADD_CHECK) continue;
	if(iint && fp->add_fld[iint - 1] > ADD_CHECK) continue;
	    
        a.x = fp->pp[0];
        a.y = fp->pp[1];
	a.z = fp->pp[2];
	    
	b.x = fp1->pp[0];
	b.y = fp1->pp[1];
	b.z = fp1->pp[2];
	    
	crosp(&a, &b, &ab);
	
	norm = sqrt(dotp(&ab, &ab));
		
	normv(&ab, norm);	
	    
        for(j=1; j < ncnt; j++){
	
	    ivf = 0;
	    
	    v1 = vcnt + j - 1;
            v2 = vcnt + j;
	    
	    crosp(v1, v2, &vv);
	    
	    norm = sqrt(dotp(&vv, &vv));
		
	    normv(&vv, norm);	
	    
	    crosp(&vv, &ab, &vf); 
	    
	    norm = sqrt(dotp(&vf, &vf));
		
	    normv(&vf, norm);
	    
	    if(intersect_t(&vf, v1, v2) && intersect_t(&vf, &a, &b)) ivf = 1;
	    
	    if(!ivf){
	    
	      crosp(&ab, &vv, &vf);  
		
	      norm = sqrt(dotp(&vf, &vf));
		
	      normv(&vf, norm);
	      
	      if(intersect_t(&vf, v1, v2) && intersect_t(&vf, &a, &b)) ivf = 1;
	      
	    }
		

	       
            if(ivf){	      
	
	       iout = 0;	   
	       if(!ic) {
	          ic = 1; 
		  iout = 1;
		  lndf = (LNDF *)calloc(1, sizeof(LNDF));
		  mem_er((lndf == NULL) ? 0 : 1, sizeof(LNDF));
		  *nlnd = 1;
		  lndf->lfdir = 0;
	       }
	       else {
	       
	          inr = 0;
		  for(k=0; k < *nlnd; k++){
		      if(track->time){
		         if((labs(timesep((lndf + *nlnd - 1)->time, fp->time)) <= TIMTOL || 
			     labs(timesep((lndf + *nlnd - 1)->time, fp1->time)) <= TIMTOL   ) ||
			     geodist(&((lndf + *nlnd - 1)->vp), &vf) <= DISTOL                 ) inr = 1;
		      }
		      else {
		         if((abs((lndf + *nlnd - 1)->tstep - fp->fr_id) <= TMSTOL || 
			     abs((lndf + *nlnd - 1)->tstep - fp1->fr_id) <= TMSTOL   ) ||
			     geodist(&((lndf + *nlnd - 1)->vp), &vf) <= DISTOL          ) inr = 1;
		      }
		      
		  } 

                  if(!inr){
		    iout = 1;
		    *nlnd += 1;
		    lndf = (LNDF *)realloc_n(lndf, *nlnd * sizeof(LNDF));
		    mem_er((lndf == NULL) ? 0 : 1, *nlnd * sizeof(LNDF));
		    (lndf + *nlnd - 1)->lfdir = 0;
		  }

	       }	  
		
		  
	       if(iout){
	          if(ilndf) iflnd = 1;
		  
		  if(lsm->ilm){
		     nlp1 = orog_test(lsm->glng, lsm->glat, lsm->lglng, lsm->lglat, lsm->ilms, fp->xf, fp->yf);
		     nlp2 = orog_test(lsm->glng, lsm->glat, lsm->lglng, lsm->lglat, lsm->ilms, fp1->xf, fp1->yf);

		     ioff = 1;
		     while(nlp1 == nlp2){
		         if(i - ioff - 1 >= 0)
			    nlp1 = orog_test(lsm->glng, lsm->glat, lsm->lglng, lsm->lglat, lsm->ilms, (fp - ioff)->xf, (fp - ioff)->yf);
			 if(i + ioff < track->num)
			    nlp2 = orog_test(lsm->glng, lsm->glat, lsm->lglng, lsm->lglat, lsm->ilms, (fp1 + ioff)->xf, (fp1 + ioff)->yf);
			  
		         ++ioff;
		         if(ioff > 3) break;
		     }
		     
		     if(nlp1 < nlp2) {(lndf + *nlnd - 1)->lfdir = 1;}
		     else if(nlp1 > nlp2) {(lndf + *nlnd - 1)->lfdir = -1;}
		     else {(lndf + *nlnd - 1)->lfdir = 0;}
		    
		  }
		      
	          yy = (FP_PI2 - acos(vf.z)) / FP_PI;
	          if(yy > 90.0) yy = 90.0;
                  else if(yy < -90.0) yy = -90.0;
	          xx = atan2(vf.y, vf.x) / FP_PI;
	          if(xx < 0.0) xx = 360.0 + xx;
		  
	          (lndf + *nlnd - 1)->lng = xx;
	          (lndf + *nlnd - 1)->lat = yy;
	          ((lndf + *nlnd - 1)->vp).x = vf.x;
	          ((lndf + *nlnd - 1)->vp).y = vf.y;
	          ((lndf + *nlnd - 1)->vp).z = vf.z;
	          (lndf + *nlnd - 1)->time = fp->time;
	          if(!iint)(lndf + *nlnd - 1)->str = fp->zf;
	          else (lndf + *nlnd - 1)->str = fp->add_fld[iint - 1];
	          (lndf + *nlnd - 1)->tstep = fp->fr_id;
	          (lndf + *nlnd - 1)->trid = track->trid; 
	       
	       }
		
            }			
		
        }   
	

        if(iflnd) break;
	
    }


    return lndf;

}

void free_tracks(struct tot_tr *alltr, int tr_count)
{
    int i, j;

    struct tot_tr *trr=NULL;
    struct fet_pt_tr *atr=NULL;

    for(i=0; i < tr_count; i++) {
       trr = alltr + i;
       for(j=0; j < trr->num; j++) {
           atr = (trr->trpt) + j;
           if(nfld) free(atr->add_fld);
       }
       free(trr->trpt);
    }

    free(alltr);
    free(nfwpos);

    return;

}

long int timesep(long int t1, long int t2)
{
    int nhr=0;
    
    if(t1 == t2) return 0;
    
    else if(t1 < t2){
    
       while(t1 < t2){
           --nhr;  
           t1 = new_time(t1, TSTEP);
       }
    
    }
    else {
       while(t2 < t1){
           ++nhr;  
           t2 = new_time(t2, TSTEP);
       }
    
    }

    return nhr * TSTEP;

}

int intersect_t(VEC *vv, VEC *v1, VEC *v2)
{
    int isect=0;
    
    double dp1=0.0, dp2=0.0;
    double dvv=0.0;   
   
    dvv = dotp(v1, v2);
    dp1 = dotp(v1, vv);
    dp2 = dotp(v2, vv);
	    
		
    if(dvv <= ((dp1 < dp2) ? dp1: dp2)) isect = 1; 
    
    return isect;

}

double geodist(VEC *lnd, VEC *vf)
{
   double ang=0.0;
  
   ang = acos(lnd->x * vf->x + lnd->y * vf->y + lnd->z * vf->z) / FP_PI;

   return ang;
}
