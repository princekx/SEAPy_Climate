#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"
#include "vec.h"
#include "geo_values.h"
#include "st_fo.h"

#define  NFSTT   64
#define  TSTEP   6
#define  ILARGE  100000000

/* program to compute landfall and properties */

typedef struct lndf {
    int tstep;
    int trid;
    long int time;
    float lng;
    float lat;
    float str;
    VEC vp;
} LNDF;


struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
LNDF *find_lndf(struct tot_tr * , VEC * , int , int * , int , int , int , int );
void free_tracks(struct tot_tr * , int );
long int new_time(long int , int );
long int timesep(long int , long int );
void lead_time_error(LNDF * , LNDF * , double * , double * , double * , double * , double * , int * , int , int , long int );
void write_stats(double * , double * , double * , double * , double * , int * , int );
void convert_track_add(struct tot_tr * , int , int , int );

extern float sum_per;
extern int iper_num;
extern int nfld, nff;
extern int *nfwpos;

int tf=4;

int main(int argc, char **argv)
{
    int i, j;
    int trnum=0, trnummn=0;
    int gpr=0, ipr=0;    
    int ncnt=0;
    int ntst=0;
    int ilndf=0;
    int idet=0;
    int ilng=0, ilat=0, iint=0;
    int upos=0, ipos=0;
    
    int nlnd=0, nlnd_an=0, nlnd_m=0;
    
    int nmem[NFSTT]={0}, nmean[NFSTT]={0}, ncntrl[NFSTT]={0}, ndet[NFSTT]={0}, nspr[NFSTT]={0};
    
    long int fcst=0;
    
    char buff[MAXCHR];
    char filin[MAXCHR];
    char filcnt[MAXCHR];
    char dirnam[MAXCHR];

    FILE *fin=NULL;
    FILE *fcnt=NULL;
    FILE *fout=NULL;
    
    float alat=0.0, alng=0.0;

    float *lngcnt=NULL, *latcnt=NULL;
    
    double phi=0.0, thet=0.0;
    double s1=0.0, s2=0.0, c1=0.0, c2=0.0;
    
    double tmem[NFSTT]={0}, tmean[NFSTT]={0}, tcntrl[NFSTT]={0}, tspr[NFSTT]={0};
    double tabmem[NFSTT]={0}, tabmean[NFSTT]={0}, tabcntrl[NFSTT]={0}, tabspr[NFSTT]={0};
    double dmem[NFSTT]={0}, dmean[NFSTT]={0}, dcntrl[NFSTT]={0}, dspr[NFSTT]={0};
    double smem[NFSTT]={0}, smean[NFSTT]={0}, scntrl[NFSTT]={0}, sspr[NFSTT]={0};
    double sabmem[NFSTT]={0}, sabmean[NFSTT]={0}, sabcntrl[NFSTT]={0}, sabspr[NFSTT]={0};
    
    struct tot_tr *tracks=NULL, *trkmn=NULL, *trr=NULL;
    
    VEC *vcnt=NULL;
    
    LNDF *lndf_an=NULL, *lndf=NULL, *lndf_m=NULL;
    
    struct dirent *entry=NULL, *entryt=NULL;
    
    DIR *dir=NULL, *dirt=NULL;
    
    if(argc != 7) {
      printf("Usage: landfall [path to dirs] [longitude id] [latitude id] [intensity id] [Match dir] [0 or 1 for Det.]\n\n");
      exit(1);
    }
    
    printf("What is the country boundary file?                        \n\n");
    scanf("%s", filcnt);
    
    printf("Do you want just first landfall or all land-ocean transitions, '0' for all and '1' for first landfall only.\n\n");
    scanf("%d", &ilndf);
    
    fcnt = fopen(filcnt, "r");
    if(!fcnt){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filcnt);
       exit(1);
    } 
    
    sscanf(argv[2], "%d", &ilng);    
    sscanf(argv[3], "%d", &ilat);
    sscanf(argv[4], "%d", &iint);
    
    sscanf(argv[6], "%d", &idet);
    if(idet < 0 || idet > 1){
       printf("****ERROR****, incorrect specifier for DET option.\n\n");
       exit(1);    
    }
    
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
    
    dir = opendir(argv[1]);
    if(dir == NULL){
       printf("****ERROR****, cannot open directory %s\n\n", argv[1]);
       exit(1);
    }
    entry = readdir( dir );
    
    fout = fopen("ann.dat", "w");
    
    while(entry != NULL){
    
       if(!strstr(entry->d_name, ".")){
       
           if(!idet)
              sprintf(dirnam, "%s/%s/%s", argv[1], entry->d_name, argv[5]);
	   else
	      sprintf(dirnam, "%s/%s/%s/DET", argv[1], entry->d_name, argv[5]);

	  
          printf("%s\n", dirnam); 
	  
	  sscanf(entry->d_name, "%ld", &fcst);
	 
       
          dirt = opendir(dirnam);
       
          entryt = readdir( dirt );
	  
       
          while(entryt != NULL){
             if(! strncmp( entryt->d_name, "trmatch", 7) && ! strstr(entryt->d_name, "mean")){
	     	     
	        sprintf(filin, "%s/%s", dirnam, entryt->d_name);
		
		printf("%s\n", filin);

	        fin = fopen(filin, "r");
                if(!fin){
                   printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
                   exit(1);
                }

                tracks = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
		
		if(! tracks->trpt->time){
		   printf("****ERROR****, real time not attached to tracks for file %s.\n\n", filin);
		   exit(1);		
		}
        
                fclose(fin); fin = NULL;
       
                convert_track_add(tracks, trnum, upos, ipos);
		
/* analysis landfall */		
    
                nlnd_an = 0;
                lndf_an = find_lndf(tracks, vcnt, ncnt, &nlnd_an, ilndf, ilng, ilat, iint);
		
		if(!nlnd_an) {entryt = readdir( dirt ); continue;}
		
		if(!ilndf){
		   for(i=0; i < nlnd_an; i++) fprintf(fout, "%f %f %ld\n", (lndf_an + i)->lng, (lndf_an + i)->lat, (lndf_an + i)->time);
		}
		else {
		   fprintf(fout, "%f %f %ld\n", lndf_an->lng, lndf_an->lat, lndf_an->time);		
		}
		
		
/* control landfall */

                if( (tracks + 1)->trid == 1 && !idet) {
		
		   nlnd = 0;
		   lndf = find_lndf(tracks + 1, vcnt, ncnt, &nlnd, ilndf, ilng, ilat, iint);
		   
/* find lead time */

                   lead_time_error(lndf_an, lndf, tcntrl, tabcntrl, dcntrl, scntrl, sabcntrl, ncntrl, nlnd, nlnd_an, fcst);

		
		   free(lndf);
		
		}
		
/* loop over all members */

                for(i=1; i < trnum; i++){
		    trr = tracks + i;
		    
		    if(trr->trid == 1 && !idet) continue;
		    
		    nlnd = 0;
		    lndf = find_lndf(trr, vcnt, ncnt, &nlnd, ilndf, ilng, ilat, iint);

                    lead_time_error(lndf_an, lndf, tmem, tabmem, dmem, smem, sabmem, nmem, nlnd, nlnd_an, fcst);

		
		    free(lndf);		     
		
		
		}

/* check for mean track file */
	  
	        sprintf(filin, "%s_mean", filin);
		printf("%s\n", filin);
	        fin = fopen(filin, "r");
                if(fin){		
                   trkmn = read_tracks(fin, &trnummn, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
                   if(! trkmn->trpt->time){
		      printf("****ERROR****, real time not attached to tracks for file %s.\n\n", filin);
		      exit(1);		
		   }
        
                   fclose(fin); fin = NULL;
       
                   convert_track_add(trkmn, trnummn, upos, ipos);
		   		 
		   nlnd_m = 0;
		   
                   lndf_m = find_lndf(trkmn, vcnt, ncnt, &nlnd_m, ilndf, ilng, ilat, iint);
		   
/* find lead time */

                   lead_time_error(lndf_an, lndf_m, tmean, tabmean, dmean, smean, sabmean, nmean, nlnd_m, nlnd_an, fcst);

		       		   
/* calculate spread */

                   for(i=1; i < trnum; i++){
		   
		       nlnd = 0;
		       lndf = find_lndf(trr, vcnt, ncnt, &nlnd, ilndf, ilng, ilat, iint);

                       lead_time_error(lndf_m, lndf, tspr, tabspr, dspr, sspr, sabspr, nspr, nlnd, nlnd_m, fcst);

		
		       free(lndf);		   
		   
		   }	
		   
		   free(lndf_m);	   	
	           	
		}
		
		
		free(lndf_an);
		
                free_tracks(tracks, trnum);
       
             }
             entryt = readdir( dirt );
          }
          
          closedir(dirt);
       
       }
       
       entry = readdir( dir );
    }
    
    fclose(fout);

    for(i=0; i < NFSTT; i++){
        if(nmem[i]) {
	   tmem[i] /= nmem[i];
	   tabmem[i] /= nmem[i];
	   dmem[i] /= nmem[i];
	   smem[i] /= nmem[i];
	   sabmem[i] /= nmem[i];
	}
	if(nmean[i]){
	   tmean[i] /= nmean[i];
	   tabmean[i] /= nmean[i];
	   dmean[i] /= nmean[i];
	   smean[i] /= nmean[i];
	   sabmean[i] /= nmean[i];
	}
	if(ncntrl[i]){
	   tcntrl[i] /= ncntrl[i];
	   tabcntrl[i] /= ncntrl[i];
	   dcntrl[i] /= ncntrl[i];
	   scntrl[i] /= ncntrl[i];
	   sabcntrl[i] /= ncntrl[i];	   
	}
	if(nspr[i]){
	   tspr[i] /= nspr[i];
	   tabspr[i] /= nspr[i];
	   dspr[i] /= nspr[i];
	   sspr[i] /= nspr[i];
	   sabspr[i] /= nspr[i];	
	}
    }
 
    if(!idet){
       write_stats(tmem, tabmem, dmem, smem, sabmem, nmem, 0);   
       write_stats(tcntrl, tabcntrl, dcntrl, scntrl, sabcntrl, ncntrl, 1);
       write_stats(tmean, tabmean, dmean, smean, sabmean, nmean, 2);
       write_stats(tspr, tabspr, dspr, sspr, sabspr, nspr, 3);
    }
    else{
       write_stats(tmem, tabmem, dmem, smem, sabmem, nmem, 4);
    }
    
    closedir(dir);
    
    free(vcnt);

    return 0;
}


LNDF *find_lndf(struct tot_tr *track, VEC *vcnt, int ncnt, int *nlnd, int ilndf, int ilng, int ilat, int iint)
{

    int i=0, j=0, k=0;
    int ic=0, iout=0, inr=0;
    int iflnd=0;
    
    float xx=0.0, yy=0.0;
    
    double dab=0.0, dp1=0.0, dp2=0.0, norm=0.0;
    double dvv=0.0, dlr=0.0;
    
    struct fet_pt_tr *fp=NULL, *fp1=NULL;
    
    struct feature_pts fpt={{0}, {0}, {0}, {0}, {0.0, 0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, {0.0, 0.0}, 0, 0, 0, 0};
    
    VEC *v1=NULL, *v2=NULL;
    VEC a, b, ab, vf;
    VEC vv1, vv2, vv;
    
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
	
	dab = dotp(&a, &b);
	    
	crosp(&a, &b, &ab);
	    
        for(j=1; j < ncnt; j++){
	    
	    v1 = vcnt + j - 1;
            v2 = vcnt + j;
		
	    vv1.x = v1->x;
	    vv1.y = v1->y;
	    vv1.z = v1->z;
	    vv2.x = v2->x;
	    vv2.y = v2->y;
	    vv2.z = v2->z;		
	    
	    dp1 = -dotp(&vv1, &ab);
	    dp2 = dotp(&vv2, &ab);
		
            mulv(&vv1, &vv1, dp2);
	    mulv(&vv2, &vv2, dp1);
		
	    addv(&vv1, &vv2, &vf);
		
	    norm = sqrt(dotp(&vf, &vf));
		
	    normv(&vf, norm);
		
	    dvv = dotp(v1, v2);
	    dp1 = dotp(v1, &vf);
	    dp2 = dotp(v2, &vf);
		
            if(dvv <= ((dp1 < dp2) ? dp1: dp2)){
	    
	       dp1 = dotp(&a, &vf);
	       dp2 = dotp(&b, &vf);
	       
	       if(dab <= ((dp1 < dp2) ? dp1: dp2)){
		   
                  iout = 0;
		   
	          if(!ic) {
		     ic = 1; 
		     iout = 1;
		     lndf = (LNDF *)calloc(1, sizeof(LNDF));
		     mem_er((lndf == NULL) ? 0 : 1, sizeof(LNDF));
		     *nlnd = 1;
	          }
	          else {
		     inr = 0;
		     for(k=0; k < *nlnd; k++){
                         if(labs(timesep((lndf + *nlnd - 1)->time, fp->time)) <= 6 || 
			    labs(timesep((lndf + *nlnd - 1)->time, fp1->time)) <= 6   ) inr = 1;
		     }
		     if(!inr) {
		        iout = 1;
		        *nlnd += 1;
		        lndf = (LNDF *)realloc_n(lndf, *nlnd * sizeof(LNDF));
		        mem_er((lndf == NULL) ? 0 : 1, *nlnd * sizeof(LNDF));
		     }
	          }	  
		
	          if(iout){
		  
		     if(ilndf) iflnd = 1;
		     
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
		     (lndf + *nlnd - 1)->tstep = i;
	             (lndf + *nlnd - 1)->trid = track->trid; 
		      
	          }
	       
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

void lead_time_error(LNDF *lndf_an, LNDF *lndf, double *tlnd, double *tablnd, double *dlnd, double *slnd, double *sablnd, int *npp, int nlnd, int nlnd_an, long int fcst)
{

    int i=0, j=0;
    int ii=0;
    
    int nbin=0;
    int nhr=0, nhm=0, nmin=0;
    int imat=0;
    
    int *ian=NULL;
    
    double dist=0.0, dstr=0.0;
    
    LNDF *l1=NULL, *l2=NULL;
    
    ian = (int *)calloc(nlnd_an, sizeof(int));
    mem_er((ian == NULL) ? 0 : 1, nlnd_an*sizeof(int));
    for(i=0; i < nlnd_an; i++) *(ian + i) = 1;

    for(i=0; i < nlnd; i++){
    
        l1 = lndf + i;

        nhr = timesep(l1->time, fcst);
	if(nhr < 0){
           printf("****ERROR****, negative forecast lead time.\n\n");
	   exit(1);
        }
        nbin = nhr / TSTEP;
		      
	nmin = ILARGE;
		      
	imat = 0;
	for(j=0; j < nlnd_an; j++){
	   nhr = timesep((lndf_an + j)->time, fcst);
	   if(nhr < 0) continue;
	   imat = 1; 
           nhr = timesep(l1->time, (lndf_an + j)->time);
	   if(abs(nhr) < nmin) {nmin = abs(nhr); nhm = nhr; ii = j;}

        }
	
	if(imat && *(ian + ii)){
	   *(ian + ii) = 0;
	   l2 = lndf_an + ii;
	   ++(npp[nbin]);
           tlnd[nbin] += (float) nhm;
           tablnd[nbin] += fabs((float) nhm);
		      
           dist = (l1->vp).x * (l2->vp).x + (l1->vp).y * (l2->vp).y + (l1->vp).z * (l2->vp).z;
	   dlnd[nbin] += acos(dist) / FP_PI;
	   
	   dstr = l1->str - l2->str;
	   slnd[nbin] += dstr;
	   sablnd[nbin] += fabs(dstr);
	}
    }
    
    free(ian);

    return;
    
}

void write_stats(double *tlnd, double *tablnd, double *dlnd, double *slnd, double *sablnd, int *npp, int ityp)
{

   int i=0;
   
   if(!ityp) printf("Ensemble Member Tracks\n\n");
   else if(ityp == 1) printf("Control Tracks\n\n");
   else if(ityp == 2) printf("Mean Tracks\n\n");
   else if(ityp == 3) printf("Ensemble Spread\n\n");
   else if(ityp == 4) printf("Deterministic Tracks\n\n");
   else {
      printf("****ERROR****, track type %d not recognised.\n\n", ityp);
      exit(1);   
   }
   
   printf("Lead Time (days)  Point Number  Time Absolute Error (hrs)   Time Bias (hrs)  Distance Error (deg.) Intensity Absolute Error  Intensity Bias\n\n");
   
   for(i=0; i < NFSTT; i++){
       printf("%6.2f %6d %10.4f %10.4f %10.4f %10.4f %10.4f\n", ((float) i) / 4, *(npp + i), *(tablnd + i), *(tlnd + i), *(dlnd + i), *(sablnd + i), *(slnd + i));   
   }

   return;
   
}
