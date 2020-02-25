#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "splice.h"
#include "mem_er.h"
#include "vec.h"

#define  FP_PI   0.017453292519943295   
#define  FP_PI2  1.570796326794896619

#define  LARGE_DIST  1.0e+12
#define  NTIM        50000
#define  MCHRL       100
#define  TOLLINE     1.0e-4
#define  DELB        1.0

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
long int new_time(long int , int );
void sincos(double , double * , double * );
int point_in_polygon_sphere(VEC * , VEC * , VEC * , int , int );
int point_in_polygon(float * , float * , float , float , float , float , int );
int orientation(VEC * , int );

int noheader=0;

int main(int argc, char *argv[])

{

    int i,j;
    int gpr, ipr, trnum;
    int inn=0;
    int negate=0;
    int itim=0;
    int tstep=6;
    int nbnd=0;
    int ntr=0;
    int igen=0;
    int ior=0;
    int ibtyp=0;
    
    int nnn=0;

    long int ntime[NTIM], ttt;

    float alng=0.0, alat=0.0;
    float clat=0.0, clng=0.0;

    float *blng=NULL, *blat=NULL;
    
    float xmn=100000.0, xmx=-100000.0, ymn=100000.0, ymx=-100000.0;

    double s1, c1, s2, c2;
    double tlng=0.0, tlat=0.0;

    VEC *vbnd=NULL;
    VEC cptt;
    VEC ptt;
    
    char line[MCHRL];
    char newf[MCHRL];
    char com[]="count [trackfile] [boundary] [lat] [lng] [Genesis (0)/Lysis (1)/Passing(2)/All Times(3)] [Negate (1)] [Sphere/Cart boundary (0)] [Start Time, YYYYMMDDHH] [tstep]";

    struct tot_tr *all_tr=NULL, *atr=NULL;
    struct fet_pt_tr *pt=NULL;

    FILE *fin=NULL;
    FILE *fout=NULL;
    FILE *fbnd=NULL;


    if(argc < 8){

       printf("%s \n", com);
       exit(1);


    }
    
    
    sscanf(argv[3], "%f", &clat);
    sscanf(argv[4], "%f", &clng);
    sscanf(argv[5], "%d", &igen);
    sscanf(argv[6], "%d", &negate);
    sscanf(argv[7], "%d", &ibtyp);


    fin = fopen(argv[1], "r");
    if(!fin){

       printf("No such track file exists, %s\n", argv[1]);
       exit(1);

    }

    all_tr = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);

    fbnd = fopen(argv[2], "r");
    if(!fbnd){

       printf("No such boundary file exists, %s\n", argv[2]);
       exit(1);

    }

    fgets(line, MCHRL, fbnd);
    sscanf(line, "%d", &nbnd);
    blng = (float *)calloc(nbnd, sizeof(float));
    mem_er((blng == NULL) ? 0 : 1, nbnd * sizeof(float));
    blat = (float *)calloc(nbnd, sizeof(float));
    mem_er((blat == NULL) ? 0 : 1, nbnd * sizeof(float));

    vbnd = (VEC *)calloc(nbnd, sizeof(VEC));
    mem_er((vbnd == NULL) ? 0 : 1, nbnd*sizeof(VEC)); 
    
    for(i=0; i < nbnd; i++){
       fgets(line, MCHRL, fbnd);
       sscanf(line, "%f %f", blng+i, blat+i); 
       
       if(*(blng+i) < xmn) xmn = *(blng+i);
       if(*(blng+i) > xmx) xmx = *(blng+i);
       
       if(*(blat+i) < ymn) ymn = *(blat+i);
       if(*(blat+i) > ymx) ymx = *(blat+i);      

       tlat = FP_PI2 - *(blat + i) * FP_PI;
       tlng = *(blng + i) * FP_PI;

       sincos(tlng, &s1, &c1);
       sincos(tlat, &s2, &c2);

       (vbnd + i)->x = s2 * c1;
       (vbnd + i)->y = s2 * s1;
       (vbnd + i)->z = c2;
       
    } 
    
    printf("%f %f %f %f\n", xmn, xmx, ymn, ymx);
       
    xmn -= DELB;
    if(xmn < 0.0) xmn = 0.0;
    xmx += DELB;
    if(xmx > 360.0) xmx = 360.0;
       
    ymn -= DELB;
    if(ymn < -90.0) ymn = -90.0;
    ymx += DELB;
    if(ymx > 90.0) ymx = 90.0;      
       
    printf("%f %f %f %f\n", xmn, xmx, ymn, ymx);  
    
    if(*blng != *(blng + nbnd - 1) || *blat != *(blat + nbnd - 1)){
       printf("****ERROR****, not a closed boundary.\n\n");
       exit(1);
    } 

    fclose(fbnd);
    
/* check polygon orientation */

    if(ibtyp) {
    
       ior = orientation(vbnd, nbnd);
    
/* convert chosen interior point to Cartesians */

       tlat = FP_PI2 - clat * FP_PI;
       tlng = clng * FP_PI;

       sincos(tlng, &s1, &c1);
       sincos(tlat, &s2, &c2);

       cptt.x = s2 * c1;
       cptt.y = s2 * s1;
       cptt.z = c2;    
    
    }

    if(igen == 3){

       if(argc >= 9){

          sscanf(argv[8], "%ld", &ttt);
          printf("%ld\n", ttt);
          ntime[0] = ttt;
          sscanf(argv[9], "%d", &tstep);

          for(i=1; i < NTIM; i++)ntime[i] = new_time(ntime[i-1], tstep);
       }

       else {

          printf("****ERROR****, to use actual times require a starting time\n\n");
          printf("%s \n", com);
          exit(1);

       }

    }


    if(igen == 3){

       for(i=0; i < trnum; i++){

          atr = all_tr + i;
          atr->time = ntime[atr->trpt->fr_id - 1];
          itim = atr->trpt->fr_id - 1;

          for(j=0; j< atr->num; j++){

              pt = atr->trpt + j;
              pt->time = ntime[itim];
              itim += 1 + pt->nfm;
          }

       }

    }



    else if(igen == 2){

       

       for(i=0; i < trnum; i++){

          atr = all_tr + i;

          inn = 0;


          for(j=0; j< atr->num; j++){

              pt = atr->trpt + j;
	      
/* pre-filter */
              if((xmx - pt->xf) * (pt->xf - xmn) >= 0.0 && (ymx - pt->yf) * (pt->yf - ymn) >= 0.0) {	      

/* point in polygon */

                 if(!ibtyp) {
	             tlat = FP_PI2 - pt->yf * FP_PI;
                     tlng = pt->xf * FP_PI;

                     sincos(tlng, &s1, &c1);
                     sincos(tlat, &s2, &c2);

                     ptt.x = s2 * c1;
	             ptt.y = s2 * s1;
	             ptt.z = c2;
	             nnn = point_in_polygon_sphere(vbnd, &cptt, &ptt, nbnd, ior);
	         }
	         else nnn = point_in_polygon(blng, blat, clng, clat, pt->xf, pt->yf, nbnd);
	      
	         if(!(nnn % 2)) {inn=1; break;}
	      
	      }

           }

           if(inn && negate) atr->num = 0;

           else if(!inn && !negate) atr->num = 0;


       }

    }

    else {

       for(i=0; i < trnum; i++){

          atr = all_tr + i;

          inn = 0;

          if(!igen) pt = atr->trpt;
          else pt = atr->trpt + atr->num - 1;
	  
/* pre-filter */
          if((xmx - pt->xf) * (pt->xf - xmn) >= 0.0 && (ymx - pt->yf) * (pt->yf - ymn) >= 0.0) {

/* point in polygon */
  
             if(!ibtyp) {
	    
	        tlat = FP_PI2 - pt->yf * FP_PI;
                tlng = pt->xf * FP_PI;
 

                sincos(tlng, &s1, &c1);
                sincos(tlat, &s2, &c2);
	  
	        ptt.x = s2 * c1;
	        ptt.y = s2 * s1;
	        ptt.z = c2;
                nnn = point_in_polygon_sphere(vbnd, &cptt, &ptt, nbnd, ior);
	  
	     }
	     else nnn = point_in_polygon(blng, blat, clng, clat, pt->xf, pt->yf, nbnd);
	  
	     if(!(nnn % 2))inn=1;
	  
	  }

          if(inn && negate) atr->num = 0;

          else if(!inn && !negate) atr->num = 0;	  


       }

    }


    strcpy(newf, argv[1]);
    strcat(newf, ".new");

    fout = fopen(newf, "w");

    meantrd(fout, all_tr, trnum, 's', gpr, ipr, alat, alng);

    fclose(fout);


    return ntr;

}

int point_in_polygon_sphere(VEC *vbnd, VEC *ptin, VEC *ptt, int nbnd, int ior)
{
    int i=0;
    int nint=0;

    
    double norm=0.0;
    double aa1=0.0, bb1=0.0, dotab1=0.0; 
    double aa2=0.0, bb2=0.0, dotab2=0.0;

    VEC cpv, cqv, cpp;
    VEC *p1, *p2;
    
    VEC t1, t2;
    
    crosp(ptin, ptt, &cqv);  
    
    norm = sqrt(dotp(&cqv, &cqv));
    normv(&cqv, norm);
    
    dotab1 = dotp(ptin, ptt);   

    for(i=0; i < nbnd-1; i++){
        p1 = vbnd + i;
	p2 = vbnd + i + 1;
	
        if(!ior) {crosp(p1, p2, &cpv);}
	else {crosp(p2, p1, &cpv);}
	
        norm = sqrt(dotp(&cpv, &cpv));
        normv(&cpv, norm);
	
	crosp(&cpv, &cqv, &cpp);
        norm = sqrt(dotp(&cpp, &cpp));
        normv(&cpp, norm);

        aa2 = dotp(&cpp, p1);
        bb2 = dotp(&cpp, p2);
        dotab2 = dotp(p1, p2);
       
        aa1 = dotp(&cpp, ptt);
        bb1 = dotp(&cpp, ptin);  
	
crosp(p1, &cpp, &t1);
norm = sqrt(dotp(&t1, &t1));
normv(&t1, norm);
crosp(&cpp, p2, &t2);
norm = sqrt(dotp(&t2, &t2));
normv(&t2, norm);
printf("%f\n", dotp(&t1, &t2));
     

        if(dotab2 <= ((aa2 < bb2) ? aa2: bb2) && dotab1 <= ((aa1 < bb1) ? aa1: bb1)) ++nint;
printf("%d %f\n", ior, acos(dotp(&cqv, &cpv))/FP_PI);
printf("%f %f %f\n", aa1,bb1,dotab1);
printf("%f %f %f\n", aa2,bb2,dotab2);
printf("%d %f %f %f\n", i, p1->x, p1->y, p1->z); printf("%d %f %f %f\n", i+1, p2->x, p2->y, p2->z);
printf("%f %f %f\n", cpp.x, cpp.y, cpp.z);
printf("%f %f\n", atan2(p1->y, p1->x)/FP_PI, 90 - acos(p1->z)/FP_PI);
printf("%f %f\n", atan2(p2->y, p2->x)/FP_PI, 90 - acos(p2->z)/FP_PI);
printf("%f %f\n", atan2(ptin->y, ptin->x)/FP_PI, 90 - acos(ptin->z)/FP_PI);
printf("%f %f\n", atan2(ptt->y, ptt->x)/FP_PI, 90 - acos(ptt->z)/FP_PI);
printf("%f %f\n", atan2(cpp.y, cpp.x)/FP_PI, 90 - acos(cpp.z)/FP_PI);
printf("%d\n", nint);
printf("\n\n");   


    }
printf("%d\n", nint);    
exit(1);

    return nint;
}

int point_in_polygon(float *blng, float *blat, float clng, float clat, float ptlng, float ptlat, int nbnd)
{
    int i=0;
    int nint=0;
    
    double pbx1=0.0, pby1=0.0, pbx2=0.0, pby2=0.0;
    double pxm=0.0, pym=0.0;
    double xy1=0.0, xy2=0.0;
    
    double m1, m2, c1, c2;
    
    double dnm=0.0;
    
    for(i=0; i < nbnd-1; i++){
        pbx1 = *(blng + i);
	pby1 = *(blat + i);
	pbx2 = *(blng + i + 1);
	pby2 = *(blat + i + 1);
	
        dnm = (pbx1 - pbx2) * (clat - ptlat) - (clng - ptlng) * (pby1 - pby2);
	
	if(fabs(dnm) > TOLLINE){
	   xy1 = pbx1 * pby2 - pby1 * pbx2;
	   xy2 = clng * ptlat - clat * ptlng;
	   pxm = (xy1 * (clng - ptlng) - xy2 * (pbx1 - pbx2)) / dnm;
	   pym = (xy1 * (clat - ptlat) - xy2 * (pby1 - pby2)) / dnm;

	   if(((pxm - pbx1) * (pbx2 - pxm) >= 0.0 && (pym - pby1) * (pby2 - pym) >= 0.0) &&
	      ((pxm - clng) * (ptlng - pxm) >= 0.0 && (pym - clat) * (ptlat - pym) >= 0.0)  ) {
	      ++nint;
	   }
	   
	} 
    
    }

    return nint;
}

int orientation(VEC *vbnd, int nbnd)
{
   int i=0;
   
   int clk=0;
   
   double dp=0.0, norm=0.0;
   double dsum=0.0;
   
   VEC vt1, vt2, vcr;

   for(i=0; i < nbnd - 2; i++) {
      dp = dotp((vbnd + i + 1), (vbnd + i));
      mulv((vbnd + i + 1), &vt1, dp);
      subv(&vt1, (vbnd + i), &vt1);
      norm = sqrt(dotp(&vt1, &vt1));
      normv(&vt1, norm);
      
      dp = dotp((vbnd + i + 2), (vbnd + i + 1));
      mulv((vbnd + i + 1), &vt2, dp);
      subv((vbnd + i + 2), &vt2, &vt2);
      norm = sqrt(dotp(&vt2, &vt2));
      normv(&vt2, norm); 
      
      crosp(&vt1, &vt2, &vcr)
      norm = sqrt(dotp(&vcr, &vcr));
      normv(&vcr, norm);
      
      dsum += dotp((vbnd + i + 1), &vcr);
      
   }

/* positive sum is anticlockwise */
   
   if(dsum > 0.0) clk = 1;

   return clk;
   
}
