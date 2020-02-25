#include <Stdio.h>
#include <stdlib.h>
#include "st_fo.h"
#include "st_obj.h"
#include "m_values.h"
#include "grid.h"
#include "zones.h"

extern int tf;
extern int x1u, y1u;

extern float period;

extern GRID *gr;

int lngreg(float , float , float );

float dmaxx(ZONE *zone, struct feature_pts *f1, struct feature_pts *f2, float dmax)
{

    static int fmsg=0;

    int i;
    int inb=0;
    float d1=0., d2=0.;
    float lat1=0., lat2=0.;
    float lng1=0., lng2=0.;

    REG *zz=NULL;

    if(!zone){
       printf("****ERROR****, no zonal displacement data to access.\n\n");
       exit(1);

    }

    if(!f1 && !f2){

       if(!fmsg){

          printf("****WARNING****, possible error, both feature points are  \r\n"
                 "                 phantom for upper bound dmax calculation.\n\n");

          fmsg = 1;

       }


       return dmax;

    }

    if(tf == 3) {

       if(f1) {
           lat1 = *(gr->ygrid + ((f1->y).ixy) + y1u - 2);
           lng1 = *(gr->xgrid + ((f1->x).ixy) + x1u - 2);
       }
       if(f2) {
           lat2 = *(gr->ygrid + ((f2->y).ixy) + y1u - 2);
           lng2 = *(gr->xgrid + ((f2->x).ixy) + x1u - 2);
       }
    } 

    else {

       if(f1) {lat1 = (f1->y).xy; lng1 = (f1->x).xy;}
       if(f2) {lat2 = (f2->y).xy; lng2 = (f2->x).xy;}

    }

    if(f1){

       inb = 0;

       for(i=0; i< zone->nz; i++){
           zz = zone->zlat + i;
           if((zz->y2 - lat1) * (lat1 - zz->y1) >= 0.0 && lngreg(zz->x1, zz->x2, lng1)) {d1 = zz->zdmax; inb = 1; break;}  
           
/*              (zz->x2 - lng1) * (lng1 - zz->x1) >= 0.0) {d1 = zz->zdmax; inb = 1; break;} */
       }

       if(!inb){

          printf("****ERROR*****, for zonal regional upper bound displacements \r\n"
                 "                zones must cover whole region.               \r\n"
                 "                Latitude %f or longitude %f not in any       \r\n"
                 "                specified zone.\n\n", lat1, lng1);

          exit(1); 

       }


    }

    else d1 = dmax;


    if(f2){

       inb = 0;

       for(i=0; i< zone->nz; i++){
           zz = zone->zlat + i;
           if((zz->y2 - lat2) * (lat2 - zz->y1) >= 0.0 && lngreg(zz->x1, zz->x2, lng2)) {d2 = zz->zdmax; inb = 1; break;}
           
/*               (zz->x2 - lng2) * (lng2 - zz->x1) >= 0.0) {d2 = zz->zdmax; inb = 1; break;} */
       }

       if(!inb){

          printf("****ERROR*****, for zonal regional upper bound displacements \r\n"
                 "                zones must cover whole region.               \r\n"
                 "                Latitude %f or longitude %f not in any       \r\n"
                 "                specified zone.\n\n", lat2, lng2);

          exit(1); 

       }

    }

    else d2 = dmax;


    return 0.5 * (d1 + d2);

}

int lngreg(float z1, float z2, float lng)
{
    int inreg=0;
    
    if(z1 > z2){
      if((z2 - lng) * lng >= 0.0 ||
         (period - lng) * (lng - z1) >= 0.0) inreg = 1;
    
    }
    else {
      if((z2 - lng) * (lng - z1) >= 0.0) inreg = 1;
    }
    
    return inreg;

}
