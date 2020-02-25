#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include "grid.h"
#include "m_values.h"

/* function to interpolate between latitudes for hole filling */

extern GRID *gr;

int missing(float, float, int);


void lat_interp(float *abuf, float dlat, float lat1, float lat2, float mval, int il1, int il2, int ix, int isin, int icmp)
{

    
    int i, j;

    float slp, cint;
    float fl1, fl2;
    float lat; 

    for(i=0; i < ix; i++){

        fl1 = *(abuf + il1 * ix + i);
        fl2 = *(abuf + il2 * ix + i);


        if(missing(fl1, mval, icmp) || missing(fl2, mval, icmp)){
           if(missing(fl1, mval, icmp)){
              printf("Latitude %d, %f contains missing data values\n\n", il1, lat1);
           }
           if(missing(fl2, mval, icmp)){
              printf("Latitude %d, %f contains missing data values\n\n", il2, lat2);
           }
           exit(1);
        }

        slp = (fl1 - fl2) / dlat;
        cint = (fl2 * lat1 - fl1 * lat2) / dlat;
        
        for(j=il1+1; j < il2; j++){

            lat = *(gr->ygrid + j);

            if(isin) lat = sin(lat * FP_PI);

            *(abuf + j * ix + i) = slp * lat + cint;


        }

        

    }



    return;

}
