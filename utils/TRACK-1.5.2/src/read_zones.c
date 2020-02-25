#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "files_in.h"
#include "file_handle.h"
#include "mem_er.h"
#include "m_values.h"
#include "zones.h"

#define  MAXCHR  200

/* read zonal upper bound displacements */

extern int tom;

ZONE *read_zones(int parn)

{

    int i;

    char cbuf[MAXCHR];
    char par_file[MAXCHR];
    char num[4];

    ZONE *zz=NULL;
    REG *rr=NULL;

    FILE *fz=NULL;

    sprintf(num, "%d", parn);
    strncpy(par_file, DATZN, MAXCHR);
    strcat(par_file, num);

    fz = open_file(par_file, "r");

    printf("****INFORMATION****, reading regional dmax values from file:\r\n"
           "                     %s\n\n", par_file);

    zz = (ZONE *)malloc_initl(sizeof(ZONE));
    mem_er((zz == NULL) ? 0 : 1, sizeof(ZONE));

    fscanf(fz, "%d\n", &(zz->nz));

    zz->zlat = (REG *)calloc(zz->nz, sizeof(REG));
    mem_er((zz->zlat == NULL) ? 0 : 1, zz->nz * sizeof(REG));

    zz->zonemax = 0.0;

    for(i=0; i<zz->nz; i++) {
        fgets(cbuf, MAXCHR, fz);
        rr = zz->zlat + i;
        sscanf(cbuf, "%f %f %f %f %f", &(rr->x1), &(rr->x2), &(rr->y1), &(rr->y2), &(rr->zdmax));
        if(tom == 'g') rr->zdmax *= FP_PI;
        if(rr->zdmax > zz->zonemax) zz->zonemax = rr->zdmax;
    }

    close_file(fz, par_file);
    return zz;

}
