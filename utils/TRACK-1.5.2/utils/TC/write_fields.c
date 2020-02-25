#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"
#include "vec.h"

#define  SMALL -1.0e+12
#define  LARGE 1.0e+12;
#define  GTOL  1.0e-4

/* program to identify tropical cyclones */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );

int main(void )
{
    int i=0;
    int ii=0;
    int ntr=0;
    int irdim=0;
    int irtrn=0, irnr=0, irnth=0, irnf=0, iptnum=0;
    int gpr=0, ipr=0;
    int trnum=0;
    int nctr=0, ifld=0, np=0;
    int ich=0, itrid=0;

    float alat=0.0, alng=0.0;
    float *sradf=NULL;
    float *slng=NULL, *slat=NULL;

    char filrad[MAXCHR], filin[MAXCHR];
    char line[MAXCHR];
    char trfile[]="rad_track.dat";

    FILE *frad=NULL, *fin=NULL;
    FILE *fout=NULL;

    off_t place1, place2, blklen;

    struct tot_tr *tracks=NULL, *atr=NULL;

    printf("What is the track file to read?\n\n");
    scanf("%s", filin);

    fin = fopen(filin, "r");
    if(!fin){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filin);
       exit(1);
    }

    tracks = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);

    printf("What is the filename for the radial field?\n\n");
    scanf("%s", filrad);
    frad = fopen(filrad, "r");
    if(!frad){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", filrad);
       exit(1);
    }
    fgets(line, MAXCHR, frad);
    sscanf(line, "%d %d %d %d %d", &irtrn, &iptnum, &irnth, &irnr, &irnf);

    printf("There are %d fields associated with the tracks, which field is required?\n\n", irnf);
    scanf("%d", &ifld);

    irdim = irnr * irnth;

    sradf = (float *)calloc(irdim, sizeof(float));
    mem_er((sradf == NULL) ? 0 : 1, irdim * sizeof(float));
    slng = (float *)calloc(irnth, sizeof(float));
    mem_er((slng == NULL) ? 0 : 1, irnth * sizeof(float));
    slat = (float *)calloc(irnr, sizeof(float));
    mem_er((slat == NULL) ? 0 : 1, irnr * sizeof(float));
    fread(slng, irnth*sizeof(float), 1, frad);
    fscanf(frad, "%*c");
    fread(slat, irnr*sizeof(float), 1, frad);
    fscanf(frad, "%*c");

    place1 = ftello(frad);
    fread(sradf, irdim*sizeof(float), 1, frad);
    fscanf(frad, "%*c");
    place2 = ftello(frad);
    blklen = place2 - place1;
    fseeko(frad, place1, SEEK_SET);

    printf("Choose track by position in the file, '0' or by track ID, '1'\n\n");
    scanf("%d", &ich);
    if(!ich){
       printf("Number of tracks is %d, what track is required.\n\n", irtrn);
       scanf("%d", &ntr);
    }
    else{
       printf("What track ID is required.\n\n");
       scanf("%d", &itrid);
       for(i=0; i < trnum; i++){
           atr = tracks + i;
           if(atr->trid == itrid){ntr = i + 1; break;}
       }
    }

    printf("****INFORMATION****, using track number %d\n\n", ntr);

    nctr = 0;

    for(i=0; i < ntr - 1; i++){
        atr = tracks + i;
        nctr += atr->num;
    }


/* open file for writing */

    fout = fopen(trfile, "w");
    if(!fin){
       printf("***ERROR***, unable to open file %s for 'r'\n\n", trfile);
       exit(1);
    }

    atr = tracks + ntr - 1;

    fprintf(fout, "%d %d %d\n", irnth, irnr, atr->num);
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

    for(i=0; i < atr->num; i++){
        fseeko(frad, ((nctr + np) * irnf + ifld - 1) * blklen + place1, SEEK_SET);
        fread(sradf, irdim*sizeof(float), 1, frad);
        fscanf(frad, "%*c");
        fprintf(fout, "FIELD %d\n", i+1);
        fwrite(sradf, irdim * sizeof(float), 1, fout);
        fprintf(fout, "\n");
        ++np;
    }

    fclose(fout);

    return 0;

}


