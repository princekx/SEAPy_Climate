#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include "splice.h"
#include "mem_er.h"

#define LENSTR 100
#define TIMINT 6
#define MAXGEN 12

#define MANG 10.0
#define FP_PI 0.017453292519943295     /* pi divided by 180 in radians */
#define  FP_PI2  1.57079632679489661923   /* pi divided by 2 in radians   */

#define NSTP 8


/* program to convert normal track files into kml formated files for use in google earth */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);
void write_kml_header(FILE * , int );
void close_kml_tags(FILE * );
int disp(struct tot_tr * , double );
void free_tracks(struct tot_tr * , int );
void convert_tracks(struct tot_tr * , int );
long int new_time(long int , int );
void write_kml_track(FILE * , struct tot_tr * , int );
void write_kml_tracks(FILE * , struct tot_tr * , int );
void process_tracks(FILE * , char * );

int aniso;
int nff, nfld;
int *nfwpos;

int main(int argc, char **argv)
{

   int ior=0;

   char trfil[LENSTR];
   char fkml[]="alltr.kml";

   FILE *fkmlw=NULL;

   if(argc != 3){
      printf("Usage: tr2kml [track file] [orientation, -1 (SH), 0 (tropics), 1 (NH)\n\n");
      exit(1);
   } 

   sscanf(argv[1], "%s", trfil);
   sscanf(argv[2], "%d", &ior);

   fkmlw = fopen(fkml, "w");

/* write header information */

   write_kml_header(fkmlw, ior);

   fprintf(fkmlw, "\t<Folder>\n");
   fprintf(fkmlw, "\t\t<name>Cyclone Tracks</name>\n");
   fprintf(fkmlw, "\t\t<open>0</open>\n");
   fprintf(fkmlw, "\t\t<visibility>1</visibility>\n"); 

/* write matched clusters and cluster means */

   process_tracks(fkmlw, trfil);

   fprintf(fkmlw, "\t</Folder>\n");

   close_kml_tags(fkmlw);

   fclose(fkmlw);


   return 0;
}


void write_kml_header(FILE *fkmlw, int ior)
{
   char line1[]="<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
   char line2[]="<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">";

   fprintf(fkmlw, "%s\n", line1);
   fprintf(fkmlw, "%s\n", line2);
   fprintf(fkmlw, "<Document>\n");
   fprintf(fkmlw, "\t<open>0</open>\n");
   fprintf(fkmlw, "\t<visibility>1</visibility>\n");
/* line styles */
   fprintf(fkmlw, "\t<Style id=\"cntlln\">\n");
   fprintf(fkmlw, "\t\t<LineStyle><color>aa0000ff</color><width>3</width></LineStyle></Style>\n");
/* end of styles */

   if(ior == -1){
      fprintf(fkmlw, "\t<LookAt>\n");
      fprintf(fkmlw, "\t\t<longitude>-30.0</longitude>\n");
      fprintf(fkmlw, "\t\t<latitude>-60.0</latitude>\n");
      fprintf(fkmlw, "\t\t<altitude>0</altitude>\n");
      fprintf(fkmlw, "\t\t<tilt>0</tilt>\n");
      fprintf(fkmlw, "\t\t<range>9044700.0</range>\n");
      fprintf(fkmlw, "\t</LookAt>\n");
   }
   else if(ior == 1) {
      fprintf(fkmlw, "\t<LookAt>\n");
      fprintf(fkmlw, "\t\t<longitude>-30.0</longitude>\n");
      fprintf(fkmlw, "\t\t<latitude>60.0</latitude>\n");
      fprintf(fkmlw, "\t\t<altitude>0</altitude>\n");
      fprintf(fkmlw, "\t\t<tilt>0</tilt>\n");
      fprintf(fkmlw, "\t\t<range>9044700.0</range>\n");
      fprintf(fkmlw, "\t</LookAt>\n");
   }
   else {
      fprintf(fkmlw, "\t<LookAt>\n");
      fprintf(fkmlw, "\t\t<longitude>-30.0</longitude>\n");
      fprintf(fkmlw, "\t\t<latitude>0.0</latitude>\n");
      fprintf(fkmlw, "\t\t<altitude>0</altitude>\n");
      fprintf(fkmlw, "\t\t<tilt>0</tilt>\n");
      fprintf(fkmlw, "\t\t<range>9044700.0</range>\n");
      fprintf(fkmlw, "\t</LookAt>\n");
   }
   

   return;
}

void close_kml_tags(FILE *fkmlw)
{

   fprintf(fkmlw, "</Document>\n");
   fprintf(fkmlw, "</kml>\n");

   return;
}

void write_kml_track(FILE *fkmlw, struct tot_tr *tr, int type)
{
   int i;

   float alt=0.0;

   if(type) alt = 0.0;
   else alt = 0.0;

   struct fet_pt_tr *atr=NULL; 

   fprintf(fkmlw, "\t\t\t\t\t<LineString><tessellate>1</tessellate><altitudeMode>clampToGround</altitudeMode><coordinates>\n"); 

   for(i=0; i < tr->num; i++){
      atr = tr->trpt + i;
      fprintf(fkmlw, "\t\t\t\t\t%f,%f,%f\n", atr->xf, atr->yf, alt);
   }
   fprintf(fkmlw, "\n");
   fprintf(fkmlw, "\t\t\t\t\t</coordinates></LineString>\n");

   return;
}

void write_kml_tracks(FILE *fkmlw, struct tot_tr *clstr, int trcount)
{
   int i;

   struct tot_tr *trr=NULL;

/* write eps members */
   fprintf(fkmlw, "\t\t\t\t<Placemark>\n"); 
   fprintf(fkmlw, "\t\t\t\t\t<description> </description>\n");
   fprintf(fkmlw, "\t\t\t\t\t<styleUrl>#cntlln</styleUrl>\n");
   fprintf(fkmlw, "\t\t\t\t\t<MultiGeometry>\n");

   for(i=0; i < trcount; i++){
      trr = clstr + i;
      write_kml_track(fkmlw, trr, 0);
   }

   fprintf(fkmlw, "\t\t\t\t\t</MultiGeometry>\n");
   fprintf(fkmlw, "\t\t\t\t</Placemark>\n");

   return;
}

void process_tracks(FILE *fkmlw, char *trfil)
{
   int tr_count=0;
   int gpr=0, ipr=0;

   float alat=0.0, alng=0.0;

   FILE *ftr=NULL;

   struct tot_tr *alltr=NULL;

   ftr = fopen(trfil, "r");
   if(!ftr){
      printf("****ERROR****, can't open file %s for 'r'\n\n", trfil);
      exit(1);
   }
   alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
   fclose(ftr);

   convert_tracks(alltr, tr_count);

   write_kml_tracks(fkmlw, alltr, tr_count);

   free_tracks(alltr, tr_count);
           
   return;

}

int disp(struct tot_tr *trr, double dispmin)
{
   int ids=0;

   double dis;
   double thet1, phi1, thet2, phi2;

   struct fet_pt_tr *at1=NULL, *at2=NULL;

   at1 = trr->trpt;
   at2 = trr->trpt + trr->num - 1;

   phi1 = at1->xf * FP_PI;
   phi2 = at2->xf * FP_PI;
   thet1 = FP_PI2 - at1->yf * FP_PI;
   thet2 = FP_PI2 - at2->yf * FP_PI;

   if(thet1 < 0.) thet1 = 0.;
   if(thet2 < 0.) thet2 = 0.;

   dis = sin(thet1) * sin(thet2) * cos(phi2 - phi1) + cos(thet1) * cos(thet2);

   if(fabs(dis) > 1.) dis = (dis < 0.) ? -1.0 : 1.0;

   dis = acos(dis);

   if(dis >= dispmin) ids = 1;

   return ids;
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
           free(atr->add_fld);
       }
       free(trr->trpt);
    }

    free(alltr);

    return;

}

void convert_tracks(struct tot_tr *alltr, int trcount)
{
    int i, j;

    struct tot_tr *trr=NULL;
    struct fet_pt_tr *atr=NULL;

    for(i=0; i < trcount; i++){
       trr = alltr + i;
       for(j=0; j < trr->num; j++){
          atr = (trr->trpt) + j;
          if(atr->xf > 180.0) atr->xf -= 360.0;
       }
    }

    return;
}
