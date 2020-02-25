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

#define NEPS  21
#define LENSTR 100
#define TIMINT 6
#define MAXGEN 12

#define MANG 10.0
#define FP_PI 0.017453292519943295     /* pi divided by 180 in radians */
#define  FP_PI2  1.57079632679489661923   /* pi divided by 2 in radians   */

#define NSTP 8


/* program to convert EPS track files into kml formated files for use in google earth */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);
void write_kml_header(FILE * , long int , int , char * );
void close_kml_tags(FILE * );
int disp(struct tot_tr * , double );
void free_tracks(struct tot_tr * , int );
void convert_tracks(struct tot_tr * , int );
long int new_time(long int , int );
void write_kml_track(FILE * , struct tot_tr * , int );
void write_kml_cluster(FILE * , struct tot_tr * , struct tot_tr * , int , int , int );
void process_clusters(FILE * , int * , int , int );

int aniso;
int nff, nfld;
int *nfwpos;

int main(int argc, char **argv)
{
   int i, j, k, l;
   int nmch=0;
   int nfil=0, nnom=0;
   int hemi=0;

   int tr_count=0, mtr_count=0;
   int gpr, ipr;
   int gprt, iprt;
   int ngtr=0;
   int ifr=0;

   int ittyp=0, iclstr=0;

   int yr, mn, dy, hr;

   int *iclass=NULL;

   float alat, alng;
   float alatt, alngt;

   double dispmin=MANG*FP_PI;

   long int datei, datef;
   long int stdate=0, ctim=0;

   char lsstr[LENSTR];
   char trfil[LENSTR], trfil2[LENSTR];
   char fkml[]="alltr.kml";
   char ffld[5];
   char date[20];

   FILE *fkmlw=NULL;
   FILE *ftr=NULL;

   struct dirent *entry;

   DIR *dir=NULL;

   struct tot_tr *alltr=NULL, *mtr=NULL, *trr=NULL;
   struct fet_pt_tr *atr;

   if(argc != 4){
      printf("Usage: eps2kml [Date] [Field: vor or mslp] [Hemisphere: NH=0; SH=1]\n\n");
      exit(1);
   } 

   sscanf(argv[1], "%ld", &datei);
   sscanf(argv[2], "%s", ffld);
   sscanf(argv[3], "%d", &hemi);

   sprintf(lsstr, "./%ld_%s", datei, ffld);

   dir = opendir(lsstr);
   if(dir == NULL){
      printf("****ERROR****, cannot open directory %s\n\n", lsstr);
      exit(1);
   }
   entry = readdir( dir );
   while(entry != NULL){
     if( strncmp( entry->d_name, "trmatch", 7 ) == 0 && ! (strstr(entry->d_name, "mean"))) ++nfil;
     else if( strncmp( entry->d_name, "trnomatch", 9 ) == 0) ++nnom;
     entry = readdir( dir );
   }
   closedir(dir);

   printf("%d %d\n", nfil, nnom);

   if(chdir(lsstr)){
     printf("****ERROR****, can't change directory to %s\n", lsstr);
     exit(1);
   }

/* partition clusters according to mean track properties */

   iclass = (int *)calloc(nfil, sizeof(int));
   mem_er((iclass == NULL) ? 0 : 1, nfil*sizeof(int));   

   for(i=0; i < nfil; i++){
      iclstr = i + 1;
      if(i+1 < 10) {
         sprintf(trfil, "trmatch_cntl_tr000%d_mean", i+1);
         sprintf(trfil2, "trmatch_cntl_tr000%d", i+1);
      }
      else if(i+1 < 100) {
         sprintf(trfil, "trmatch_cntl_tr00%d_mean", i+1);
         sprintf(trfil2, "trmatch_cntl_tr00%d", i+1);
      }
      else if(i+1 < 1000) {
         sprintf(trfil, "trmatch_cntl_tr0%d_mean", i+1);
         sprintf(trfil2, "trmatch_cntl_tr0%d", i+1);
      }
      else {
         printf("****ERROR****, not configured for this many clusters.\n\n");
         exit(1);
      }

      ftr = fopen(trfil, "r");
      if(!ftr){
        *(iclass + i) = -1; 
        continue; 
      }

      fclose(ftr);

      ftr = fopen(trfil2, "r");
      if(!ftr){
         printf("****ERROR****, file does not exist.\n\n");
         exit(1);
      }

      alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
      fclose(ftr);

      if(!ifr){gprt = gpr; iprt = ipr; alatt = alat; alngt = alng; ifr=1;}
      else{
         if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
            printf("****ERROR****, track files are incompatable.\n\n");
            exit(1);
         }
      }

      if(alltr->trid){
         printf("****ERROR****, first track file is not control.\n\n");
         exit(1);
      }

/* determine genesis start */

      ngtr = 1;
      ctim = datei;

      while(alltr->trpt->time != ctim) {++ngtr; ctim = new_time(ctim, TIMINT);}

      if(disp(alltr, dispmin) && alltr->num >= NSTP && ngtr <= MAXGEN) *(iclass + i) = 1;
      else *(iclass + i) = 0;

      free_tracks(alltr, tr_count);
   }

   fkmlw = fopen(fkml, "w");

/* write header information */

   write_kml_header(fkmlw, datei, hemi, ffld);

   fprintf(fkmlw, "\t<Folder>\n");
   fprintf(fkmlw, "\t\t<name>Clustered</name>\n");
   fprintf(fkmlw, "\t\t<open>0</open>\n");
   fprintf(fkmlw, "\t\t<visibility>1</visibility>\n"); 
   fprintf(fkmlw, "\t\t<description>Matches against control</description>\n");

/* write matched clusters and cluster means */

   fprintf(fkmlw, "\t\t<Folder>\n");
   fprintf(fkmlw, "\t\t\t<name>Mobile storms with mean tracks</name>\n");
   fprintf(fkmlw, "\t\t\t<open>0</open>\n");
   fprintf(fkmlw, "\t\t\t<visibility>1</visibility>\n");
   fprintf(fkmlw, "\t\t\t<description><![CDATA[<p>Mean track: lifetime &gt; 2days and disp. &gt; 1000km. and Genesis within first 3 days.</p>]]></description>\n");

   process_clusters(fkmlw, iclass, 1, nfil);

   fprintf(fkmlw, "\t\t</Folder>\n");

/* matched storms of limited mobility. */
   fprintf(fkmlw, "\t\t<Folder>\n");
   fprintf(fkmlw, "\t\t\t<name>Less mobile storms with mean tracks.</name>\n");
   fprintf(fkmlw, "\t\t\t<open>0</open>\n");
   fprintf(fkmlw, "\t\t\t<visibility>0</visibility>\n");
   fprintf(fkmlw, "\t\t\t<description>i<![CDATA[<p>Mean track: lifetime &lt; 2days or disp. &lt; 1000km or Genesis &gt; 3 days.</p>]]></description>\n"); 

   process_clusters(fkmlw, iclass, 0, nfil);

   fprintf(fkmlw, "\t\t</Folder>\n");

/* matched storms without mean tracks */

   fprintf(fkmlw, "\t\t<Folder>\n");
   fprintf(fkmlw, "\t\t\t<name>Clusters with too few tracks for mean.</name>\n");
   fprintf(fkmlw, "\t\t\t<open>0</open>\n");
   fprintf(fkmlw, "\t\t\t<visibility>0</visibility>\n");
   fprintf(fkmlw, "\t\t\t<description>Clusters with fewer than 5 members.</description>\n");

   process_clusters(fkmlw, iclass, -1, nfil);

   fprintf(fkmlw, "\t\t</Folder>\n");

   fprintf(fkmlw, "\t</Folder>\n");

/* non matched storms */
   fprintf(fkmlw, "\t<Folder>\n");
   fprintf(fkmlw, "\t\t<name>Non Clustered</name>\n");
   fprintf(fkmlw, "\t\t<open>0</open>\n");
   fprintf(fkmlw, "\t\t<visibility>0</visibility>\n");
   fprintf(fkmlw, "\t\t<description>Non matched storms</description>\n");

   for(i=0; i < nnom; i++){

      fprintf(fkmlw, "\t\t<Folder>\n");
      fprintf(fkmlw, "\t\t<name>Non-matched for Ensemble %d</name>\n", i+1);
      fprintf(fkmlw, "\t\t<open>0</open>\n");
      fprintf(fkmlw, "\t\t<visibility>0</visibility>\n");
      fprintf(fkmlw, "\t\t\t<Style><ListStyle><listItemType>checkHideChildren</listItemType></ListStyle></Style>\n");
      if(i+1 < 10) {
         sprintf(trfil2, "trnomatch_ens000%d", i+1);
      }
      else if(i+1 < 100) {
         sprintf(trfil2, "trnomatch_ens00%d", i+1);
      }
      else if(i+1 < 1000) {
         sprintf(trfil2, "trnomatch_ens0%d", i+1);
      }
      else {
         printf("****ERROR****, not configured for this many files.\n\n");
         exit(1);
      } 

      ftr = fopen(trfil2, "r");
      if(!ftr){
         printf("****ERROR****, file does not exist.\n\n");
         exit(1);
      }

      alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
      fclose(ftr);

      convert_tracks(alltr, tr_count);

      if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
         printf("****ERROR****, track files are incompatable.\n\n");
         exit(1);
      }

      fprintf(fkmlw, "\t\t\t\t<Placemark>\n"); 
      fprintf(fkmlw, "\t\t\t\t\t<styleUrl>#epsln</styleUrl>\n");
      fprintf(fkmlw, "\t\t\t\t\t<MultiGeometry>\n"); 
      for(j=0; j < tr_count; j++){
         trr = alltr + j;
         write_kml_track(fkmlw, trr, 0);
      }

      fprintf(fkmlw, "\t\t\t\t\t</MultiGeometry>\n");
      fprintf(fkmlw, "\t\t\t\t</Placemark>\n"); 
      fprintf(fkmlw, "\t\t</Folder>\n");

      free_tracks(alltr, tr_count);

   }


   fprintf(fkmlw, "\t</Folder>\n");

   close_kml_tags(fkmlw);

   fclose(fkmlw);

   free(iclass);

   return 0;
}


void write_kml_header(FILE *fkmlw, long int date, int hemi, char *ffld)
{
   char line1[]="<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
   char line2[]="<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">";

   fprintf(fkmlw, "%s\n", line1);
   fprintf(fkmlw, "%s\n", line2);
   fprintf(fkmlw, "<Document>\n");
   if(!hemi) fprintf(fkmlw, "\t<name>NH_eps_%ld_%s</name>\n", date, ffld);
   else fprintf(fkmlw, "\t<name>SH_eps_%ld_%s</name>\n", date, ffld);
   fprintf(fkmlw, "\t<open>0</open>\n");
   fprintf(fkmlw, "\t<visibility>1</visibility>\n");
/* line styles */
   fprintf(fkmlw, "\t<Style id=\"epsln\">\n");
   fprintf(fkmlw, "\t\t<LineStyle><color>55ffff77</color><width>2</width></LineStyle>\n"); 
   fprintf(fkmlw, "\t</Style>\n");
   fprintf(fkmlw, "\t<Style id=\"cntlln\">\n");
   fprintf(fkmlw, "\t\t<LineStyle><color>aa0000ff</color><width>3</width></LineStyle>\n");
   fprintf(fkmlw, "\t</Style>\n");
   fprintf(fkmlw, "\t<Style id=\"meanln\">\n");
   fprintf(fkmlw, "\t\t<LineStyle><color>aa00aa00</color><width>3</width></LineStyle>\n");
   fprintf(fkmlw, "\t</Style>\n");
   fprintf(fkmlw, "\t<Style id=\"detln\">\n");
   fprintf(fkmlw, "\t\t<LineStyle><color>aa000000</color><width>3</width></LineStyle>\n");
   fprintf(fkmlw, "\t</Style>\n");
/* end of styles */

   if(!hemi){
      if(!strcmp(ffld, "vor")) {
        fprintf(fkmlw, "\t<description><![CDATA[<p>Northern Hemisphere Vorticity Tracks</p>");
      }
      else if(!strcmp(ffld, "mslp")) {
        fprintf(fkmlw, "\t<description><![CDATA[<p>Northern Hemisphere MSLP Tracks</p>");
      }
      else {
        printf("****ERROR****, field type unknown.\n\n");
        exit(1);
      }
      fprintf(fkmlw, "<p>Blue lines - Ensemble members</p><p>Red lines - Control</p><p>Green lines - Ensemble mean</p><p>Black lines - Determinstic</p>]]></description>\n");

      fprintf(fkmlw, "\t<LookAt>\n");
      fprintf(fkmlw, "\t\t<longitude>-30.0</longitude>\n");
      fprintf(fkmlw, "\t\t<latitude>60.0</latitude>\n");
      fprintf(fkmlw, "\t\t<altitude>0</altitude>\n");
      fprintf(fkmlw, "\t\t<tilt>0</tilt>\n");
      fprintf(fkmlw, "\t\t<range>9044700.0</range>\n");
      fprintf(fkmlw, "\t</LookAt>\n");
   }
   else {
      if(!strcmp(ffld, "vor")) {
        fprintf(fkmlw, "\t<description><![CDATA[<p>Southern Hemisphere Vorticity Tracks</p>");
      }
      else if(!strcmp(ffld, "mslp")) {
        fprintf(fkmlw, "\t<description><![CDATA[<p>Southern Hemisphere MSLP Tracks</p>");
      }
      else {
        printf("****ERROR****, field type unknown.\n\n");
        exit(1);
      }
      fprintf(fkmlw, "<p>Blue lines - Ensemble members</p><p>Red lines - Control</p><p>Green lines - Ensemble mean</p><p>Black lines - Determinstic</p>]]></description>\n");

      fprintf(fkmlw, "\t<LookAt>\n");
      fprintf(fkmlw, "\t\t<longitude>-30.0</longitude>\n");
      fprintf(fkmlw, "\t\t<latitude>-60.0</latitude>\n");
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

void write_kml_cluster(FILE *fkmlw, struct tot_tr *clstr, struct tot_tr *mtr, int trcount, int icld, int type)
{
   int i;
   int ivis=0;

   struct tot_tr *trr=NULL;

   if(type == 1) ivis = 1;
   else ivis = 0; 

   fprintf(fkmlw, "\t\t\t<Folder>\n");
   fprintf(fkmlw, "\t\t\t<visibility>%d</visibility>\n", ivis);
   fprintf(fkmlw, "\t\t\t\t<name>Storm Cluster %d with %d members, control, and mean.</name>\n", icld, trcount-1);
   fprintf(fkmlw, "\t\t\t\t<Style><ListStyle><listItemType>checkHideChildren</listItemType></ListStyle></Style>\n");

/* write eps members */
   fprintf(fkmlw, "\t\t\t\t<Placemark>\n"); 
   fprintf(fkmlw, "\t\t\t\t\t<description> </description>\n");
   fprintf(fkmlw, "\t\t\t\t\t<styleUrl>#epsln</styleUrl>\n");
   fprintf(fkmlw, "\t\t\t\t\t<MultiGeometry>\n");
   if(mtr){
      fprintf(fkmlw, "\t\t\t\t\t<Point><coordinates>%f, %f, %f</coordinates></Point>\n", mtr->trpt->xf, mtr->trpt->yf, mtr->trpt->zf);
   }
   else {
      fprintf(fkmlw, "\t\t\t\t\t<Point><coordinates>%f, %f, %f</coordinates></Point>\n", clstr->trpt->xf, clstr->trpt->yf, clstr->trpt->zf);
   }

   if((clstr + 1)->trid == 1){
      for(i=2; i < trcount; i++){
         trr = clstr + i;
         write_kml_track(fkmlw, trr, 0);
      }
   }
   else {
      for(i=1; i < trcount; i++){
         trr = clstr + i;
         write_kml_track(fkmlw, trr, 0);
      }

   }

   fprintf(fkmlw, "\t\t\t\t\t</MultiGeometry>\n");
   fprintf(fkmlw, "\t\t\t\t</Placemark>\n");

/* write control */
   fprintf(fkmlw, "\t\t\t\t<Placemark>\n");
   fprintf(fkmlw, "\t\t\t\t\t<styleUrl>#cntlln</styleUrl>\n");
   write_kml_track(fkmlw, clstr, 1);
   fprintf(fkmlw, "\t\t\t\t</Placemark>\n");
/* write eps mean */
   if(mtr){
      fprintf(fkmlw, "\t\t\t\t<Placemark>\n");
      fprintf(fkmlw, "\t\t\t\t\t<styleUrl>#meanln</styleUrl>\n");
      write_kml_track(fkmlw, mtr, 1);
      fprintf(fkmlw, "\t\t\t\t</Placemark>\n");
   }
/* write deterministic */
   if((clstr + 1)->trid == 1){
      fprintf(fkmlw, "\t\t\t\t<Placemark>\n");
      fprintf(fkmlw, "\t\t\t\t\t<styleUrl>#detln</styleUrl>\n");
      write_kml_track(fkmlw, clstr+1, 1);
      fprintf(fkmlw, "\t\t\t\t</Placemark>\n");
   }
   fprintf(fkmlw, "\t\t\t</Folder>\n");

   return;
}

void process_clusters(FILE *fkmlw, int *iclass, int type, int nfil)
{
   int i;
   int ifr=0;
   int tr_count=0, mtr_count=0;
   int gpr=0, ipr=0;
   int gprt=0, iprt=0;
   int iclstr=0;

   float alat=0.0, alng=0.0;
   float alatt=0.0, alngt=0.0;

   char trfil[LENSTR];

   FILE *ftr=NULL;

   struct tot_tr *alltr=NULL, *mtr=NULL;


   for(i=0; i < nfil; i++){
      if(*(iclass + i) == type){
         iclstr = i + 1;
         if(i+1 < 10) sprintf(trfil, "trmatch_cntl_tr000%d", i+1);
         else if(i+1 < 100) sprintf(trfil, "trmatch_cntl_tr00%d", i+1);
         else if(i+1 < 1000) sprintf(trfil, "trmatch_cntl_tr0%d", i+1);

         ftr = fopen(trfil, "r");
         if(!ftr){
            printf("****ERROR****, can't open file %s for 'r'\n\n", trfil);
            exit(1);
         }
         alltr = read_tracks(ftr, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
         fclose(ftr);

         convert_tracks(alltr, tr_count);

         if(!ifr){gprt = gpr; iprt = ipr; alatt = alat; alngt = alng; ifr=1;}
         else{
            if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
               printf("****ERROR****, track files are incompatable.\n\n");
               exit(1);
            }
         }

/* read corresponding mean track file */

         strcat(trfil, "_mean");

         ftr = NULL;
         mtr = NULL;
         ftr = fopen(trfil, "r");
         if(ftr){
            mtr = read_tracks(ftr, &mtr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
            fclose(ftr);

            convert_tracks(mtr, 1);

            if(gprt != gpr || iprt != ipr || alatt != alat || alngt != alng) {
               printf("****ERROR****, track files are incompatable.\n\n");
               exit(1);
            }

         }

         write_kml_cluster(fkmlw, alltr, mtr, tr_count, iclstr, type);

         free_tracks(alltr, tr_count);
         if(mtr) free_tracks(mtr, mtr_count);
      }

   }
           
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
